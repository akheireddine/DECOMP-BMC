#include "NuSMV.hh"
#include "../painless/painless-src/utils/Parameters.h"
#include <unordered_set>

NuSMV::NuSMV(const char *inputFile, EnvBMC *envi)
{
	modelFile_.assign(inputFile);
	env = envi;

	/* these are for the core */
	FP_V_E iq_fns[][2] =
		{
			{AddonsCore_Init, AddonsCore_Quit},
		};

	/* Initializes data such as tool name, tool version, email.. */
	NuSMVCore_init_data();

	env_nusmv = NuSMVEnv_create();
	/* Initializes all packages, having the list of init/quit mfunctions */
	NuSMVCore_init(env_nusmv, iq_fns, sizeof(iq_fns) / sizeof(iq_fns[0]));

	/* Adds the command line options of NuSMV */
	NuSMVCore_init_cmd_options(env_nusmv);

	/*************************SOURCE FILE***************************************************/
	read_model();

	// SET rbc2cnf conversion with tseitin algorithm
	Rbc_2CnfAlgorithm rbc_rbc2cnf_algorithm = RBC_TSEITIN_CONVERSION;
	if (envi->convertor == 1)
		rbc_rbc2cnf_algorithm = RBC_SHERIDAN_CONVERSION;

	OptsHandler_ptr opt = OPTS_HANDLER(NuSMVEnv_get_value(env_nusmv, ENV_OPTS_HANDLER));

	set_rbc2cnf_algorithm(opt, rbc_rbc2cnf_algorithm);

	set_verbose_level(opt, 0);
}

NuSMV::~NuSMV()
{
	NuSMVCore_quit(env_nusmv);
	NuSMVEnv_destroy(env_nusmv);
}

void NuSMV::read_model()
{
	ErrorMgr_ptr errmgr;
	OptsHandler_ptr opts;
	StreamMgr_ptr streams;

	opts = OPTS_HANDLER(NuSMVEnv_get_value(env_nusmv, ENV_OPTS_HANDLER));
	errmgr = ERROR_MGR(NuSMVEnv_get_value(env_nusmv, ENV_ERROR_MANAGER));
	streams = STREAM_MGR(NuSMVEnv_get_value(env_nusmv, ENV_STREAM_MANAGER));
	int res;

	set_bmc_mode(opts);

	/* Necessary to have standard behavior in the batch mode */
	ErrorMgr_reset_long_jmp(errmgr);
	CATCH(errmgr)
	{

		{ /* 1: Read the model */
			char *fname;
			if (modelFile_.find(".smv") != std::string::npos)
				fname = strdup(modelFile_.c_str());
			else
				fname = strdup((modelFile_ + ".smv").c_str()); // get_input_file(opts);
			if (NULL == fname)
			{
				StreamMgr_print_error(streams,
									  "Input file is (null). You must set the input file before.\n");
				goto batch_exit_fail;
			}
			res = Parser_read_model(env_nusmv, fname);
			free(fname);
			if (res)
				goto batch_exit_fail;
		}

		{ /* 2: Flatten hierarchy */
			res = CompileFlatten_flatten_smv(env_nusmv, true, false);
			if (res)
				goto batch_exit_fail;
		}

		{ /* 3: Builds the encodings */
			res = Compile_encode_variables(env_nusmv,
										   NULL /*input_order_file_name*/,
										   false /*bdd_enc_enum_only*/);
			if (res)
				goto batch_exit_fail;
		}

		{ /* 4: Builds the flat FSMs */
			res = Compile_create_flat_model(env_nusmv);
			if (res)
				goto batch_exit_fail;
		}
		if (Parameters::isSet("show-var"))
		{
			FILE *outstream = StreamMgr_get_output_stream(streams);
			OStream_ptr ostream = OStream_create(outstream);
			Compile_show_vars(env_nusmv, false, false, true, true, false, true, ostream, true);
			// exit(0);
		}

		/* Start BMC-SAT part */
		if (opt_bmc_mode(opts))
		{
			/* build_boolean_model may have been already called if the output
				boolean model was specified in the argument list. */
			if (Compile_check_if_bool_model_was_built(env_nusmv, NULL, false))
			{
				res = Compile_create_boolean_model(env_nusmv);
				if (res)
					goto batch_exit_fail;
			}

			/* Initializes the bmc package, and commits both the model and the
				determinization layers: */
			res |= Bmc_Pkg_bmc_setup(env_nusmv, false);
			if (res)
				goto batch_exit_fail;

			/* exits */
			goto batch_exit_success;
		}
	} /* end catch */
	FAIL(errmgr)
	{
		StreamMgr_print_error(streams, "\n%s terminated by a signal\n",
							  NuSMVCore_get_tool_name());
		goto batch_exit_fail;
	}

batch_exit_success:
	return;

batch_exit_fail:
	StreamMgr_print_error(streams, "\nAborting batch mode\n");
	ErrorMgr_nusmv_exit(errmgr, 1);
	return;
}

void NuSMV::generate_cnf_formula(int k_val, int prop_value)
{
	set_k(k_val);
	OptsHandler_ptr opts;
	PropDb_ptr prop_db;
	NodeMgr_ptr nodemgr;
	BddEnc_ptr bdd_enc;
	node_ptr bltlspec; /* Its booleanization */
	BeFsm_ptr be_fsm;  /* The corresponding be fsm */
	BeEnc_ptr be_enc;
	Be_Manager_ptr be_mgr;
	Be_Cnf_ptr cnf; /* The CNFed be problem */
	Prop_Rewriter_ptr rewriter = NULL;
	int step_phy_index = -1;
	const StreamMgr_ptr streams =
		STREAM_MGR(NuSMVEnv_get_value(env_nusmv, ENV_STREAM_MANAGER));
	FILE *errstream = StreamMgr_get_error_stream(streams);

	opts = OPTS_HANDLER(NuSMVEnv_get_value(env_nusmv, ENV_OPTS_HANDLER));
	prop_db = PROP_DB(NuSMVEnv_get_value(env_nusmv, ENV_PROP_DB));
	nodemgr = NODE_MGR(NuSMVEnv_get_value(env_nusmv, ENV_NODE_MGR));
	bdd_enc = BDD_ENC(NuSMVEnv_get_value(env_nusmv, ENV_BDD_ENCODER));
	be_fsm = BE_FSM(NuSMVEnv_get_value(env_nusmv, ENV_BE_FSM));

	const Be_CnfAlgorithm cnf_alg = get_rbc2cnf_algorithm(opts);
	Prop_ptr ltlprop = PROP(NULL);
	int rel_loop;
	FILE *dimacsfile;
	char szLoop[16]; /* to keep loopback string */
	char *str_formula;
	int l;
	char *dump_fname_template = strdup((modelFile_ + "_k" + std::to_string(k_)).c_str());
	std::string dump_fname_extension(dump_fname_template);
	dump_fname_extension = dump_fname_extension + ".dimacs";

	Bmc_DumpType dump_type = BMC_DUMP_DIMACS;

	set_bmc_dimacs_filename(opts, dump_fname_template);

	set_bmc_pb_length(opts, k_);

	if (Bmc_check_if_model_was_built(env_nusmv, errstream, false))
	{
		std::cerr << "Model not build " << std::endl;
		free(dump_fname_template);
		return;
	}

	rel_loop = Bmc_Utils_ConvertLoopFromString(get_bmc_pb_loop(opts), NULL);

	/* Checks ltlspecs */

	// DO not call prop
	if (prop_value == 0)
	{
		str_formula = (char *)"FALSE";

		SymbTable_ptr symb_table = SYMB_TABLE(NuSMVEnv_get_value(env_nusmv, ENV_SYMB_TABLE));

		ltlprop = PROP(NULL);
		// Parse and add new LTL property
		int idx = PropDb_prop_parse_and_add(prop_db,
											symb_table,
											str_formula, Prop_Ltl, Nil);
		/* index is ok */
		if (idx != -1)
		{
			nusmv_assert(ltlprop == PROP(NULL));
			ltlprop = PropDb_get_prop_at_index(prop_db, idx);
		}
	}
	else if (prop_value == 1)
	{
		// generate random using SPOT (TODO!)
	}
	else if (prop_value == 2)
	{
		ltlprop = PropDb_get_prop_at_index(prop_db, 0);
	}
	else
	{
		std::cerr << "c Problem in converting LTL spec\n";
		exit(2);
	}

	/* checks that a property was selected: */
	nusmv_assert(ltlprop != PROP(NULL));

	if (Prop_get_status(ltlprop) != Prop_Unchecked)
	{
		free(dump_fname_template);
		return;
	}

	/* ---------------------------------------------------------------------- */
	/* At this point a property was selected                                  */
	/* ---------------------------------------------------------------------- */

	if (opt_cone_of_influence(opts) == true)
	{
		Prop_apply_coi_for_bmc(env_nusmv, ltlprop);
	}

	be_fsm = Prop_get_be_fsm(ltlprop);

	if (be_fsm == (BeFsm_ptr)NULL)
	{
		Prop_set_environment_fsms(env_nusmv, ltlprop);
		be_fsm = Prop_get_be_fsm(ltlprop);
		nusmv_assert(be_fsm != (BeFsm_ptr)NULL);
	}

	rewriter = Prop_Rewriter_create(env_nusmv, ltlprop,
									WFF_REWRITE_METHOD_DEADLOCK_FREE,
									WFF_REWRITER_REWRITE_INPUT_NEXT,
									FSM_TYPE_BE, bdd_enc);

	ltlprop = Prop_Rewriter_rewrite(rewriter);
	be_fsm = Prop_get_be_fsm(ltlprop);

	/* booleanized, negated and NNFed formula: */
	bltlspec = Wff2Nnf(env_nusmv, Wff_make_not(nodemgr, Compile_detexpr2bexpr(bdd_enc,
																			  Prop_get_expr_core(ltlprop))));
	be_enc = BeFsm_get_be_encoding(be_fsm);
	be_mgr = BeEnc_get_be_manager(be_enc);

	/* Start problems generations: */

	/* the loopback value could be depending on the length
		if it were relative: */
	l = Bmc_Utils_RelLoop2AbsLoop(rel_loop, k_);

	/* this is for verbose messages */
	Bmc_Utils_ConvertLoopFromInteger(rel_loop, szLoop, sizeof(szLoop));

	/* checks for loopback vs k compatibility */
	if (Bmc_Utils_IsSingleLoopback(l) && ((l >= k_) || (l < 0)))
	{
		std::cerr << "Error: check loopback k compatibility " << std::endl;
		free(dump_fname_template);
		return;
	}

	/* generates the problem: */

	be_ptr prob = Bmc_Gen_LtlProblem(be_fsm, bltlspec, k_, l, &step_phy_index);
	prob = Bmc_Utils_apply_inlining(be_mgr, prob); /* inline if needed */
	/* Problem is cnf-ed */
	cnf = (Be_Cnf_ptr)NULL;
	/* Problem dumping: */
	if (dump_type != BMC_DUMP_NONE)
	{
		cnf = Be_ConvertToCnf(be_mgr, prob, 1, cnf_alg);

		Bmc_Dump_WriteProblem(be_enc, cnf, ltlprop, k_, l, step_phy_index,
							  dump_type, dump_fname_template);
	}

	dimacsfile = fopen(dump_fname_extension.c_str(), "a");
	Compile_write_model_flat_bool(env_nusmv, dump_fname_template, dimacsfile);

	free(dump_fname_template);

	Prop_Rewriter_update_original_property(rewriter);
	Prop_Rewriter_destroy(rewriter);
	rewriter = NULL;

	if (cnf != (Be_Cnf_ptr)NULL)
	{
		Be_Cnf_Delete(cnf);
		cnf = (Be_Cnf_ptr)NULL;
	}
}

void NuSMV::run_NuSMV(std::string source)
{
	int argc = 3;
	int status;
	std::string tmpfile = "runNuSMV" + std::to_string(rand() % 1000000) + ".txt";
	char tmpfilename[tmpfile.size() + 1];
	strcpy(tmpfilename, tmpfile.c_str());
	// char* tmpfilename = stoc(tmpfile);
	char *argv_tmp[] = {(char *)"./bsaltic", (char *)"-source", tmpfilename};
	char **argv = argv_tmp;
	std::string all_source = source + "quit\n";
	std::ofstream tmpflux(tmpfile.c_str());

	if (tmpflux)
		tmpflux << all_source;
	tmpflux.close();
	/* Finally, call the main function */
	NuSMVCore_main(env_nusmv, argc, argv, &status);
	remove(tmpfilename);
}

void NuSMV::print_properties()
{
	std::string printProp_source = "show_property ";
	printProp_source += "-s -F tabular \n";
	run_NuSMV(printProp_source);
}
