
#include <string>

#include "painless/painless-src/utils/Parameters.h"

#include "src/utils/EnvBMC.hh"
#include "painless/painless-src/painless.h"
#include "src/mc/NuSMV.hh"
#include "src/utils/Dimacs.hh"

std::string fname;
EnvBMC *env;

/*
	Convert SMV program to CNF
	- Use Tseitin convertion:  set conv parameter to 0
	- Use Sheridan convertion: set conv parameter to 1
*/
bool run_Conversion()
{
	// filename with .smv extension
	if (env->mode == "c")
	{
		Parameters::removeExtensionFromFilename(".smv");
		// Create Model
		NuSMV *mctool = new NuSMV(Parameters::getFilename(), env);
		// Convert SMV instance to CNF
		clock_t tm = clock();
		mctool->generate_cnf_formula(env->k, 2);
		tm = clock() - tm;
		std::cout << "c CONVERSION TO CNF USING " << (env->convertor == 0 ? "TSEITIN" : "SHERIDAN")
				  << ": " << ((float)tm) / CLOCKS_PER_SEC << std::endl;
		delete mctool;
		return true;
	}
	return false;
}

bool run_Painless()
{
	if (env->mode == "p" && ((fname.find(".cnf") != std::string::npos) ||
							 (fname.find(".dimacs") != std::string::npos)))
	{
		/* Read CNF to collect additional information about
		 variables (property vars, time step, ...) */
		parse_NuSMV_DIMACS_main(fname.c_str(), *env);

		// Run DeSAT for formula decomposition through Painless
		run_painless_strategyDeSAT(env);
		return true;
	}

	return false;
}

int main(int argc, char **argv)
{
	// Init parameters
	Parameters::init(argc, argv);

	if (argc < 2 || Parameters::isSet("h"))
	{
		std::cerr << "Missing parameters. At least provide : " << std::endl;
		std::cout << "USAGE: " << argv[0] << " [parameters] filename" << std::endl;
		std::cout << "Parameters:" << std::endl;
		std::cout << "\tfilename  \t      : Filename problem {.smv,.cnf,.dimacs}." << std::endl;
		std::cout << "\t-t=<int>  \t      : timeout (default none)." << std::endl;
		std::cout << "\t-max-memory=<int>     : memory limit in Go (default none)." << std::endl;
		std::cout << "\t-mode={c,p} \t      : two possible modes (c : Convert SMV program to DIMACS. " << std::endl;
		std::cout << "\t\t\t\t\t\t    p : Run solving). Default is p" << std::endl;
		std::cout << std::endl;
		std::cout << "\t Only for 'c' mode: " << std::endl;
		std::cout << "\t    -k=<int> \t        : Bound for the transition relation (bound 0 included : k+1)." << std::endl;
		std::cout << "\t    -conv={0,1} \t: Conversion algorithm (0: Tseitin, 1: Sheridan). Default is 1." << std::endl;
		std::cout << std::endl;

		std::cout << "\t Only for 'p' mode: " << std::endl;
		std::cout << "\t    -no-model  \t                   : The model will not be printed." << std::endl;
		std::cout << "\t    -log  \t                   : Enable to print interpolants in output file." << std::endl;
		std::cout << "\t    -s={maple, minisat-old, minisat, desat} \t   : Specify solver name. Default is desat." << std::endl;
		std::cout << "\t    -decomp={batch, rand, bmc} \t   : Decomposition strategy when running desat solver. Default is batch." << std::endl;
		std::cout << "\t    -nleafs=<int> \t           : Number of leaves running desat solver. Default is 0." << std::endl;
		std::cout << "\t    -nsteps=<int> \t           : Number of steps grouped in one partition on desat solver with 'bmc' decomp. Default is 0." << std::endl;

		std::cout << "\t    -c=<int> \t                   : Number of SAT solvers in the portfolio. Default is 1." << std::endl;
		std::cout << "\t    -bmc-on \t                   : Enable DeSat solver when running portfolio approach." << std::endl;
		std::cout << "\t    -shr-root \t                   : Enable gSolver of DeSat to share conflict clauses with others." << std::endl;
		std::cout << "\t    -imp-root \t                   : Enable gSolver of DeSat to receive conflict clauses from others." << std::endl;
		std::cout << "\t    -shr-strat=<int> \t           : Sharing strategy (0: No sharing, 1: SimpleSharing, 2: Hordesat). Default is 0." << std::endl;
		std::cout << "\t    -lbd-limit=<int> \t           : Limit of clause's LBD to share. Default is 4." << std::endl;
		std::cout << "\t    -lbd-limit-itp=<int> \t   : Limit of ITP clause's LBD to share. Default is 4." << std::endl;
		std::cout << "\t    -diversify \t                   : Enable diversification (SparseRandom)." << std::endl;
		std::cout << std::endl;
		return 1;
	}

	// Initialize variables
	env = new EnvBMC();
	fname = Parameters::getFilename();

	if (run_Conversion())
		std::cout << "c Conversion done successfully.\n";
	else
	{
		if (run_Painless())
			std::cout << "c Run DECOMP-BMC done successfully.\n";
		else
			std::cerr << "Run nothing..." << std::endl;
	}

	std::cout << "c Done." << std::endl;

	delete env;
	return 0;
}
