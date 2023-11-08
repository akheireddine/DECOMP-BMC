#include <cerrno>
#include <string>

#include "painless/painless-src/utils/Parameters.h"

#include "src/utils/EnvBMC.hh"
#include "painless/painless-src/painless.h"
#include "src/utils/Dimacs.hh"

// Global main variables
string fname;
EnvBMC *env;

bool read_Filename()
{
	parse_NuSMV_DIMACS_main(fname.c_str(), *env);
	return true;
}

bool run_Painless()
{
	if ((fname.find(".cnf") != std::string::npos) || (fname.find(".dimacs") != std::string::npos))
	{
		read_Filename();
		/// Empty variables
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
		std::cout << "\t-t=<int> int \t   : timeout (default none)." << std::endl;
		std::cout << "\t-max-memory=<int> : memory limit in Go (default none)." << std::endl;
		std::cout << "\t-no-model  \t   : The model will not be printed." << std::endl;
		std::cout << "\t-log  \t        : Enable to print interpolants in output file." << std::endl;
		std::cout << "\t-s={minisat, minisat-old, maple, desat} \t   : Special solver name, default is minisat." << std::endl;
		std::cout << "\t-decomp={batch, rand, bmc} \t   : Decomposition strategy when running desat solver, default is batch." << std::endl;
		std::cout << "\t-nleafs=int \t   : Number of leaves running desat solver with bmc decomp, default is 0." << std::endl;
		std::cout << "\t-nsteps=int \t   : Number of steps grouped in one partition on desat solver with bmc decomp, default is 0." << std::endl;

		std::cout << "\t-c=<int> \t     : Number of cpus for portfolio parallelism, default is 1." << std::endl;
		std::cout << "\t-nc=int \t     : Number of cpus for DeSat partitions, default is 1." << std::endl;
		std::cout << "\t-bmc-on \t     : Enable DeSat solver when running portfolio approach." << std::endl;
		std::cout << "\t-shr-root \t   : Enable globalSolver of DeSat to share conflict clauses with others." << std::endl;
		std::cout << "\t-imp-root \t   : Enable globalSolver of DeSat to receive conflict clauses from others." << std::endl;
		std::cout << "\t-shr-strat=int \t     : Sharing strategy (0: No sharing, 1: SimpleSharing, 2: Hordesat), default is 0." << std::endl;
		std::cout << "\t-lbd-limit=int \t     : Limit of clause's LBD to share, default is 4." << std::endl;
		std::cout << "\t-lbd-limit-itp=int \t : Limit of ITP clause's LBD to share, default is 4." << std::endl;
		std::cout << "\t-d=int \t             : diversification method (0: None, 7: SparseRandom+native), default is 0." << std::endl;
		std::cout << std::endl;
		return 0;
	}

	// Initialize variables
	env = new EnvBMC();
	fname = Parameters::getFilename();

	if (run_Painless())
		std::cout << "c Run BSaLTic done successfully.\n";
	else
		std::cerr << "Run nothing..." << std::endl;
	// delete env;
	std::cout << "c Done." << std::endl;
	return finalResult;
}
