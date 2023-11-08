// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
// RM Martins, 2012

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <time.h>

#include <iostream>
#include <iomanip>

#include <mpi.h>

#include "mpi_commands.h"
#include "mpdesat_master.h"
#include "mpdesat_worker.h"
#include "mpdesat_dynworker.h"
#include "simple_worker.h"
#include "simplest_worker.h"
#include "returnvalue.h"
#include "version.h"

#include "static_worker.h"

#include <Windows.h>

const char *filename=NULL;
int cores = 1; // These are local cores; negative value means that all available cores will be used.
int verbose = 0;
int stats = 0;
//InterpolationMode interpolator = NONE;
InterpolationMode interpolator = MCMILLAN;
//InterpolationMode interpolator = PUDLAK;
int leafs = 2;
int depth = 1;
int inner = 2;
int green = 0;
int time_limit = 0;
int clauses = 0;
int conflicts = 0;
int memory_limit = INT_MAX; // 1048576 = 1MB, 1073741824 = 1GB
int working_nodes = 0;

void treeSplit(int depth, int leafs, int inner, std::vector<int> &master, std::vector< std::vector<int> > &worker)
{

	assert(leafs >= 2);

	int id = 0;
	std::vector<int> worker_list;

	master.push_back(0);
	id++;

	std::vector<int> lvl_list;
	std::vector< std::vector<int> > tree;

	lvl_list.push_back(0);
	tree.push_back(lvl_list);


	//tree generation
	for(int i = 0 ; i < depth; i++)
	{
		lvl_list.clear();
		for(unsigned j = 0; j < tree[i].size(); j++){

			worker_list.clear();
			for(unsigned k = 0; k < (unsigned)inner; k++){

				worker_list.push_back(id);
				lvl_list.push_back(id++);
				master.push_back(tree[i][j]);

			}
			worker.push_back(worker_list);
		}

		if(!lvl_list.empty()) tree.push_back(lvl_list);
	}

	//leaf generation
	lvl_list.clear();
	for(unsigned i=0; i < tree[depth].size(); i++)
	{
		worker_list.clear();
		for(int j=0; j < leafs; j++)
		{
			worker_list.push_back(id);
			lvl_list.push_back(id++);
			master.push_back(tree[depth][i]);
		}
		if(leafs!=0) worker.push_back(worker_list);
		if(!lvl_list.empty()) tree.push_back(lvl_list);
	}

	worker_list.clear();

	for(unsigned i=0; i < tree[depth].size()*leafs; i++){
		worker.push_back(worker_list);	
	}

}

ReturnValue get_options(int argc, char ** argv)
{

	for (int i=1; i<argc; i++)
	{
		if (argv[i][0]!='-')
		{
			if (filename!=NULL)
			{
				std::cout << "Error: Multiple filenames given." << std::endl;
				return S_ERROR;
			}
			else
				filename=argv[i];
		}
		else
		{      
			if (strcmp(argv[i], "-n")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -n requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				cores = atoi(argv[i]);
			}
			else if (strcmp(argv[i], "-i")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -i requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				if (strcmp(argv[i], "m")==0)
					interpolator = MCMILLAN;
				else if (strcmp(argv[i], "im")==0)
					interpolator = INVERSE_MCMILLAN;
				else if (strcmp(argv[i], "p")==0)
					interpolator = PUDLAK;
				else if (strcmp(argv[i], "n")==0)
					interpolator = NONE;
				else 
				{
					std::cout << "Error: Unknown interpolation mode `" << argv[i] << "'." << std::endl;
					return S_ERROR;
				}
			}
			else if(strcmp(argv[i], "-s")==0)
			{
				stats=1;
			}
			else if (strcmp(argv[i], "-v")==0)
			{
				verbose=1;
			}
			else if(strcmp(argv[i], "-g")==0)
			{
				green=1;
			}
			else if (strcmp(argv[i], "-vv")==0)
			{
				verbose=2;
			}
			else if (strcmp(argv[i], "-l")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -l requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				leafs = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-d")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -d requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				depth = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-t")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -t requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				time_limit = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-k")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -k requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				inner = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-cf")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -rc requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				conflicts = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-cl")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -k requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				clauses = atoi(argv[i]);
			}
			else if(strcmp(argv[i], "-m")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -m requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				memory_limit = atoi(argv[i])*1048576; //memory limit expressed in MB
			}
			else if(strcmp(argv[i], "-w")==0)
			{
				if (i+1==argc)
				{
					std::cout << "Error: -w requires an argument." << std::endl;
					return S_ERROR;
				}
				i++;
				working_nodes = atoi(argv[i]);
			}
			else
			{
				std::cout << "Unknown option: " << argv[i] << std::endl;
				return S_ERROR;
			}
		}
	}

	if (!filename)
	{
		std::cout << "No filename given" << std::endl;
		return S_ERROR;
	}

	return S_UNDEFINED;
}

void print_options()
{
	std::cout << "OPTIONS: " << std::endl << std::endl;
	std::cout << "-n [1..n] (number of nodes)" << std::endl;
	std::cout << "-i {m,p,im,n} (interpolation mode (McMillan, Pudlak, Inverse McMillan, None)" << std::endl;
	//std::cout << "-s (more statistics)" << std::endl;
	std::cout << "-v (verbose level 1)" << std::endl;
	std::cout << "-vv (verbose level 2)" << std::endl;
	std::cout << "-g (memory aware version - flat structure)" << std::endl;
	std::cout << "-g -w [1..n] (memory aware version - starts with [1..n] nodes)" << std::endl;
	std::cout << "-g -m [1..n] (memory aware version - limits the memory to [1..n] MB)" << std::endl;
	std::cout << "-l [2..n] (number of leafs in the tree structure)" << std::endl;
	std::cout << "-k [1..n] (number of inner nodes per level in the tree structure)" << std::endl;
	std::cout << "-d [0..n] (depth of the tree structure)" << std::endl;
	//std::cout << "-p [1..100] (percentage of the timelimit allowed during the sequential and problem splitting stages)" << std::endl;
	//std::cout << "-t [1..n] (timelimit for the solving of the formula)" << std::endl << std::endl;
}

int master(int argc, char ** argv, int rank, int size, int masterId, std::vector<int> workersId)
{

	clock_t before = clock();  

	ReturnValue res = get_options(argc, argv);

	ExpressionManager *m = new ExpressionManager();  
	MPDeSATMaster *masterWorker = new MPDeSATMaster(*m,size,rank,filename,cores,interpolator, masterId, workersId);
	masterWorker->setVerbose(verbose);
	masterWorker->setStats(stats);

	if (res==S_UNDEFINED)
	{
		try
		{ 
			bool r = masterWorker->solve();
			res = (r) ? S_SAT : S_UNSAT ;

		}
		catch (std::bad_alloc &)
		{
			std::cout << "MEMORY EXHAUSTED" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
		catch (std::exception &e)
		{
			std::cout << "Exception caught: " << e.what() << std::endl;   
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
		catch (...)
		{
			std::cout << "UNKNOWN EXCEPTION CAUGHT" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
	}

	std::cout << ((res==10)?"SAT":(res==20)?"UNSAT":"UNKNOWN") << std::endl;;
	fflush(stdout);

#ifdef _DEBUG
	delete masterWorker;
	if (m) delete m;
#endif

	return res;
}

int worker(int argc, char ** argv, int rank, int size, int masterId, std::vector<int> workersId)
{
	ReturnValue res = get_options(argc, argv);

	ExpressionManager *m = new ExpressionManager();  
	MPDeSATWorker *worker = new MPDeSATWorker(*m,size,rank,filename,cores,interpolator, masterId, workersId);
	worker->setVerbose(verbose);
	worker->setStats(stats);

	if (res==S_UNDEFINED)
	{
		try
		{
			worker->solve();

		}
		catch (std::bad_alloc &)
		{
			std::cout << "Node " << rank << ": MEMORY EXHAUSTED" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
		catch (std::exception &e)
		{
			std::cout << "Node " << rank << ": Caught exception: " << e.what() << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}    
		catch (...)
		{
			std::cout << "Node " << rank << ": UNKNOWN EXCEPTION." << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}  
	}

#ifdef _DEBUG
	delete worker;
	if (m) delete m;
#endif

	return res;
}

int main(int argc, char ** argv)
{
	int rank, size, res;

	res = MPI_Init(&argc, &argv);  
	if( res != MPI_SUCCESS )
	{
		std::cout << "Error starting MPI program. Terminating." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, res);
	}

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	srand(rank);

	char name[MAX_COMPUTERNAME_LENGTH+1];
	DWORD bufsize = MAX_COMPUTERNAME_LENGTH+1;
	GetComputerNameA(name, &bufsize);
	DWORD pid = GetCurrentProcessId();
	std::cout << "[" << rank << "] running on " << name << ":" << pid << std::endl;

	std::vector<int> master_list;
	std::vector< std::vector<int> > worker_list;

	if (argc == 1 )
	{
		if (rank == 0) print_options();

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return S_UNDEFINED;
	}

	depth = 0;
	leafs = size-1;
	inner = 1;

	get_options(argc, argv); 	

	if (green) {
		if (rank == 0) std::cout << "[0] File: " << filename << std::endl;
		if (working_nodes == 0)
			working_nodes = size;
		SimplestWorker * worker = new SimplestWorker(filename, size, working_nodes, memory_limit);
		worker->setVerbose(verbose);	    	    

		try
		{ 
			bool r = worker->solve();
			res = (r) ? S_SAT : S_UNSAT ;
		}
		catch (std::bad_alloc &)
		{
			std::cout << "[" << rank << "] MAIN MEMORY EXHAUSTED" << std::endl;			
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
		catch (std::exception &e)
		{
			std::cout << "[" << rank << "] Exception caught: " << e.what() << std::endl;   
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}
		catch (...)
		{
			std::cout << "[" << rank << "] UNKNOWN EXCEPTION CAUGHT" << std::endl;
			MPI_Abort(MPI_COMM_WORLD,S_ERROR);
			return S_ERROR;
		}

		if (worker->isMaster())
		{
			worker->getStats();
			std::cout << ((res==10)?"SAT":(res==20)?"UNSAT":"UNKNOWN") << std::endl;;

		}
		fflush(stdout);

#ifdef _DEBUG
		delete worker;
#endif
	}
	else 
	{
		//generates the tree topology based on the number of leafs, leafs per inner nodes and depth of the tree
		//for tree based topologies use: -d n -k 2 -l 2
		//for tree-flat based topologies use: -d m -k n -l n
		//for flat based topologies use: -d 0 -k 1 -l (n-1)
		treeSplit(depth, leafs, inner, master_list, worker_list);

		double num_nodes = 0;
		for(int i = 0; i <= (int)depth; i++)
		{
			num_nodes += pow((double)inner,i);
		}
		num_nodes += pow((double)inner,(int)depth)*leafs;

		if (num_nodes > size)
		{
			if (rank==MPI_MASTER)
			{
				std::cout << "Error: not enough nodes for this tree structure (" << num_nodes << " nodes required)." << std::endl;
				MPI_Abort(MPI_COMM_WORLD,S_ERROR);
				return S_ERROR;
			}
		}


		filename = NULL; //why is the command line being read two times?

		if (rank == MPI_MASTER)
			res = master(argc, argv, rank, size, master_list[rank], worker_list[rank]);
		else    
			worker(argc, argv, rank, size, master_list[rank], worker_list[rank]);


		/* StaticWorker is under development: refactoring of the tree and tree-flat topologies

		ExpressionManager *m = new ExpressionManager(); 
		StaticWorker * worker = new StaticWorker(*m,filename, size, interpolator, master_list[rank], worker_list[rank]);

		try
		{ 
		bool r = worker->solve();
		res = (r) ? S_SAT : S_UNSAT ;
		}
		catch (std::bad_alloc &)
		{
		std::cout << "[" << rank << "] MAIN MEMORY EXHAUSTED" << std::endl;			
		MPI_Abort(MPI_COMM_WORLD,S_ERROR);
		return S_ERROR;
		}
		catch (std::exception &e)
		{
		std::cout << "[" << rank << "] Exception caught: " << e.what() << std::endl;   
		MPI_Abort(MPI_COMM_WORLD,S_ERROR);
		return S_ERROR;
		}
		catch (...)
		{
		std::cout << "[" << rank << "] UNKNOWN EXCEPTION CAUGHT" << std::endl;
		MPI_Abort(MPI_COMM_WORLD,S_ERROR);
		return S_ERROR;
		}

		if (worker->isMaster())
		{
		std::cout << ((res==10)?"SAT":(res==20)?"UNSAT":"UNKNOWN") << std::endl;;
		}
		fflush(stdout);
		*/

	}


MPI_Finalize();

#ifdef _DEBUG
// _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
_CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
_CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
_CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
_CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
_CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);  

//_CrtDumpMemoryLeaks();
//_CrtMemState s1;
//_CrtMemCheckpoint( &s1 );
//_CrtMemDumpStatistics( &s1 );
#endif    

return res;
}