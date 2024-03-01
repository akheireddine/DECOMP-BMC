// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifdef _DEBUG
  #define _CRTDBG_MAP_ALLOC
  #define _CRTDBG_MAP_ALLOC_NEW
  #include <stdlib.h>
  #include <crtdbg.h>
#endif

#include <time.h>

#include <iostream>
#include <iomanip>

#include <ExpressionManager.h>
#include <interpolation_mode.h>
#include <decomposition_mode.h>

#include "../libDeSAT/desat.h"
#include <math.h>
#include <cstring>


//#include <desat_old.h>
//typedef DeSAT_old DeSAT;

#include "version.h"

// Return values
typedef enum { S_SAT=10, S_UNSAT=20, S_ERROR=30, S_UNDEFINED=40 } ReturnValue;

// #include "../MiniSat/Solver.h"

using namespace Desat;

clock_t before = 0;
clock_t load_time = 0;
DeSAT *desat = NULL;


static unsigned intWidth(int i)
{
   if (i == 0)
      return 1;
   
   return (i < 0) + 1 + (unsigned) log10(fabs(i));
}

void printModel_(std::vector<int> & model)
{
  model.push_back(0);
  
  unsigned usedWidth = 0;
  
  for (unsigned i = 0; i < model.size(); i++) {
     if (usedWidth + 1 + intWidth(model[i]) > 80) {
        printf("\n");
        usedWidth = 0;
     }
  
     if (usedWidth == 0) {
        usedWidth += printf("v");
     }
  
     usedWidth += printf(" %d", model[i]);
  }

  printf("\n");
}


void showTime(clock_t totalTime, const char *name, clock_t time)
{
  std::cout << "c " << name << "Time: " << std::setprecision(3) <<
    time/(double)CLOCKS_PER_SEC << " sec" <<
    " (" << std::setw(3) << std::setfill(' ') <<
    100 * time/(double)totalTime << " %)" << std::endl;
}

void showTimes(const clock_t &before, const clock_t &load_time, const DeSAT *desat, const ReturnValue &res)
{
  clock_t totalTime = clock()-before;

  showTime(totalTime, "Load ", load_time);

  if (desat)
  {
    std::cout << "c Interpolants: " << desat->interpolants_imported << std::endl;
	std::cout << "c Solutions: " << desat->solutions_imported << std::endl;
	std::cout << "c All sat: " << desat->all_sat_found << std::endl;
    std::cout << "c Refinement iterations: " << desat->rounds << std::endl;
    showTime(totalTime, "Global ", desat->globalTime);
    showTime(totalTime, "Partitions ", desat->partitionsTime);
    showTime(totalTime, "Import ", desat->importTime);
	showTime(totalTime, "Last Iteration",desat->lastIterationTime);
  }

  showTime(totalTime, "", totalTime);

  std::cout << ((res==10)?"s SATISFIABLE":(res==20)?"s UNSATISFIABLE":"s UNKNOWN") << std::endl;
}

/*
BOOL WINAPI ConsoleHandler(DWORD CEvent)
{
    switch(CEvent)
    {
      case CTRL_BREAK_EVENT:
        {
          std::cout << "c *** INTERMEDIATE STATISTICS ***" << std::endl;
          ReturnValue temp = S_UNDEFINED;
          showTimes(before, load_time, desat, temp);
          break;
        }
      case CTRL_C_EVENT:
      case CTRL_CLOSE_EVENT:
      case CTRL_LOGOFF_EVENT:
      case CTRL_SHUTDOWN_EVENT:
        {
          std::cout << "c *** INTERRUPT ***" << std::endl;
          ReturnValue temp = S_UNDEFINED;
          showTimes(before, load_time, desat, temp);
          exit(S_UNDEFINED);
        }
    }
    return TRUE;
}
*/

void show_help(void) {
  std::cout << std::endl << "Usage: " << DESAT_VERSION_NAME << " [Options] <filename>" << std::endl;
  std::cout << std::endl << "Options:" << std::endl;
  std::cout << "-n <n>           number of partitions" << std::endl;
  std::cout << "-c <n>           number of cores to utilize" << std::endl;
  std::cout << "-i {m,p,im,n}    interpolator (m=McMillan, p=Pudlak, im=Inverse McMillan, n=None)" << std::endl;
  std::cout << "-d {l,r,c,v,b,n} decompositor (l=lazy, r=random, c=cyclic, v=variable, b=BMC, n=None)" << std::endl;
  std::cout << "-s               sequentialize partition evaluation" << std::endl;
  std::cout << "-v, -vv          verbosity" << std::endl;
}

int main(int argc, const char ** argv)
{
  before = clock();

  /*
  if (SetConsoleCtrlHandler( (PHANDLER_ROUTINE)ConsoleHandler,TRUE)==FALSE)
  {
      printf("Unable to install console handler!\n");
      return 30;
  }
  */

  std::cout << "c " << DESAT_VERSION_NAME << " " <<
               DESAT_VERSION_NUMBER << " (compiled " << __DATE__ << " " << __TIME__ ") " <<
               DESAT_COPYRIGHT <<
               std::endl;

  const char *filename=NULL;
  int cores = -1;
  unsigned partitions = 1;
  int verbose = 0;
  bool sequential = false;
  InterpolationMode interpolator = MCMILLAN;
  DecompositionMode decomposition = NODECOMP;
  for (int i=1; i<argc; i++)
  {
    if (argv[i][0]!='-')
    {
      if (filename!=NULL)
      {
        std::cout << "c Error: Multiple filenames given." << std::endl;
        return S_ERROR;
      }
      else
        filename=argv[i];
    }
    else
    {
      if (strcmp(argv[i], "-h") == 0)
      {
        show_help();
        return S_UNDEFINED;
      }
      else if (strcmp(argv[i], "-n")==0)
      {
        if (i+1==argc)
        {
          std::cout << "Error: -n requires an argument." << std::endl;
          return S_ERROR;
        }
        i++;
        partitions = atoi(argv[i]);
      }
      else if (strcmp(argv[i], "-c")==0)
      {
        if (i+1==argc)
        {
          std::cout << "Error: -c requires an argument." << std::endl;
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
      else if (strcmp(argv[i], "-s")==0)
      {
        sequential = true;
      }
      else if (strcmp(argv[i], "-v")==0)
      {
        verbose=1;
      }
      else if (strcmp(argv[i], "-vv")==0)
      {
        verbose=2;
      }
      else if (strcmp(argv[i], "-d")==0)
      {
        if (i+1==argc)
        {
          std::cout << "Error: -d requires an argument." << std::endl;
          return S_ERROR;
        }
        i++;
        if (strcmp(argv[i], "l")==0)
          decomposition = BATCH;
        else if (strcmp(argv[i], "r")==0)
          decomposition = RANDOM;
        else if (strcmp(argv[i], "c")==0)
          decomposition = CYCLE;
        else if (strcmp(argv[i], "v")==0)
          decomposition = VARIABLE;
        else if (strcmp(argv[i], "b")==0)
          decomposition = BMC;
        else if (strcmp(argv[i], "n")==0)
          decomposition = NODECOMP;
        else
        {
          std::cout << "Error: Unknown decomposition mode `" << argv[i] << "'." << std::endl;
          return S_ERROR;
        }
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
    std::cout << "Error: No filename given." << std::endl;
    return S_ERROR;
  }

  ReturnValue res = S_UNDEFINED;
  ExpressionManager *m = new ExpressionManager();

  if (cores<=0)
  {
    std::cout << "Warning: Using all available cores." << std::endl;
    desat = new DeSAT(*m, partitions, decomposition, omp_get_num_procs(), filename);
  }

  else
  {
#ifdef WIN32
    if (cores>64)
    {
      std::cout << "Warning: Number of concurrent threads capped at 64 by the system." << std::endl;
      cores = 64;
    }
#endif
    desat = new DeSAT(*m, partitions, decomposition, cores, filename);
  }

  desat->setVerbose(verbose);
  if(interpolator != NONE) desat->setInterpolator(interpolator);
  clock_t load_time=clock();

  std::vector<int> model;
  try
  {


    if (!desat->readDimacsFile(filename))
      res = S_UNSAT;

	//desat->readDimacsFile(filename);

    load_time = clock() - load_time;

    if (res == S_UNDEFINED)
    {
      bool r = desat->solve();

      res = r ? S_SAT : S_UNSAT;
      
      if (res == S_SAT)
      {
        model = desat->getFinalModel();
      }
    }

  }
  catch (std::bad_alloc &)
  {
    std::cout << "c MEMORY EXHAUSTED" << std::endl;
  }
  catch (std::exception &e)
  {
    std::cout << "c Exception caught: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cout << "c UNKNOWN EXCEPTION CAUGHT" << std::endl;
  }

  showTimes(before, load_time, desat, res);


  if (res == S_SAT )
  {
    printModel_(model);
  }

  #ifdef _DEBUG
  delete m;
  delete desat; // Save some time in release mode.

  // _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
  _CrtSetReportFile(_CRT_WARN, _CRTDBG_FILE_STDOUT);
  _CrtSetReportFile(_CRT_ERROR, _CRTDBG_FILE_STDOUT);
  _CrtSetReportFile(_CRT_ASSERT, _CRTDBG_FILE_STDOUT);
  _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_FILE);
  _CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_FILE);
  _CrtSetReportMode(_CRT_ASSERT, _CRTDBG_MODE_FILE);

  _CrtDumpMemoryLeaks();
  _CrtMemState s1;
  _CrtMemCheckpoint( &s1 );
  _CrtMemDumpStatistics( &s1 );
  #endif

  fflush(stdout);

  return res;
}
