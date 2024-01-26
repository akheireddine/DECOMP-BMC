// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// This file is part of PaInleSS.
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

#include "painless.h"

#include "clauses/ClauseManager.h"
#include "sharing/SimpleSharing.h"
#include "sharing/HordeSatSharing.h"
#include "sharing/Sharer.h"
#include "solvers/SolverFactory.h"
#include "utils/Logger.h"
#include "utils/Parameters.h"
#include "utils/SatUtils.h"
#include "utils/System.h"
#include "working/SequentialWorker.h"
#include "working/Portfolio.h"
#include "solvers/SolverInterface.h"
#include "working/WorkingStrategy.h"

#include <unistd.h>
#include <thread>

/// Is it the end of the search
extern atomic<bool> globalEnding;

/// Final result
extern SatResult finalResult;

/// Model for SAT instances
extern vector<int> finalModel;

/// Array of sharers
extern Sharer **sharers;

/// Size of the array of sharers
extern int nSharers;

// -------------------------------------------
// Declaration of global variables
// -------------------------------------------
atomic<bool> globalEnding(false);

SatResult finalResult = UNKNOWN;

vector<int> finalModel;

Sharer **sharers = NULL;

int nSharers = 0;

// -------------------------------------------

using namespace std;


void run_painless_strategyDeSAT(EnvBMC *env_bmc)
{
   vector<int> cube;
   clock_t t1, t2;
   vector<SolverInterface *> solvers;
   string solverType = Parameters::getParam("s", "desat");
   int nSolvers = Parameters::getIntParam("c", 1);

   if (solverType == "desat")
   {
      // Partition the formula and solved through multiple MiniSat1.14
      solvers.push_back(SolverFactory::createDeSATSolver(env_bmc));
   }
   else if (solverType == "minisat-old")
   {
      // Flat resolution through MiniSat1.14
      SolverFactory::createDeSATSolvers(nSolvers, env_bmc, solvers);
   }
   else if (solverType == "maple")
   {
      // Flat resolution through MapleCOMSPS
      SolverFactory::createMapleSolvers(nSolvers, solvers);
   }
   else
   {
      // Flat resolution through MiniSat2.2.0
      SolverFactory::createMiniSatSolvers(nSolvers, solvers);
   }

   if(Parameters::isSet("diversify")){
      SolverFactory::sparseRandomDiversification(solvers);
   }

   // For parallelism only
   if (Parameters::isSet("bmc-on"))
   {
      solvers.push_back(SolverFactory::createDeSATSolver(env_bmc));
      nSolvers++;
   }

   switch (env_bmc->shrValue)
   {
   case 1:
      sharers = new Sharer *[1];
      sharers[0] = new Sharer(nSolvers, new SimpleSharing(), solvers,
                                       // BMC-D or LZY-D is not a consumer
                                       std::vector<SolverInterface *>(solvers.begin(), solvers.end() - 1));
      break;
   case 2:
      sharers = new Sharer *[1];
      sharers[0] = new Sharer(nSolvers, new HordeSatSharing(), solvers,
                                       // BMC-D or LZY-D is not a consumer
                                       std::vector<SolverInterface *>(solvers.begin(), solvers.end() - 1));
      break;
   default:
      break;
   }

   WorkingStrategy *working = new Portfolio();
   for (size_t i = 0; i < solvers.size(); i++)
      working->addSlave(new SequentialWorker(solvers[i]));

   // Init the management of clauses
   ClauseManager::initClauseManager();
   working->solve(cube);

   int timeout = Parameters::getIntParam("t", -1);
   t1 = clock();
   while (globalEnding == false)
   {
      sleep(1);
      t2 = clock() - t1;
      if (timeout > 0 && (((float)t2) / CLOCKS_PER_SEC >= timeout))
      {
         globalEnding = true;
         std::cout << "c TIMEOUT : " << ((float)t2) / CLOCKS_PER_SEC << endl;
         working->setInterrupt();
      }
   }

   t2 = clock() - t1;
   std::cout << "c TIME SOLVING : " << ((float)t2) / CLOCKS_PER_SEC << endl;

   // Print solver stats
   SolverFactory::printStats(solvers);

   // Print the result and the model if SAT
   if (finalResult == SAT)
   {
      std::cout << "s SATISFIABLE" << endl;

      if (Parameters::isSet("no-model") == false)
         printModel(finalModel);
   }
   else if (finalResult == UNSAT)
      std::cout << "s UNSATISFIABLE" << endl;
   else
      std::cout << "s UNKNOWN" << endl;
}