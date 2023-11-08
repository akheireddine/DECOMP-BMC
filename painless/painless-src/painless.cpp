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
#include "working/CubeAndConquer.h"

#include "../../painless_wrapper.hh"

#include <unistd.h>

using namespace std;

void diversification(vector<SolverInterface *> &solvers)
{
   int diversification = Parameters::getIntParam("d", 0);

   switch (diversification)
   {
   case 1:
      SolverFactory::sparseDiversification(solvers);
      break;

   case 2:
      SolverFactory::binValueDiversification(solvers);
      break;

   case 3:
      SolverFactory::randomDiversification(solvers, 2015);
      break;

   case 4:
      SolverFactory::nativeDiversification(solvers);
      break;

   case 5:
      SolverFactory::sparseDiversification(solvers);
      SolverFactory::nativeDiversification(solvers);
      break;

   case 6:
      SolverFactory::sparseRandomDiversification(solvers);
      break;

   case 7:
      SolverFactory::sparseRandomDiversification(solvers);
      SolverFactory::nativeDiversification(solvers);
      break;

   case 0:
      break;
   }
}

void run_painless_strategyDeSAT(EnvBMC *env_bmc)
{
   vector<int> cube;
   clock_t t1, t2;
   vector<SolverInterface *> solvers;
   string solverType = Parameters::getParam("s", "desat");
   int nSolvers = Parameters::getIntParam("c", 1);
   bool desatOn = Parameters::isSet("bmc-on"); // for parallel only

   if (solverType == "maple")
   {
      SolverFactory::createMapleSolvers(nSolvers, solvers);
   }
   else if (solverType == "minisat-old")
   {
      SolverFactory::createDeSATSolvers(nSolvers, env_bmc, solvers);
   }
   else if (solverType == "desat")
      solvers.push_back(SolverFactory::createDeSATSolver(env_bmc));
   else
      SolverFactory::createMiniSatSolvers(nSolvers, solvers);

   diversification(solvers);

   if (desatOn)// for parallel setting
   {
      solvers.push_back(SolverFactory::createDeSATSolver(env_bmc));
      nSolvers++;
   }

   switch (env_bmc->shrValue)
   {
   case 1:
      sharers = new Sharer *[1];
      /// Only std solvers are consumers
      sharers[0] = new Sharer(nSolvers, new SimpleSharing(), solvers, solvers);
      break;
   case 2:
      sharers = new Sharer *[1];
      sharers[0] = new Sharer(nSolvers, new HordeSatSharing(), solvers, solvers);
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
   sharers[0]->printStats();
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
