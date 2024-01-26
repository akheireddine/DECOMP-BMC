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

#pragma once

#include "../solvers/Maple.h"
#include "../solvers/MiniSat.h"
#include "../solvers/DeSATSolver.h"
#include "../solvers/SolverInterface.h"
#include "../../../src/utils/EnvBMC.hh"

#include <vector>

using namespace std;

static atomic<int> currentIdSolver(0);

/// Factory to create solvers.
class SolverFactory
{
public:
   /// Instantiate and return a Minisat solver.
   static SolverInterface * createMiniSatSolver();

   /// Instantiate and return a Maple solver.
   static SolverInterface * createMapleSolver();
                     
   static SolverInterface * createDeSATSolver(EnvBMC* env_bmc, bool seq=false);

   /// Instantiate and return a group of Maple solvers.
   static void createDeSATSolvers(int groupSize, EnvBMC* env_bmc,
                                  vector<SolverInterface *> & solvers);

   /// Instantiate and return a group of Minisat solvers.
   static void createMiniSatSolvers(int groupSize,
                                    vector<SolverInterface *> & solvers);

   /// Instantiate and return a group of Maple solvers.
   static void createMapleSolvers(int groupSize,
                                  vector<SolverInterface *> & solvers);

   /// Print stats of a groupe of solvers.
   static void printStats(const vector<SolverInterface *> & solvers);

   /// Apply a sparse diversification on solvers.
   static void sparseDiversification(const vector<SolverInterface *> & solvers);

   /// Apply a random diversification on solvers.
   static void randomDiversification(const vector<SolverInterface *> & solvers,
                                     unsigned int seed);

   /// Apply a sparse and random diversification on solvers.
   static void sparseRandomDiversification(const
                                           vector<SolverInterface *> & solvers);

   /// Apply a native diversification on solvers.
   static void nativeDiversification(const vector<SolverInterface *> & solvers);

   /// Apply a binary value diversification on solvers.
   static void binValueDiversification(const
                                       vector<SolverInterface *> & solvers);
};
