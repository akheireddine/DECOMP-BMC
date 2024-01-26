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

#include "../clauses/ClauseBuffer.h"
#include "../solvers/SolverInterface.h"
#include "../utils/Threading.h"
#include "../../../src/utils/EnvBMC.hh"

#include "libDeSAT/desat.h"
#include "libDeSAT/decomposition_mode.h"

using namespace std;

using namespace Desat;

/// Instance of a DeSAT solver
class DeSATSolver : public SolverInterface
{
public:
   /// Load formula from a given dimacs file, return false if failed.
   bool loadFormula(const char *filename);

   /// Get the number of variables of the current resolution.
   int getVariablesCount();

   /// Get the number of clauses of the current resolution.
   int getClausesCount();

   /// Get a variable suitable for search splitting.
   int getDivisionVariable();

   /// Set initial phase for a given variable.
   void setPhase(const int var, const bool phase);

   /// Bump activity of a given variable.
   void bumpVariableActivity(const int var, const int times);

   /// Interrupt resolution, solving cannot continue until interrupt is unset.
   void setSolverInterrupt();

   /// Remove the SAT solving interrupt request.
   void unsetSolverInterrupt();

   /// Solve the formula with a given cube.
   SatResult solve(const vector<int> &cube);

   /// Add a permanent clause to the formula.
   void addClause(ClauseExchange *clause);

   bool addClause(vector<int> &cls);

   /// Add a list of permanent clauses to the formula.
   void addClauses(const vector<ClauseExchange *> &clauses);

   /// Add a list of initial clauses to the formula.
   void addInitialClauses(const vector<ClauseExchange *> &clauses);

   bool addInitialClauses();

   // Clause distribution between partitions
   bool addInitialClausesToPartitions();

   /// Add a learned clause to the formula.
   void addLearnedClause(ClauseExchange *clause);

   /// Add a list of learned clauses to the formula.
   void addLearnedClauses(const vector<ClauseExchange *> &clauses);

   /// Get a list of learned clauses.
   void getLearnedClauses(vector<ClauseExchange *> &clauses);

   /// Request the solver to produce more clauses.
   void increaseClauseProduction();

   /// Request the solver to produce less clauses.
   void decreaseClauseProduction();

   /// Get solver statistics.
   SolvingStatistics getStatistics();

   /// Return the model in case of SAT result.
   vector<int> getModel();

   /// Native diversification.
   void diversify(int id);

   /// Filter Interpolants to be shared
   bool selectInterpolant(int partition);

   /// Constructor.
   DeSATSolver(int id, EnvBMC *env, bool seq = false);

   /// Destructor.
   virtual ~DeSATSolver();

   vector<vector<int> > getFinalAnalysis() { return {}; };

   bool emptyDataClause() { return env_bmc->nb_clauses <= 0; };

   bool importClauses();

   void update();

protected:
   /// Pointer to a Maple solver.
   DeSAT *solver;

   EnvBMC *env_bmc;

   ExpressionManager *m;
   DecompositionMode d;

   int n_partitions;
   int exported_itp_size;
   int exported_itp_lbd;
   int imported_cls;
   int exported_confl_cls;

   // disable partitioning (flat resolution)
   bool seq;

   /// Buffer used to import clauses (units included).
   ClauseBuffer clausesToImport;
   // TODO! import unit clauses if enabled for gSolver
   ClauseBuffer unitsToImport;

   /// Buffer used to export clauses (units included).
   ClauseBuffer clausesToExport;

   /// Size limit used to share clauses.
   atomic<int> lbdLimit, lbdLimitItp;

   /// Used to stop or continue the resolution.
   atomic<bool> stopSolver;

   /// Callback to export clauses.
   friend void desatExportPartitionClause(void *, std::vector<int> &);
   friend void desatExportgSolverClause(void *, int, std::vector<int> &);
   friend bool desatImportClause(void *, std::vector<int> &, int &);
};
