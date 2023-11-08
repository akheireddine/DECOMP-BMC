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

using namespace std;

// Some forward declatarations for MapleCOMSPS
namespace MapleCOMSPS
{
   class SimpSolver;
   class Lit;
   template <class T>
   class vec;
}

/// Instance of a Maple solver
class Maple : public SolverInterface
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

   bool addInitialClauses_(const vector<ClauseExchange *> &clauses);

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

   /// Constructor.
   Maple(int id, EnvBMC *env = NULL);

   /// Destructor.
   virtual ~Maple();

   void createNewVar(bool polarity = true, bool dvar = true);

   void set_simplification(bool value);

   void disable_callback();

   vector<vector<int>> getFinalAnalysis();

   bool emptyDataClause();

   std::vector<int> getAssignments();

   MapleCOMSPS::SimpSolver *getSolver() { return solver; };

   bool selectConflitClause(MapleCOMSPS::vec<MapleCOMSPS::Lit> &clause_confl);

   bool selectConflitClause_LBDI_OptimalConfig(MapleCOMSPS::vec<MapleCOMSPS::Lit> &clause_confl, int lbd);

protected:
   /// Pointer to a Maple solver.
   MapleCOMSPS::SimpSolver *solver;

   /// Environment variable of BMC information
   EnvBMC *env_bmc;

   /// Buffer used to import clauses (units included).
   ClauseBuffer clausesToImport;
   ClauseBuffer unitsToImport;

   /// Buffer used to export clauses (units included).
   ClauseBuffer clausesToExport;

   /// Buffer used to add permanent clauses.
   ClauseBuffer clausesToAdd;

   /// Size limit used to share clauses.
   atomic<int> lbdLimit;

   atomic<int> shrStrat;

   atomic<int> memLimit;

   atomic<int> code_filter_confl;

   atomic<int> exportClauses;
   /// Used to stop or continue the resolution.
   atomic<bool> stopSolver;

   /// Callback to export/import clauses.
   friend MapleCOMSPS::Lit mapleImportUnit(void *);
   friend bool mapleImportClause(void *, int *, MapleCOMSPS::vec<MapleCOMSPS::Lit> &);
   friend void mapleExportClause(void *, int, MapleCOMSPS::vec<MapleCOMSPS::Lit> &);
};
