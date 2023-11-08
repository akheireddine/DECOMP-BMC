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

using namespace std;


// Some forward declatarations for Glucose
namespace Glucose
{
	class ParallelSolver;
	class Lit;
	template<class T> class vec;
   class Clause;
}

/// Instance of a Glucose solver
class GlucoseSyrup : public SolverInterface
{
public:
   /// Load formula from a given dimacs file, return false if failed.
   bool loadFormula(const char* filename);

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
   SatResult solve(const vector<int> & cube);

   /// Add a permanent clause to the formula.
   void addClause(ClauseExchange * clause);

   bool addClause(vector<int>& cls){return false;};

   /// Add a list of permanent clauses to the formula.
   void addClauses(const vector<ClauseExchange *> & clauses);
   
   /// Add a list of initial clauses to the formula.
   void addInitialClauses(const vector<ClauseExchange *> & clauses);

   /// Add a learned clause to the formula.
   void addLearnedClause(ClauseExchange * clause);
   
   /// Add a list of learned clauses to the formula.
   void addLearnedClauses(const vector<ClauseExchange *> & clauses);

   /// Get a list of learned clauses.
   void getLearnedClauses(vector<ClauseExchange *> & clauses);

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
   GlucoseSyrup(int id);
   
   /// Copy constructor.
   GlucoseSyrup(const GlucoseSyrup & other, int id);

   /// Destructor.
   virtual ~GlucoseSyrup();

   void createNewVar(bool polarity = true, bool dvar = true);

   void set_simplification(bool value);

   vector< vector<int> > getFinalAnalysis() { return {}; };

   bool emptyDataClause();

protected:
   /// Pointer to a Glucose solver.
   Glucose::ParallelSolver * solver;

   /// Buffer used to import units.
   ClauseBuffer unitsToImport;

   /// Buffer used to import clauses.
   ClauseBuffer clausesToImport;

   /// Buffer used to import clauses (units included).
   ClauseBuffer clausesToExport;

   /// Buffer used to add permanent clauses.
   ClauseBuffer clausesToAdd;

   /// LBD limit used to share clauses.
   atomic<int> glueLimit;

   /// Callback to export unit clauses.
   friend void glucoseExportUnary(void*, Glucose::Lit &);

   /// Callback to export clauses.
   friend void glucoseExportClause(void *, Glucose::Clause &);

   /// Callback to import unit clauses.
   friend Glucose::Lit glucoseImportUnary (void *);

   /// Callback to import clauses.
   friend bool glucoseImportClause(void *, int *, Glucose::vec<Glucose::Lit> &);
};
