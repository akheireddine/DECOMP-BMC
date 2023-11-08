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

#include "../clauses/ClauseExchange.h"
#include "../utils/DebugUtils.h"
#include "../utils/Parameters.h"
// #include "../../../src/utils/StatsBMC.hh"


#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <set>


using namespace std;


/// Code for SAT result
enum SatResult
{
	SAT     = 10,
	UNSAT   = 20,
	UNKNOWN = 0
};


/// Code  for the type of solvers
enum SolverType
{
	GLUCOSE   = 0,
	LINGELING = 1,
	MAPLE     = 2,
	MINISAT   = 3,
   PERIPLO   = 4,
   HAMMADI   = 5,
   DESAT     = 6
};


/// Structure for solver statistics
struct SolvingStatistics
{ 
   /// Constructor
	SolvingStatistics()
   {
      propagations = 0;
      decisions    = 0;
      conflicts    = 0;
      restarts     = 0;
      memPeak      = 0;
   }

	unsigned long propagations; ///< Number of propagations.
	unsigned long decisions;    ///< Number of decisions taken.
	unsigned long conflicts;    ///< Number of reached conflicts.
	unsigned long restarts;     ///< Number of restarts.
	double        memPeak;      ///< Maximum memory used in Ko.
};


/// Interface of a solver that provides standard features.
class SolverInterface
{
public:
   /// Load formula from a given dimacs file, return false if failed.
   virtual bool loadFormula(const char* filename) = 0;

   /// Get the number of variables of the current resolution.
   virtual int getVariablesCount() = 0;

   /// Get the number of clauses of the current resolution.
   virtual int getClausesCount() = 0;

   /// Get a variable suitable for search splitting.
   virtual int getDivisionVariable() = 0;

   /// Set initial phase for a given variable.
   virtual void setPhase(const int var, const bool phase) = 0;

   /// Bump activity of a given variable.
   virtual void bumpVariableActivity(const int var, const int times) = 0;

   /// Interrupt resolution, solving cannot continue until interrupt is unset.
   virtual void setSolverInterrupt() = 0;

   /// Remove the SAT solving interrupt request.
   virtual void unsetSolverInterrupt() = 0;

   /// Solve the formula with a given cube.
   virtual SatResult solve(const vector<int> & cube) = 0;

   /// Add a permanent clause to the formula.
   virtual void addClause(ClauseExchange * clause) = 0;

   /// Add a permanent clause to the formula.
   virtual bool addClause(vector<int>& cls) = 0;

   /// Add a list of permanent clauses to the formula.
   virtual void addClauses(const vector<ClauseExchange *> & clauses) = 0;

   /// Add a list of initial clauses to the formula.
   virtual void addInitialClauses(const vector<ClauseExchange *> & clauses) = 0;

   /// Add a learned clause to the formula.
   virtual void addLearnedClause(ClauseExchange * clauses) = 0;
   
   /// Add a list of learned clauses to the formula.
   virtual void addLearnedClauses(const vector<ClauseExchange *> & clauses) = 0;

   /// Get a list of learned clauses.
   virtual void getLearnedClauses(vector<ClauseExchange *> & clauses) = 0;

   /// Request the solver to produce more clauses.
   virtual void increaseClauseProduction() = 0;
   
   /// Request the solver to produce less clauses.
   virtual void decreaseClauseProduction() = 0;

   /// Get solver statistics.
   virtual SolvingStatistics getStatistics() = 0;

   /// Return the model in case of SAT result.
   virtual vector<int> getModel() = 0;

   /// Native diversification.
   virtual void diversify(int id) = 0;

   virtual vector< vector<int> > getFinalAnalysis() = 0;

   virtual bool emptyDataClause() = 0;

   /// Constructor.
   SolverInterface(int solverId, SolverType solverType)
   {
      id    = solverId;
      type  = solverType;
      nRefs = 1;
      //statistics_bmc = new StatsBMC(Parameters::getIntParam("k",-1) + 1);
      nClauses = 0;
   }

   /// Destructor.
   virtual ~SolverInterface()
   {
      //delete statistics_bmc;
   }

   /// Increase the counter of references of this solver.
   void increase()
   {
      nRefs++;
   }

   /// Decrease the counter of references of this solver, delete it if needed.
   void release()
   {
      int oldValue = nRefs.fetch_sub(1);

      if (oldValue - 1 == 0) {
         delete this;
      }
   }

   int getInitialClausesCount()
   {
      return nClauses;
   }

   void printHeader(ofstream& opf)
   {
      opf << "Node,";
      opf << "Time,";
      opf << "Iteration,";
      opf << "DecisionVars,";
      opf << "NClausesBeforeSolve,";
      opf << "NClausesAfterSolve,";
      opf << "NReceived/NSended ITP,";
      opf << "NLearntClauses,";
      opf << "AVGLength ITP,";
      opf << "NNewFreshVars" << std::endl;
   }

   /// Id of this solver.
   int id;

   /// Type of this solver.
   SolverType type;

   /// Number of references pointing on this solver.
   atomic<int> nRefs;

   //std::unordered_set< int > my_vars_set;

   vector< bool > vars_bool;

   int nClauses;

   //StatsBMC * statistics_bmc;

};
