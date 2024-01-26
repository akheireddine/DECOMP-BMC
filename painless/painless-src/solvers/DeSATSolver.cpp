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

#include <iomanip>

// DeSAT includes
#include "libDeSAT/desat.h"

#include "../utils/Logger.h"
#include "../utils/System.h"
#include "../utils/Parameters.h"
#include "../clauses/ClauseManager.h"
#include "DeSATSolver.h"

using namespace Desat;

std::ofstream opf_itp;

void showTime_DeSAT(clock_t totalTime, const char *name, clock_t time)
{
  std::cout << "c " << name << "Time: " << std::setprecision(3) << time / (double)CLOCKS_PER_SEC << " sec"
            << " (" << std::setw(3) << std::setfill(' ') << 100 * time / (double)totalTime << " %)" << std::endl;
}

void showTimes_DeSAT(const clock_t &before, const DeSAT *desat, const SatResult &res,
                     int exported_itp, int exported_itp_size, int exported_itp_lbd, int exported_conflicts_cls, int imported_conflict_cls)
{
  clock_t totalTime = clock() - before;

  std::cout << "c Interpolants: " << desat->interpolants_imported << std::endl;
  std::cout << "c Exported ITP " << exported_itp << std::endl;
  if (exported_itp > 0)
  {
    std::cout << "c Exported ITP cls AVG size " << std::setprecision(3) << exported_itp_size / (double)exported_itp << std::endl;
    std::cout << "c Exported ITP cls AVG LBD " << std::setprecision(3) << exported_itp_lbd / (double)exported_itp << std::endl;
  }
  std::cout << "c Exported CONFLICT clauses " << exported_conflicts_cls << std::endl;
  std::cout << "c Imported CONFLICT clauses " << imported_conflict_cls << std::endl;
  std::cout << "c Solutions: " << desat->solutions_imported << std::endl;
  std::cout << "c All sat: " << desat->all_sat_found << std::endl;
  std::cout << "c Refinement iterations: " << desat->rounds << std::endl;
  showTime_DeSAT(totalTime, "Global ", desat->globalTime);
  showTime_DeSAT(totalTime, "Partitions ", desat->partitionsTime);
  showTime_DeSAT(totalTime, "Import ", desat->importTime);
  showTime_DeSAT(totalTime, "Last Iteration", desat->lastIterationTime);

  showTime_DeSAT(totalTime, "", totalTime);
}

void desatExportPartitionClause(void *issuer, std::vector<int> &cls)
{
  DeSATSolver *ms = (DeSATSolver *)issuer;

  int lbd = ms->solver->getGlobalSolver()->computeLBD(cls);

  if (lbd > ms->lbdLimitItp)
    return;

  ClauseExchange *ncls = ClauseManager::allocClause(cls.size());
  ms->exported_itp_lbd += lbd;

  // Fake glue value
  ncls->lbd = lbd;

  for (int i = 0; i < cls.size(); i++)
    ncls->lits[i] = cls[i];

  ncls->from = ms->id;

  ms->clausesToExport.addClause(ncls);
}

void desatExportgSolverClause(void *issuer, int lbd, std::vector<int> &cls)
{
  DeSATSolver *ms = (DeSATSolver *)issuer;

  if (lbd > ms->lbdLimit)
    return;

  ClauseExchange *ncls = ClauseManager::allocClause(cls.size());

  ncls->lbd = lbd;

  for (int i = 0; i < cls.size(); i++)
    ncls->lits[i] = cls[i];

  ncls->from = ms->id;

  ms->clausesToExport.addClause(ncls);

  ms->exported_confl_cls++;
}

bool desatImportClause(void *issuer, std::vector<int> &mcls, int &sender)
{
  DeSATSolver *ms = (DeSATSolver *)issuer;

  ClauseExchange *cls = NULL;

  if (ms->clausesToImport.getClause(&cls) == false)
    return false;
  mcls.clear();
  mcls.resize(cls->size);
  for (int i = 0; i < cls->size; i++)
    mcls[i] = cls->lits[i];

  sender = cls->from;
  ClauseManager::releaseClause(cls);
  ms->imported_cls++;

  return true;
}

DeSATSolver::DeSATSolver(int id, EnvBMC *env, bool seq_) : SolverInterface(id, DESAT),
                                                           stopSolver(false),
                                                           exported_itp_size(0),
                                                           exported_itp_lbd(0),
                                                           exported_confl_cls(0),
                                                           imported_cls(0),
                                                           env_bmc(env),
                                                           seq(seq_)
{
  lbdLimit = Parameters::getIntParam("lbd-limit", 4);
  lbdLimitItp = Parameters::getIntParam("lbd-limit-itp", 4);
  n_partitions = seq ? 0 : env_bmc->nb_leafs;

  d = RANDOM;
  if (env->decompstrat == "bmc")
    d = BMC;
  else if (env->decompstrat == "batch")
    d = BATCH;

  m = new ExpressionManager();

  string filename = Parameters::getFilename();
  std::string save_inteprolants_filename = filename + "_n" + std::to_string(n_partitions) + "_d" + env->decompstrat + ".txt";

  solver = new DeSAT(*m, n_partitions, d, Parameters::getIntParam("nc", 1), save_inteprolants_filename, env_bmc->nb_clauses);
  solver->setClauseMax(env_bmc->nb_clauses);
  solver->setVariableMax(env_bmc->nb_variables);
  solver->setInterpolator(MCMILLAN);

  // Share conflict clauses coming from globalSolver
  SATSolver *globalSolver = solver->getGlobalSolver();
  globalSolver->setIssuer(this);

  // Enable gSolver (flat or DECOMP) to export its conflict clauses in a parallel context
  if (env_bmc->shrValue > 0 && (Parameters::isSet("shr-root") || seq))
    globalSolver->setCallbackExportRoot(&desatExportgSolverClause);

  // Enable gSolver (flat or DECOMP) to import external clauses in a parallel context
  if (env_bmc->shrValue > 0 && (Parameters::isSet("imp-root") || seq))
    globalSolver->setCallbackImportClauses(&desatImportClause);
}

DeSATSolver::~DeSATSolver()
{
  delete m;
  delete solver;
}

bool DeSATSolver::loadFormula(const char *filename)
{
  return false;
}

// Get the number of variables of the formula
int DeSATSolver::getVariablesCount()
{
  return env_bmc->nb_variables;
}

int DeSATSolver::getClausesCount()
{
  return env_bmc->nb_clauses;
}

// Get a variable suitable for search splitting
int DeSATSolver::getDivisionVariable()
{
  return -1;
}

// Set initial phase for a given variable
void DeSATSolver::setPhase(const int var, const bool phase)
{
  solver->setPhase(var, phase);
}

// Bump activity for a given variable
void DeSATSolver::bumpVariableActivity(const int var, const int times)
{
}

// Interrupt the SAT solving, so it can be started again with new assumptions
void DeSATSolver::setSolverInterrupt()
{
  stopSolver = true;
  if (solver != NULL)
    solver->setInterrupt();
}

void DeSATSolver::unsetSolverInterrupt()
{
  stopSolver = false;
  solver->unsetInterrupt();
}

bool DeSATSolver::importClauses()
{
  SATSolver *globalSolver = solver->getGlobalSolver();
  ClauseExchange *cls = NULL;
  vector<signed> cls_desat;

  while (clausesToImport.getClause(&cls) != false && !stopSolver)
  {
    for (int i = 0; i < cls->size; i++)
      cls_desat.emplace_back(cls->lits[i]);

    ClauseManager::releaseClause(cls);

    if (globalSolver->addClause(cls_desat) == false)
    {
      std::cout << "c unsat when adding imported cls" << endl;
      return false;
    }
    cls = NULL;
    imported_cls++;
  }

  return true;
}

// Diversify the solver
void DeSATSolver::diversify(int id)
{
}

void DeSATSolver::update()
{
  // std::cout << "c Solving with " << getVariablesCount() << " variables, "
  //           << getClausesCount() << " clauses and "
  //           << n_partitions << " partitions on "
  //           << "1 cores." << std::endl;

  SATSolver *globalSolver = solver->getGlobalSolver();

  solver->showDistribution();

  std::cout << "c DeSAT ID (" << id << ")\n";
  assert(globalSolver);

  solver->sharedVariables.update();
  solver->rounds = 0;
}

// Solve the formula with a given set of assumptions
// return 10 for SAT, 20 for UNSAT, 0 for UNKNOWN
SatResult DeSATSolver::solve(const vector<int> &cube)
{
  SATSolver *globalSolver = solver->getGlobalSolver();
  clock_t before = clock();
  int all_sat = 0, exported_itp = 0;
  SatResult res = UNKNOWN;
  bool stop_first_unsat = false;
  std::vector<std::vector<int> > itp_clauses;

  // Flat resolution
  if (seq && n_partitions == 0)
  {
    globalSolver->setVariableMax(env_bmc->nb_variables);
    if (!addInitialClauses())
      return UNSAT;
    bool r = globalSolver->solve();
    solver->globalTime += clock() - before;
    return r ? SAT : UNSAT;
  }

  if (!addInitialClausesToPartitions())
    return UNSAT;

  update();

  while (res == UNKNOWN && !stopSolver)
  {
    solver->rounds++;

    // Run globalSolver
    if (!solver->solveGlobals())
      res = UNSAT;
    else
    {
      stop_first_unsat = false;
      solver->setAssumption();
      // Run Partitions
      for (int i = 0; i < n_partitions && !stop_first_unsat; i++)
      {
        itp_clauses.clear();
        if (stopSolver)
        {
          showTimes_DeSAT(before, solver, res, exported_itp, exported_itp_size, exported_itp_lbd, exported_confl_cls, imported_cls);
          return UNKNOWN;
        }

        // If Partition i is SAT
        if (solver->solvePartition(i))
          all_sat++;
        // If partition is UNSAT
        else
        {
          if (!solver->importInterpolants(itp_clauses, i))
            res = SAT;
          if (env_bmc->log || env_bmc->log_itp_only)
          {
            for (auto cls : itp_clauses)
            {
              for (auto v : cls)
                env_bmc->logFile << v << " ";
              env_bmc->logFile << "0\n";
            }
            env_bmc->logFile.flush();
          }
          if (env_bmc->stop_first_unsat)
            stop_first_unsat = true;
          //&& env_bmc->stop_at_first_unsat;
        }
        // Export clauses interpolants if not empty
        if (!itp_clauses.empty())
        {
          exported_itp++;
          /// EXPORT CLAUSES
          for (auto aclause : itp_clauses)
            desatExportPartitionClause(this, aclause);
          exported_itp_size += itp_clauses.size();
        }
      }
      // When all partitions return SAT
      if (all_sat == n_partitions)
      {
        if (!solver->findDisagreement())
          res = SAT;
        solver->all_sat_found++;
      }
      all_sat = 0;
    }
  }
  showTimes_DeSAT(before, solver, res, exported_itp, exported_itp_size, exported_itp_lbd, exported_confl_cls, imported_cls);
  return res;
}

void DeSATSolver::addClause(ClauseExchange *clause)
{
  std::cerr << "DeSAT NI : addClause" << std::endl;
}

bool DeSATSolver::addClause(vector<int> &cls)
{
  std::cerr << "DeSAT NI : addClause" << std::endl;
  return false;
}

void DeSATSolver::addLearnedClause(ClauseExchange *clause)
{
  if (clause->size == 1)
  {
    unitsToImport.addClause(clause);
  }
  else
  {
    clausesToImport.addClause(clause);
  }
}

void DeSATSolver::addClauses(const vector<ClauseExchange *> &clauses)
{
  std::cerr << "DeSAT NI : addClauses" << std::endl;
}

bool DeSATSolver::addInitialClauses()
{
  assert(n_partitions == 0 || seq);
  for (int p = 0; p < env_bmc->nb_leafs; p++)
  {
    for (auto cls : env_bmc->clauses_partition[p])
    {
      if (solver->addClause(cls) == false)
      {
        std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
        return false;
      }
    }
  }
  for (size_t ind = 0; ind < env_bmc->clauses_g.size(); ind++)
  {
    if (solver->addClause(env_bmc->clauses_g[ind]) == false)
    {
      std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
      return false;
    }
  }
  return true;
}

bool DeSATSolver::addInitialClausesToPartitions()
{
  if (d == BMC)
  {
    // Initialize partitions
    for (int p = 0; p < n_partitions; p++)
    {
      for (auto cls : env_bmc->clauses_partition[p])
      {
        if (d == BMC)
        {
          if (solver->addClause(cls, p) == false)
          {
            std::cout << "c unsat when adding initial cls to partition" << p << std::endl;
            return false;
          }
        }
        else if (d == BATCH || d == RANDOM)
        {
          if (solver->addClause(cls) == false)
          {
            std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
            return false;
          }
        }
      }
    }
    // Initialize gSolver clauses
    for (size_t ind = 0; ind < env_bmc->clauses_g.size(); ind++)
    {
      vector<signed> cls = env_bmc->clauses_g[ind];
      if (solver->addClause(cls, -1) == false)
      {
        std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
        return false;
      }
    }
  }
  if (d == BATCH || d == RANDOM)
  {
    // BATCH and RANDOM decomposition is managed by DeSAT
    for (int p = 0; p < n_partitions; p++)
    {
      for (auto cls : env_bmc->clauses_partition[p])
      {
        if (solver->addClause(cls) == false)
        {
          std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
          return false;
        }
      }
    }
    for (size_t ind = 0; ind < env_bmc->clauses_g.size(); ind++)
    {
      vector<signed> cls = env_bmc->clauses_g[ind];
      if (solver->addClause(cls) == false)
      {
        std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
        return false;
      }
    }
  }
  return true;
}

// bool DeSATSolver::addInitialClausesToPartitions_()
// {
//   // Initialize partitions
//   for (int p = 0; p < n_partitions; p++)
//   {
//     for (auto cls : env_bmc->clauses_partition[p])
//     {
//       if (d == BMC)
//       {
//         if (solver->addClause(cls, p) == false)
//         {
//           std::cout << "c unsat when adding initial cls to partition" << p << std::endl;
//           return false;
//         }
//       }
//       else if (d == BATCH || d == RANDOM)
//       {
//         if (solver->addClause(cls) == false)
//         {
//           std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
//           return false;
//         }
//       }
//     }
//   }
//   // Root gSolver clauses
//   for (size_t ind = 0; ind < env_bmc->clauses_g.size(); ind++)
//   {
//     vector<signed> cls = env_bmc->clauses_g[ind];
//     // N-ARY TREE WITH RANDOM DECOMPOSITION (HAMMADI's APPROACH)
//     if (d == BATCH || d == RANDOM)
//     {
//       if (solver->addClause(cls) == false)
//       {
//         std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
//         return false;
//       }
//     }
//     else
//     {
//       if (solver->addClause(cls, -1) == false)
//       {
//         std::cout << "c unsat when adding initial cls to gSolver" << std::endl;
//         return false;
//       }
//     }
//   }
//   return true;
// }

void DeSATSolver::addInitialClauses(const vector<ClauseExchange *> &clauses)
{
  std::cerr << "DeSAT NI : addInitialClauses" << std::endl;
}

void DeSATSolver::addLearnedClauses(const vector<ClauseExchange *> &clauses)
{
  for (size_t i = 0; i < clauses.size(); i++)
  {
    addLearnedClause(clauses[i]);
  }
}

void DeSATSolver::getLearnedClauses(vector<ClauseExchange *> &clauses)
{
  clausesToExport.getClauses(clauses);
}

void DeSATSolver::increaseClauseProduction()
{
  lbdLimit++;
}

void DeSATSolver::decreaseClauseProduction()
{
  if (lbdLimit > 2)
  {
    lbdLimit--;
  }
}

SolvingStatistics
DeSATSolver::getStatistics()
{
  SolvingStatistics stats;
  solver->getGlobalSolver()->getStats(stats.conflicts, stats.propagations, stats.restarts, stats.decisions);
  stats.memPeak = Desat::memUsedPeak();
  return stats;
}

std::vector<int>
DeSATSolver::getModel()
{
  std::vector<int> model;
  model = solver->getFinalModel();
  return model;
}
