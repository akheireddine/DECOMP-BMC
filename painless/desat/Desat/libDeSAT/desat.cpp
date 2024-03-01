// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifdef _DEBUG
#define _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC_NEW
#include <stdlib.h>
#include <crtdbg.h>
#endif

#include <stdlib.h>

#include <sstream>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <omp.h>

#include "cadical.h"
#include "decomposition_batch.h"
#include "decomposition_cycle.h"
#include "decomposition_rnd.h"
#include "decomposition_variables.h"
#include "decomposition_BMC.h"
#include "interpolator_m.h"
#include "interpolator_p.h"
#include "interpolator_im.h"

#include "desat.h"

using namespace Desat;

// #include <psapi.h> FIXME
#pragma comment(lib, "psapi.lib") //added

typedef unsigned long DWORD_PTR;


DeSAT::DeSAT(ExpressionManager &m, unsigned partitions, DecompositionMode decomposition, 
             unsigned cs, std::string filename) : SATSolver(m),
                                                  globalTime(0),
                                                  partitionsTime(0),
                                                  importTime(0),
                                                  lastIterationTime(0),
                                                  maxVar(0),
                                                  interpolants_imported(0),
                                                  d(NULL),
                                                  globalSolver(NULL),
                                                  early_stop(false),
                                                  solutions_imported(0),
                                                  all_sat_found(0),
                                                  filename(filename),
                                                  decompositionMode(decomposition)
{
  init(partitions, cs);
}

void DeSAT::init(unsigned ps, unsigned cs)
{
  n_partitions = ps;
  n_cores = cs;

  if (n_partitions == 0)
      globalSolver = new Cadical(m, false);
  else
  {
    switch (decompositionMode)
    {
    case BATCH:
      d = new BatchDecomposition();
      break;
    case CYCLE:
      d = new CycleDecomposition();
      break;
    case RANDOM:
      d = new RandomDecompositon();
      break;
    case VARIABLE:
      d = new VariableDecomposition();
      break;
    case BMC:
      d = new BMCDecomposition();
      break;
    default:
      throw std::runtime_error("DecompositionMode: init");
      break;
    }
    d->setPartitions(n_partitions);
    assumptions.resize(n_partitions);
    interpolants.resize(n_partitions, m.mkNil());
    for (unsigned i = 0; i < n_partitions; i++)
      partitions.push_back(new Partition(*(new ExpressionManager()), sharedVariables, i, verbosity));
    
    globalSolver = new Cadical(m, sharedVariables, false);
  }
  // FIXME not working on macos
  omp_set_num_threads(n_cores);
}

DeSAT::~DeSAT(void)
{
  if (n_partitions > 1)
    for (unsigned i = 0; i < n_partitions; i++)
    {
      delete &partitions[i]->em();
      delete partitions[i];
    }
  if (d)
    delete d;
  if (globalSolver)
    delete globalSolver;
}

void DeSAT::setVerbose(int v)
{
  for (unsigned i = 0; i < partitions.size(); i++)
    partitions[i]->setVerbose(v);
  if (globalSolver)
    globalSolver->setVerbose(v);
  SATSolver::setVerbose(v);
}

void DeSAT::setInterpolator(InterpolationMode i)
{
  interpolationMode = i;
  if (n_partitions <= 1)
    return;

  for (unsigned j = 0; j < n_partitions; j++)
    partitions[j]->setInterpolationMode(i);
}

void DeSAT::setDecompositor(DecompositionMode i)
{
  decompositionMode = i;
}

void DeSAT::setConstraintMode(ConstraintMode constr)
{
  if (globalSolver)
    globalSolver->setConstraintMode(constr);
  for (unsigned i = 0; i < partitions.size(); i++)
    partitions[i]->setConstraintMode(constr);
}

bool DeSAT::addClause(const std::vector<signed> &literals)
{
  assert(n_cores != 0);
  if (n_partitions == 0){
    return globalSolver->addClause(literals);
  }
  int w = d->where(literals);

  if (w == -1)
  {
    // assert(false); // This should only happen with a VariableDecomposition!
    return globalSolver->addClause(literals);
  }
  else
    return partitions[w]->addClause(literals);
}

bool DeSAT::addClause(const std::vector<signed> &literals, int w)
{
  assert(n_cores != 0);

  if (n_partitions == 0)
    return globalSolver->addClause(literals);

  if (w == -1)
    return globalSolver->addClause(literals);
  else
    return partitions[w]->addClause(literals);
}

bool DeSAT::addUnit(signed l)
{
  throw std::runtime_error("NYI: addUnit");
}

void DeSAT::setVariableMax(unsigned n)
{
  if (n < maxVar)
    return;

  maxVar = n;

  assert(n_cores != 0);

  globalSolver->setVariableMax(n);

  if (n_partitions > 1)
  {
    d->setVariableMax(n);
    for (unsigned i = 0; i < partitions.size(); i++)
      partitions[i]->setVariableMax(n);
  }
}

void DeSAT::setClauseMax(unsigned n)
{
  //if (n<n_partitions)
  //{
  //  std::cout << "Warning: less clauses than cores; reducing the number of partitions to " << n << std::endl;
  //  n_partitions = n;
  //  d->setPartitions(n_partitions);
  //  while (partitions.size() > n)
  //  {
  //    delete partitions.back();
  //    partitions.pop_back();
  //  }
  //}

  if (n_partitions > 1)
  {
    d->setClauseMax(n);
    for (unsigned i = 0; i < partitions.size(); i++)
      partitions[i]->setClauseMax(n);
  }
}

signed DeSAT::addVar(void)
{
  throw std::runtime_error("NYI: adding variables");
  return 0;
}

unsigned DeSAT::numClauses(void) const
{
  unsigned res = 0;
  for (unsigned i = 0; i < partitions.size(); i++)
    res += partitions[i]->numClauses();
  res += globalSolver->numClauses();
  return res;
}

unsigned DeSAT::numVars(void)
{
  return maxVar;
}

bool DeSAT::solve(const std::vector<signed> &assumptions)
{
  throw std::runtime_error("NYI: solving under assumptions.");
  return false;
}

ModelValue DeSAT::get(signed l)
{
  throw std::runtime_error("NYI: model extraction");
  return M_UNDEF;
}

Expression DeSAT::getInterpolant(const std::vector<signed> &A)
{
  throw std::runtime_error("NYI: getInterpolant");
  return m.mkTrue();
}

void DeSAT::printInterpolant(int partition, std::ofstream &op)
{
  if (partition > n_partitions)
    throw std::runtime_error("ERR: printInterpolant");
  CExpression itp = interpolants[partition];
  std::string t = m.toString(itp);
  if (op)
    op << "c " << t.c_str() << std::endl;
}

std::vector<int> DeSAT::getFinalModel(void)
{
  std::vector<int> model;
  std::vector<int> model_g_solver = globalSolver->getModelVector();

  if (n_partitions == 0)
  {
    for (unsigned i = 0; i < model_g_solver.size(); i++)
    {
      int v = std::abs(model_g_solver[i]);
      ModelValue mv = globalSolver->get(v);
      if (mv != M_UNDEF)
        model.emplace_back(model_g_solver[i]); //mv == M_TRUE? v : (-v) );
    }
    return model;
  }
  std::vector<bool> visited_variable(maxVar + 1, false);
  for (unsigned i = 0; i < model_g_solver.size(); i++)
  {
    int v = std::abs(model_g_solver[i]);
    if ((visited_variable[v] == false)) //&& sharedVariables.isShared(v))
    {
      visited_variable[v] = true;
      ModelValue mv = globalSolver->get(v);
      if (mv != M_UNDEF)
        model.emplace_back(model_g_solver[i]); //mv == M_TRUE? v : (-v) );
    }
  }
  for (int i = 0; i < n_partitions; i++)
  {
    SATSolver *solver_partition = partitions[i]->getSolver();
    std::vector<int> model_partition = solver_partition->getModelVector();
    for (unsigned j = 0; j < model_partition.size(); j++)
    {
      int v = std::abs(model_partition[j]);
      if ((visited_variable[v] == false)) //&& !sharedVariables.isShared(v))
      {
        visited_variable[v] = true;
        ModelValue mv = solver_partition->get(v);
        if (mv != M_UNDEF)
          model.emplace_back(model_partition[j]); //mv == M_TRUE? v : (-v) );
      }
    }
  }
  std::sort(model.begin(), model.end());
  return model;
}

void DeSAT::clearNewClauses(void)
{
  throw std::runtime_error("NYI: clearNewClauses");
}

void DeSAT::addNewClauses(void)
{
  throw std::runtime_error("NYI: addNewClauses");
}

void DeSAT::getConflict(std::vector<signed> &out)
{
  throw std::runtime_error("NYI: getConflict");
}

std::vector<std::vector<signed>> DeSAT::getNewClauses(void)
{
  throw std::runtime_error("NYI: getNewClauses");
  return {};
}

std::vector<int> DeSAT::getModelVector(void)
{
  throw std::runtime_error("NYI: model extraction");
  return {};
}

Expression DeSAT::getModel(void)
{
  throw std::runtime_error("NYI: model extraction");
  return m.mkFalse();
}

bool DeSAT::addConstraint(CExpression &e)
{
  throw std::runtime_error("NYI: addConstraint");
  return false;
}

signed DeSAT::addExtension(CExpression &e)
{
  throw std::runtime_error("NYI: addExtension");
  return 0;
}

void DeSAT::showDistribution(void) const
{
  /*if (verbosity==0)
    return;*/

  std::cout << "c Distribution clauses: ";

  for (unsigned i = 0; i < partitions.size(); i++)
    std::cout << partitions[i]->numClauses() << " ";
  std::cout << " [" << globalSolver->numClauses() << "]" << std::endl;

  std::cout << "c Solving with " << maxVar << " variables, " << numClauses() << " clauses and " << n_partitions << " partitions on " << n_cores << " cores. " << std::endl;

  if (n_partitions > 16)
    return;

  int shared_count_all = n_partitions + 1;
  std::vector<std::vector<unsigned>> shared_count(shared_count_all);
  unsigned total = maxVar;
  int global_solver_id = n_partitions;
  bool sharing_variables_w_gs = globalSolver->numClauses() > 0;

  for (unsigned v = 1; v <= maxVar; v++)
  {
    bool exclusive = true;

    for (unsigned i = 0; i < n_partitions; i++)
    {
      for (unsigned j = 0; j < n_partitions; j++)
      {
        if (i == j)
          continue;

        if (shared_count[i].size() != shared_count_all)
          shared_count[i].resize(shared_count_all, 0);

        if (sharedVariables.isShared(v, i, j))
        {
          exclusive = false;
          shared_count[i][j]++;
        }
      }
    }
    // Shared variables with globalSolver
    if (sharing_variables_w_gs)
    {
      if (shared_count[global_solver_id].size() != shared_count_all)
        shared_count[global_solver_id].resize(shared_count_all, 0);

      for (unsigned i = 0; i < n_partitions; i++)
      {
        if (shared_count[i].size() != shared_count_all)
          shared_count[i].resize(shared_count_all, 0);

        if (sharedVariables.isShared(v, i, global_solver_id))
        {
          exclusive = false;
          shared_count[i][global_solver_id]++;
          shared_count[global_solver_id][i]++;
        }
      }
    }
    if (exclusive)
      total--;
  }

  printf("c Global variables: %d (%.2f%%)\n", total, 100 * total / (double)maxVar);

  printf("c Sharing matrix (%%):\n");
  printf("c Partition ");
  for (unsigned i = 0; i < n_partitions; i++)
    printf("    %02d", i);

  if (sharing_variables_w_gs)
    printf("    0G");

  printf("\n");

  for (unsigned i = 0; i < n_partitions; i++)
  {
    printf("c       %02d  ", i);
    for (unsigned j = 0; j < n_partitions; j++)
    {
      printf(" % 3.2f", (100 * shared_count[i][j]) / (double)maxVar);
    }
    if (sharing_variables_w_gs)
      printf(" % 3.2f", (100 * shared_count[i][global_solver_id]) / (double)maxVar);
    printf("\n");
  }

  if (sharing_variables_w_gs)
  {
    printf("c       0G  ");
    for (unsigned i = 0; i < shared_count_all; i++)
      printf(" % 3.2f", (100 * shared_count[global_solver_id][i]) / (double)maxVar);
    printf("\n");
  }

  fflush(stdout);
}

void DeSAT::setPhase(const int var, const bool phase)
{
  globalSolver->setPhase(var, phase);
  for (int pid = 0; pid < (signed)n_partitions; pid++)
    partitions[pid]->setPhase(var, phase);
}

bool DeSAT::solve(void)
{
  std::cout << "c Solving with " << maxVar << " variables, " << numClauses() << " clauses and " << n_partitions << " partitions on " << n_cores << " cores. " << std::endl;
  clock_t before = clock();
  std::vector<std::vector<int>> cls_list;

  if (n_partitions == 0)
  {
    bool r = globalSolver->solve();
    globalTime += clock() - before;
    return r;
  }

  // if( verbosity > 0 )
  showDistribution();

  assert(globalSolver);

  sharedVariables.update();

  rounds = 0;
  ModelValue res = M_UNDEF;

  while ((res == M_UNDEF) && !early_stop)
  {
    if (verbosity == 2)
    {
      print("Press any key to continue...");
      ::getchar();
    }

    rounds++;

    //PROCESS_MEMORY_COUNTERS pmc;
    //GetProcessMemoryInfo( GetCurrentProcess(), &pmc, sizeof(pmc));

    //std::cout << " WorkingSetSize: " << pmc.WorkingSetSize / 1048576 << " MB" << std::endl;

    // for (unsigned i = 0; i < 79; i++)
    //   print("=");
    // print("\nRound %d\n", rounds);
    // for (unsigned i = 0; i < 79; i++)
    //   print("-");
    // print("\n");
    if (!solveGlobals())
      res = M_FALSE;
    else
    {
      if (solvePartitions())
      {
        all_sat_found++;
        if (!findDisagreement())
          res = M_TRUE;
        importInterpolants(cls_list); // deletes the interpolant objects
      }
      else if (!importInterpolants(cls_list))
        res = M_TRUE; // SAT!
    }
  }

  //  for(signed i = 0; i < 1000; i++)
  //std::cout << "i: " << i << " " << m.toString(i) << std::endl;

  return (early_stop && (res != M_TRUE)) ? false : res == M_TRUE;
}

void DeSAT::setInterrupt()
{
  globalSolver->setInterrupt();

  for (int i = 0; i < n_partitions; i++)
    partitions[i]->setInterrupt();
}

void DeSAT::unsetInterrupt()
{
  globalSolver->unsetInterrupt();

  for (int i = 0; i < n_partitions; i++)
    partitions[i]->unsetInterrupt();
}

bool DeSAT::solveGlobals(void)
{
  clock_t before = clock();
  print("Finding global assignment...\n");
  print("Trail size: %d\n", trail.size());

  for (unsigned i = 0; i < trail.size(); i++)
  {
    print(" %d", trail[i]);
    if (reasons[i])
      print("!");
  }
  print("\n");

  bool r = globalSolver->solve(trail);

  while (!r && trail.size() > 0 && !early_stop)
  {
    std::vector<signed> temp;
    ((Cadical *)globalSolver)->getConflict(temp);

    signed last = 0;
    bool reason = true;
    bool inConflict = false;

    print("GLOBAL CONFLICT:");
    if (temp.size() == 0)
      print(" EMPTY");
    else
      for (unsigned i = 0; i < temp.size(); i++)
        print(" %d", temp[i]);
    print("\n");

    print("TRAIL:");
    for (unsigned i = 0; i < trail.size(); i++)
    {
      print(" %d", trail[i]);
      if (reasons[i])
        print("!");
    }
    print("\n");

    inConflict = std::find(temp.begin(), temp.end(), -last) != temp.end();
    while (trail.size() > 0 && !inConflict)
    {
      last = trail.back();
      trail.pop_back();
      reason = reasons.back();
      reasons.pop_back();
      inConflict = std::find(temp.begin(), temp.end(), -last) != temp.end();
    }

    if (inConflict && !reason)
    {
      trail.push_back(-last);
      reasons.push_back(true);
    }

    if (trail.size() == 0 && reason)
      r = false;
    else
    {
      print("TRAIL:");
      for (unsigned i = 0; i < trail.size(); i++)
      {
        print(" %d", trail[i]);
        if (reasons[i])
          print("!");
      }
      print("\n");
      r = globalSolver->solve(trail);
    }
  }
  if (!early_stop)
  {
    globalTime += clock() - before;
    lastIterationTime = clock() - before;
  }

  if (!r)
  {
    print("No global assignment.\n");
  }
  else if (verbosity > 1)
  {
    printf("Global Assignment:");
    for (int v = 1; v <= (signed)maxVar; v++)
      if (sharedVariables.isShared(v))
      {
        ModelValue mv = globalSolver->get(v);
        if (mv != M_UNDEF)
          printf(" %d", mv == M_TRUE ? v : -v);
      }
    printf("\n");
  }

  return r;
}

void DeSAT::setAssumption(void)
{
  print("Extending assignment...\n");

  assumptions.clear();
  for (signed v = 1; v <= (signed)maxVar; v++)
  {
    if (sharedVariables.isShared(v))
    {
      ModelValue mv = globalSolver->get(v);
      if (mv == M_UNDEF)
        continue;
      assumptions.push_back(mv == M_TRUE ? v : -v);
    }
  }

  if (verbosity > 0)
  {
    print("Assumptions: ");
    for (unsigned i = 0; i < assumptions.size(); i++)
      print(" %d", assumptions[i]);
    print("\n");
  }
}

bool DeSAT::solvePartition(int pid)
{
  clock_t before = clock();
  bool all_sat = true;
  have_error = false;
  have_bad_alloc = false;

  if (n_cores == 1)
  {
    interpolants[pid] = m.mkNil();
    Expression t = partitions[pid]->getInterpolant(assumptions);

    // std::cout << "T = " << t << " (" << partitions[pid]->em().get(t).left << "," << partitions[pid]->em().get(t).right << ")" << std::endl;

    if (!partitions[pid]->em().isTrue(t))
    {
      //std::cout << "duplicate" << std::endl;
      interpolants[pid] = m.duplicate(t, partitions[pid]->em());
      //std::cout << "done" << std::endl;
      all_sat = false;
      if (early_stop)
        return false;
    }
  }
  else
  {

    // #pragma omp parallel for default(shared) num_threads(n_cores)
    if (early_stop)
      return false;

    interpolants[pid] = m.mkNil();
    interpolants[pid] = partitions[pid]->getInterpolant(assumptions);

    Expression t = interpolants[pid];
    if (!partitions[pid]->em().isTrue(t))
    {
      try
      {
        //std::cout << "duplicate" << std::endl;
        interpolants[pid] = m.duplicate(t, partitions[pid]->em());
        //std::cout << "done" << std::endl;
        all_sat = false;
      }
      catch (std::bad_alloc &e)
      {
        have_bad_alloc = true;
        ba_exception = e;
      }
      catch (std::runtime_error &e)
      {
        have_error = true;
        exception = e;
      }
      catch (...)
      {
        std::cout << "UNRECOVERABLE ERROR" << std::endl;
      }
    }
  }

  partitionsTime += clock() - before;

  if (have_error)
    throw exception;
  else if (have_bad_alloc)
    throw ba_exception;

  // printf("Solution is %d  (",all_sat);
  // Cadical *s = (Cadical *) partitions[pid]->getSolver();
  // printf("status %d)\n", s->status());

  return all_sat;
}

bool DeSAT::solvePartitions(void)
{
  clock_t before = clock();
  print("Extending assignment...\n");
  bool all_sat = true;

  assumptions.clear();
  for (signed v = 1; v <= (signed)maxVar; v++)
  {
    if (sharedVariables.isShared(v))
    {
      ModelValue mv = globalSolver->get(v);
      if (mv == M_UNDEF)
        continue;
      assumptions.push_back(mv == M_TRUE ? v : -v);
    }
  }

  if (verbosity > 0)
  {
    print("Assumptions: ");
    for (unsigned i = 0; i < assumptions.size(); i++)
      print(" %d", assumptions[i]);
    print("\n");
  }

  have_error = false;
  have_bad_alloc = false;

  if (n_cores == 1)
  {

    for (int pid = 0; pid < (signed)n_partitions; pid++)
    {
      interpolants[pid] = m.mkNil();
      Expression t = partitions[pid]->getInterpolant(assumptions);

      // std::cout << "T = " << t << " (" << partitions[pid]->em().get(t).left << "," << partitions[pid]->em().get(t).right << ")" << std::endl;

      if (!partitions[pid]->em().isTrue(t))
      {
        //std::cout << "duplicate" << std::endl;
        interpolants[pid] = m.duplicate(t, partitions[pid]->em());
        //std::cout << "done" << std::endl;
        all_sat = false;
        if (early_stop)
          break;
      }
    }
  }
  else
  {

#pragma omp parallel for default(shared) num_threads(n_cores)
    for (int pid = 0; pid < (signed)n_partitions; pid++)
    {
      // FIXME Not working like that in linux
      //DWORD_PTR mask = (0x01 << omp_get_thread_num());
      //SetThreadAffinityMask( GetCurrentThread(), mask );

      interpolants[pid] = m.mkNil();
      interpolants[pid] = partitions[pid]->getInterpolant(assumptions);

      //std::cout << "T = " << t << " (" << partitions[pid]->em().get(t).left << "," << partitions[pid]->em().get(t).right << ")" << std::endl;

      //  if (!partitions[pid]->em().isTrue(t))
      //  {
      //    #pragma omp critical
      //    {
      //      // Exceptions cannot pass a critical section, which can
      //      // cause deadlocks. Therefore we have to catch them manually.
      //      try
      //      {
      //        //std::cout << "duplicate" << std::endl;
      //        interpolants[pid] = m.duplicate(t, partitions[pid]->em());
      //        //std::cout << "done" << std::endl;
      //        all_sat = false;
      //}
      //      catch (std::bad_alloc &e)
      //      {
      //        have_bad_alloc = true;
      //        ba_exception = e;
      //      }
      //      catch (std::runtime_error &e)
      //      {
      //        have_error = true;
      //        exception = e;
      //      }
      //      catch (...)
      //      {
      //        std::cout << "UNRECOVERABLE ERROR" << std::endl;
      //      }
      //    }
      //  }
    }

    for (int i = 0; i < (signed)n_partitions; i++)
    {
      if (early_stop)
        return false;
      Expression t = interpolants[i];
      if (!partitions[i]->em().isTrue(t))
      {
        // Exceptions cannot pass a critical section, which can
        // cause deadlocks. Therefore we have to catch them manually.
        try
        {
          //std::cout << "duplicate" << std::endl;
          interpolants[i] = m.duplicate(t, partitions[i]->em());
          //std::cout << "done" << std::endl;
          all_sat = false;
        }
        catch (std::bad_alloc &e)
        {
          have_bad_alloc = true;
          ba_exception = e;
        }
        catch (std::runtime_error &e)
        {
          have_error = true;
          exception = e;
        }
        catch (...)
        {
          std::cout << "UNRECOVERABLE ERROR" << std::endl;
        }
      }
    }
  }

  partitionsTime += clock() - before;

  if (have_error)
    throw exception;
  else if (have_bad_alloc)
    throw ba_exception;

  return all_sat;
}

bool DeSAT::importInterpolants(std::vector<std::vector<int>> &itp_cls_export)
{
  clock_t before = clock();
  bool did_something = false;

  for (unsigned i = 0; i < n_partitions; i++)
  {
    CExpression itp = interpolants[i];
    if (m.isNil(itp))
    {
      solutions_imported++;
      continue;
    }

    if (!m.isTrue(itp))
    {
      // if (verbosity>0)
      {
        std::string t = m.toString(itp);
        // print("Global interpolant import: %s\n", t.c_str());
        //std::string filename = "filename";
        //std::stringstream convert; // stringstream used for the conversion
        //convert << rounds;
        //convert << i;
        //filename += convert.str();
        //m.toDot(itp,filename);
      }
      globalSolver->clearNewClauses();
      globalSolver->addConstraint(itp);
      std::vector<std::vector<signed>> itp_clauses = globalSolver->getNewClauses();
      itp_cls_export = itp_clauses; // TODO CADICAL ? est-ce necessaire

      // for (unsigned i = 0; i < itp_clauses.size(); i++)
      // {
      //   std::vector<int> cls;
      //   for (unsigned j = 0; j < itp_clauses[i].size(); j++)
      //   {
      //     cls.emplace_back(itp_clauses[i][j]);
      //   }
      //   itp_cls_export.emplace_back(cls);
      // }

      globalSolver->addNewClauses();
      interpolants_imported++;
      did_something = true;
    }
    interpolants[i] = m.mkNil();
  }

  clock_t after = clock();
  importTime += after - before;
  return did_something;
}

std::unordered_set<signed> DeSAT::getInterpolantsVariables(int i)
{
  CExpression itp = interpolants[i];
  std::unordered_set<signed> variables_itp;

  if (!m.isTrue(itp) && !m.isNil(itp))
  {
    std::string t = m.toString(itp);
    // std::cout << "c " << t.c_str() << std::endl;
    extractVariablesITP(itp, variables_itp);
  }
  return variables_itp;
}

void DeSAT::extractVariablesITP(CExpression itp, std::unordered_set<signed> &vars)
{
  if (m.isLiteral(itp))
  {
    signed l = std::abs(m.getLiteral(itp));
    vars.insert(l);
    return;
  }
  size_t sz = m.nChildren(itp);
  for (int i = 0; i < sz; i++)
  {
    CExpression c = m.getChild(itp, i);
    extractVariablesITP(c, vars);
  }
}

bool DeSAT::importInterpolants(std::vector<std::vector<int>> &itp_cls_export, int i)
{
  clock_t before = clock();
  bool did_something = false;
  itp_cls_export.clear();

  CExpression itp = interpolants[i];
  if (m.isNil(itp))
  {
    solutions_imported++;
    return did_something;
  }

  if (!m.isTrue(itp))
  {
    // if (verbosity>0)
    {
      std::string t = m.toString(itp);
    }
    globalSolver->clearNewClauses();
    globalSolver->addConstraint(itp);
    std::vector<std::vector<signed>> itp_clauses = globalSolver->getNewClauses();
    itp_cls_export = itp_clauses; // TODO CADICAL ? est-ce necessaire

    for (unsigned i = 0; i < itp_clauses.size(); i++)
    {
      std::vector<int> cls;
      for (unsigned j = 0; j < itp_clauses[i].size(); j++)
      {
        cls.emplace_back(itp_clauses[i][j]);
      }
      itp_cls_export.emplace_back(cls);
    }

    globalSolver->addNewClauses();
    interpolants_imported++;
    did_something = true;
  }
  interpolants[i] = m.mkNil();

  clock_t after = clock();
  importTime += after - before;
  return did_something;
}

bool DeSAT::findDisagreement(void)
{
  bool res = false;

  std::vector<signed> new_trail;

  for (unsigned i = 0; i < n_partitions; i++)
  {
    for (signed v = 1; v <= (signed)maxVar; v++)
    {
      if (sharedVariables.occurs(v, i))
      {
        ModelValue mvi = partitions[i]->get(v);
        if (mvi != M_UNDEF)
        {
          for (unsigned j = i + 1; j < n_partitions; j++)
          {
            if (sharedVariables.occurs(v, j))
            {
              ModelValue mvj = partitions[j]->get(v);
              if (mvj != M_UNDEF && mvj != mvi)
              {
                signed x = (mvi == M_TRUE) ? v : -v;

                if (std::find(trail.begin(), trail.end(), -x) != trail.end())
                  throw std::runtime_error("Global assignment not respected");

                if (std::find(trail.begin(), trail.end(), x) == trail.end() &&
                    std::find(new_trail.begin(), new_trail.end(), x) == new_trail.end())
                  new_trail.push_back(x);
              }
            }
          }
        }
      }
    }
  }

  for (unsigned i = 0; i < new_trail.size(); i++)
  {
    //print("Forcing model agreement on %d\n", new_trail[i]);
    trail.push_back(new_trail[i]);
    reasons.push_back(false);
  }

  return new_trail.size() > 0;
}
