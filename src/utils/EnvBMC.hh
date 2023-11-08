
#pragma once

#include <vector>
#include <string>
#include <unordered_set>

#include "../../painless/painless-src/utils/Parameters.h"

typedef struct variableBelonging
{
  int num_step = -2;
  int num_partition = -1;
  bool in_property = false;
} VB;

class EnvBMC
{
public:
  EnvBMC()
  {
    decompstrat = Parameters::getParam("decomp", "batch");
    cpus = Parameters::getIntParam("c", 1); // std::thread::hardware_concurrency());
    timeout = Parameters::getIntParam("t", -1);
    memory = Parameters::getIntParam("max-memory", -1) * 1024 * 1024;
    k =  0;
    solver_type = Parameters::getParam("s", "desat");
    nb_leafs = Parameters::getIntParam("nleafs", 0);
    nb_grp_steps = Parameters::getIntParam("nsteps", 0);
    nb_clauses = 0, nb_variables = 0;
    log = Parameters::isSet("log");
    log_itp_only = Parameters::isSet("log-itp");
    shrValue = Parameters::getIntParam("shr-strat",0);
    stop_first_unsat = Parameters::isSet("stop-first");
    if (log || log_itp_only)
      logFile.open(Parameters::getFilenameStr() + "_aug_itp_"+decompstrat+"_"+std::to_string(nb_leafs)+".dimacs");
  }

  ~EnvBMC() {}

  /**************************** BSaLTic main parameters ****************************/
  std::unordered_set<std::string> ltl_varnames;
  std::vector<VB> info_variables;
  std::vector<int> rename_loop;
  std::vector<int> rename_step;
  std::vector<std::vector<std::vector<int> > > clauses_partition;
  std::vector<std::vector<int> > clauses_g;
  int idmax_model;

  int nb_clauses, nb_variables;
  std::string decompstrat;
  bool log, log_itp_only;
  bool stop_first_unsat;
  int shrValue;
  std::ofstream logFile;
  int cpus;
  int timeout;
  int memory;
  int k;
  int nb_leafs, nb_grp_steps;
  std::string solver_type;
};
