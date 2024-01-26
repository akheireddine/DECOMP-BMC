
#pragma once

#include <fstream>
#undef max
#include <vector>

//#include "nusmv-config.h"
#include "cudd/util.h"

// NuSMV is C library. Here we include it as external C library.
// All functions mandatory in our tool are explicitly imported
// in this header file.
extern "C"
{

#include "cudd/util.h"
#include "nusmv/shell/cinit/cinit.h"
#include "nusmv/core/cinit/cinitInt.h"
#include "nusmv/addons_core/addonsCore.h"
#include "nusmv/addons_core/compass/compassCmd.h"
#include "nusmv/core/cinit/NuSMVEnv.h"
#include "nusmv/core/parser/parser.h"
#include "nusmv/core/utils/ErrorMgr.h"
#include "nusmv/core/utils/Logger.h"
#include "nusmv/core/prop/propPkg.h"
#include "nusmv/core/bmc/bmc.h"
#include "nusmv/core/prop/Prop_Rewriter.h"
#include "nusmv/core/utils/defs.h"
#include "nusmv/core/wff/wff.h"
#include "nusmv/core/wff/w2w/w2w.h"
#include "nusmv/core/compile/symb_table/SymbTable.h"
#include "nusmv/core/enc/be/BeEnc.h"

  void NuSMVCore_init_data();
  int CommandGoBmc(NuSMVEnv_ptr env, int argc, char **argv);
  int CommandShowProperty(NuSMVEnv_ptr env, int argc, char **argv);
  int Bmc_CommandBmcSetup(NuSMVEnv_ptr env, int argc, char **argv);
  int Bmc_CommandGenLtlSpecBmcOnePb(NuSMVEnv_ptr env, int argc, char **argv);

  int CommandAddProperty(NuSMVEnv_ptr env, int argc, char **argv);

  int Parser_read_model(NuSMVEnv_ptr env, char *ifile);
};
