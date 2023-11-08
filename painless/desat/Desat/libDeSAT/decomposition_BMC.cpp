// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <cassert>
#include <unordered_set>
#include <assert.h>
#include <algorithm>
#include <math.h>
#include <fstream>

#include "decomposition_BMC.h"
#include "decomposition_mode.h"
#include "utils/dimacs.h"


BMCDecomposition::BMCDecomposition(void) :
  Decomposition(),
  count(0)
{
}

BMCDecomposition::~BMCDecomposition(void)
{
}



int BMCDecomposition::where(const std::vector<signed> &clause)
{
    return clause.at(count++);
}
