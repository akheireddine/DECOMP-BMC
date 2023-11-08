// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <cassert>

#include "decomposition_cycle.h"

CycleDecomposition::CycleDecomposition(void) : Decomposition()
{
  last = 0;
}

CycleDecomposition::~CycleDecomposition(void)
{
}

int CycleDecomposition::where(const std::vector<signed> &clause)
{
  last = (last + 1) % partitions;
  return last;
}