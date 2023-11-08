// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <cassert>
#include <iostream>
#include "decomposition_batch.h"

BatchDecomposition::BatchDecomposition(void) :
  Decomposition(),
  count(0)
{
}

BatchDecomposition::~BatchDecomposition(void)
{
}

int BatchDecomposition::where(const std::vector<signed> &clause)
{
  unsigned x = (partitions * count++) / clauses;
  if (x >= partitions) 
    x = partitions - 1;
  return x;
}