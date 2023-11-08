// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <stdlib.h>
#include <cassert>

#include "decomposition_rnd.h"

RandomDecompositon::RandomDecompositon(void) : Decomposition()
{
  srand(0);
}

RandomDecompositon::~RandomDecompositon(void)
{
}

int RandomDecompositon::where(const std::vector<signed> &clause)
{
  return rand() % partitions;
}