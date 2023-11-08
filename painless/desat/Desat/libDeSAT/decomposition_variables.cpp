// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <cassert>
#include <iostream>

#include "decomposition_variables.h"

VariableDecomposition::VariableDecomposition(void) :
  Decomposition(),
  partition_size(0)
{
}

VariableDecomposition::~VariableDecomposition(void)
{
}

void VariableDecomposition::setPartitions(unsigned n) 
{ 
  partitions=n;
}

void VariableDecomposition::setVariableMax(unsigned n) 
{ 
  variables=n;
  partition_size = variables/partitions;
  std::cout << "VARIABLES: " << variables<< std::endl;
  std::cout << "PARTITION SIZE: " << partition_size << std::endl;
}

int VariableDecomposition::where(const std::vector<signed> &clause)
{
  unsigned last = -1; 

  for (unsigned i=0; i<clause.size(); i++)
  {
    signed x = clause[i];
    unsigned v = (x<0) ? -x : x;
    unsigned p = v / partition_size;
    if (i>0 && last!=p) return -1;
    last = p;
  }  

  return (last==partitions) ? partitions-1 : last;
}