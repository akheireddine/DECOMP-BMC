// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_VARIABLES_H_
#define _DECOMPOSITION_VARIABLES_H_

#include "decomposition.h"

class VariableDecomposition : public Decomposition
{
public:
  VariableDecomposition(void);
  virtual ~VariableDecomposition(void);

  virtual void setPartitions(unsigned n);
  virtual void setVariableMax(unsigned n);

  virtual int where(const std::vector<signed> &clause);    

protected:  
  unsigned partition_size;
};

#endif