// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_BMC_H_
#define _DECOMPOSITION_BMC_H_

#include "decomposition.h"


class BMCDecomposition : public Decomposition
{
public:
  BMCDecomposition(void);
  virtual ~BMCDecomposition(void);

  virtual int where(const std::vector<signed> &clause);

protected:
  unsigned count;
};

#endif