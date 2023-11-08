// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_RND_H
#define _DECOMPOSITION_RND_H

#include "decomposition.h"

class RandomDecompositon : public Decomposition
{
public:
  RandomDecompositon(void);
  virtual ~RandomDecompositon(void);

  virtual int where(const std::vector<signed> &clause);
};

#endif