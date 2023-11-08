// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _DECOMPOSITION_H_
#define _DECOMPOSITION_H_

#include <vector>

class Decomposition
{
public:
  Decomposition(void) : partitions(0), variables(0), clauses(0) {}
  virtual ~Decomposition(void) {}

  virtual void setPartitions(unsigned n) { partitions=n; };
  virtual void setVariableMax(unsigned n) { variables=n; };
  virtual void setClauseMax(unsigned n) { clauses=n; };

  virtual int where(const std::vector<signed> &clause) = 0;

protected:
  unsigned partitions;
  unsigned variables;
  unsigned clauses;
};

#endif
