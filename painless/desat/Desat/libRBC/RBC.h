// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _RBC_H_
#define _RBC_H_

#include <cassert>

typedef enum { R_NIL, R_TRUE, R_FALSE, R_LITERAL, R_AND, R_OR, R_EQ } RBCOperator;

class RBC
{
public:
  RBC(void) {}
  RBC(RBCOperator op, signed left, signed right, bool isClause, bool isCube) : op(op),left(left),right(right),isClause(isClause),isCube(isCube) {}
  RBC(const RBC &other) : op(other.op),left(other.left),right(other.right),isClause(other.isClause),isCube(other.isCube) {}
  ~RBC(void) {}

  RBC &operator=(const RBC &other)
  {
    op = other.op;
    left = other.left;
    right = other.right;
    isClause = other.isClause;
    isCube = other.isCube;
    return *this;
  }    

  RBCOperator op;
  signed left;
  signed right;
  bool isClause, isCube;
};

#endif