// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _TRIVIALEXPRESSION_H_
#define _TRIVIALEXPRESSION_H_

#include <cassert>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

class TrivialExpression
{
public:
  typedef enum { AND, OR, LITERAL, TRUE, FALSE, NIL } E_TYPE;
  
  TrivialExpression(void) : type(NIL), value(0) {}
  TrivialExpression(E_TYPE t) : type(t), value(0) {}
  TrivialExpression(bool v) : value(0) { if (v) type=TRUE; else type=FALSE; }
  TrivialExpression(signed x) { type=LITERAL; value=x; }
  TrivialExpression(const TrivialExpression &other) { (*this) = other; }

  ~TrivialExpression(void) 
  {
  }
  
  bool isNil(void) const { return type==NIL; }
  bool isTrue(void) const { return type==TRUE; }
  bool isFalse(void) const { return type==FALSE; }
  bool isLiteral(void) const { return type==LITERAL; }
  bool isAnd(void) const { return type==AND; }
  bool isOr(void) const { return type==OR; }
  bool isEq(void) const { assert(false); return false; }

  unsigned nChildren(void) const { return children.size(); }
  const TrivialExpression *getChild(unsigned i) { return children[i]; }
  const TrivialExpression *getChild(unsigned i) const 
  {
    assert(i<children.size());
    return children[i]; 
  }
  
  signed getLiteral(void) const { return value; }

  
  E_TYPE type;
  signed value; // for literals
  std::vector<const TrivialExpression*> children;

  bool isCNF(void) const
  {
    if (type!=AND)
      return false;

    size_t sz = children.size();
    for (size_t i=0; i<sz; i++)
    if (!children[i]->isClause())
      return false;
  
    return true;
  }

  bool isClause(void) const
  {
    if (type==LITERAL)
      return true;

    if (type!=OR)
      return false;

    size_t sz = children.size();
    for (size_t i=0; i<sz; i++)
    {
      if (!children[i]->isLiteral())
        return false;
    }

    return true;
  }

  bool isCube(void) const
  {
    if (type==LITERAL)
      return true;

    if (type!=AND)
      return false;

    size_t sz = children.size();
    for (size_t i=0; i<sz; i++)
    {
      if (!children[i]->isLiteral())
        return false;
    }

    return true;
  }

  TrivialExpression &operator=(const TrivialExpression &other)
  {
    type = other.type;
    value = other.value;
    children = other.children;    
    return *this;
  }

  friend bool operator<(const TrivialExpression &a, const TrivialExpression &b);
  friend bool operator==(const TrivialExpression &a, const TrivialExpression &b);

  std::string toString(void) const;
  void toBuffer(std::vector<signed> &buffer) const;
  void fromBuffer(std::vector<signed> &buffer, unsigned inx);

  void addChild(const TrivialExpression *c)
  {
    if (find(children.begin(), children.end(), c)==children.end())
      children.push_back(c);
  }

  void addChildrenOf(const TrivialExpression *e)
  {
    for (unsigned i=0; i<e->children.size(); i++)
      addChild(e->children[i]);
  }

  //void negate(void)
  //{
  //  switch(type)
  //  {
  //  case TRUE: type = FALSE; break;
  //  case FALSE: type = TRUE; break;
  //  case LITERAL: value = -value; break;
  //  case AND: 
  //    type = OR; 
  //    for (unsigned i=0; i<children.size(); i++)
  //      children[i].negate();
  //    break;
  //  case OR:
  //    type = AND; 
  //    for (unsigned i=0; i<children.size(); i++)
  //      children[i].negate();
  //    break;
  //  default:      
  //    throw std::exception("Unexpected TrivialExpression type.");
  //  }
  //}  

  void simplify(void);
  bool set(signed l);
  TrivialExpression duplicate();
};

class AndTE : public TrivialExpression
{
public:
  AndTE(const std::vector<TrivialExpression*> &c) : TrivialExpression(AND) 
  {    
    for (unsigned i=0; i<c.size(); i++)    
    {
      if (c[i]->isAnd())
        addChildrenOf(c[i]);
      else
        addChild(c[i]);
    }
  }

  AndTE(const TrivialExpression *a, const TrivialExpression *b) : TrivialExpression(AND) 
  { 
    if (a->isFalse() || b->isFalse())
      type = FALSE;    
    else if (a->isTrue())
    {
      type = b->type;
      value = b->value;      
      addChildrenOf(b);
    }
    else if (b->isTrue())
    {      
      type = a->type;
      value = a->value;
      addChildrenOf(a);
    }
    else if (a==b)
    {
      type = a->type;
      value = a->value;
      addChildrenOf(a);
    }
    else    
    {
      assert(children.size()==0);

      if (a->type==AND)
        addChildrenOf(a);
      else
        addChild(a);

      if (b->type==AND)
        addChildrenOf(b);
      else
        addChild(b);
    }
  }
};

class OrTE : public TrivialExpression
{
public:
  OrTE(const std::vector<TrivialExpression*> &c) : TrivialExpression(OR) 
  {      
    for (unsigned i=0; i<c.size(); i++)
    {
      if (c[i]->isOr())
        addChildrenOf(c[i]);
      else
        addChild(c[i]);
    }
  }

  OrTE(const TrivialExpression *a, const TrivialExpression *b) : TrivialExpression(OR) 
  { 
    if (a->isTrue() || b->isTrue())
      type = TRUE;
    else if (a->isFalse())
    {
      type = b->type;
      value = b->value;
      addChildrenOf(b);
    }
    else if (b->isFalse())
    {
      type = a->type;
      value = a->value;
      addChildrenOf(a);
    }
    else if (a==b)
    {
      type = a->type;
      value = a->value;
      addChildrenOf(a);
    }
    else    
    {
      assert(children.size()==0);

      if (a->type==OR)
        addChildrenOf(a);
      else
        addChild(a);

      if (b->type==OR)
        addChildrenOf(b);
      else
        addChild(b);
    }
  }
};

class LiteralTE : public TrivialExpression
{
public:
  LiteralTE(signed x) : TrivialExpression(x) { }
};

class NilTE : public TrivialExpression
{
public:
  NilTE(void) : TrivialExpression(NIL) {}
};

class FalseTE : public TrivialExpression
{
public:
  FalseTE(void) : TrivialExpression(FALSE) {}
};

class TrueTE : public TrivialExpression
{
public:
  TrueTE(void) : TrivialExpression(TRUE) {}
};


// Inlines

inline bool operator==(const TrivialExpression &a, const TrivialExpression &b)
{
  return &a == &b || 
         (a.type == b.type &&
          a.children == b.children &&
          (a.type != TrivialExpression::LITERAL || a.value == b.value));
}

inline bool operator<(const TrivialExpression &a, const TrivialExpression &b)
{
  return &a < &b;
}

#endif

