// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include "TE.h"

std::string TrivialExpression::toString(void) const
{
  std::string res;

  switch(type)
  {
  case NIL: res+="NIL"; break;
  case AND: res+="(AND"; break;
  case OR: res+="(OR"; break;
  case TRUE: res+="TRUE"; break;
  case FALSE: res+="FALSE"; break;
  case LITERAL:
  {
    std::stringstream s("");
    s << value;
    res+=s.str(); break;
  }
  default:
    throw std::runtime_error("Unexpected TrivialExpression type.");
  }

  if (type==AND || type==OR)
  {
    for (unsigned i=0; i<children.size(); i++)
    {
      res+=" " + children[i]->toString();
    }

    res+=")";
  }

  return res;
}

//void TrivialExpression::simplify(void)
//{
//  if (type==AND)
//  {
//    for (unsigned i=0; i<children.size(); i++)
//      if (children[i]->isLiteral())
//      {
//        for (unsigned j=0; j<children.size(); j++)
//          if (j!=i) children[j]->set(children[i]->value);
//      }
//  }
//  else if (type==OR)
//  {
//    for (unsigned i=0; i<children.size(); i++)
//      if (children[i]->isLiteral())
//      {
//        for (unsigned j=0; j<children.size(); j++)
//          if (j!=i) children[j]->set(-children[i]->value);
//      }
//  }
//
//  for (unsigned i=0; i<children.size(); i++)
//    children[i]->simplify();
//
//  switch (type)
//  {
//  case AND:
//    {
//      unsigned shift=0;
//      for (unsigned i=0; i<children.size(); i++)
//      {
//        TrivialExpression c = children[i];
//        if (c.isFalse())
//        {
//          type = FALSE;
//          children.clear();
//          goto ANDdone;
//        }
//        else if (c.isTrue())
//        {
//          shift++;
//        }
//        else if (c.isAnd())
//        {
//          for (unsigned i=0; i<c.children.size(); i++)
//            children.push_back(c.children[i]);
//          shift++;
//        }
//        else if (shift!=0)
//          children[i-shift] = c;
//      }
//      children.resize(children.size()-shift);
//      if (children.size()==0)
//        type = TRUE;
//      else if (children.size()==1)
//      {
//        TrivialExpression c0 = children[0];
//        type = c0.type;
//        value = c0.value;
//        children.clear();
//        for (unsigned i=0; i<c0.children.size(); i++)
//          children.push_back(c0.children[i]);
//      }
//      ANDdone: ;
//    }
//    break;
//  case OR:
//    {
//      unsigned shift=0;
//      for (unsigned i=0; i<children.size(); i++)
//      {
//        TrivialExpression c = children[i];
//        if (c.isTrue())
//        {
//          type = TRUE;
//          children.clear();
//          goto ORdone;
//        }
//        else if (c.isFalse())
//        {
//          shift++;
//        }
//        else if (c.isOr())
//        {
//          for (unsigned i=0; i<c.children.size(); i++)
//            children.push_back(c.children[i]);
//          shift++;
//        }
//        else if (shift!=0)
//          children[i-shift] = c;
//      }
//      children.resize(children.size()-shift);
//      if (children.size()==0)
//        type = FALSE;
//      else if (children.size()==1)
//      {
//        TrivialExpression c0 = children[0];
//        type = c0.type;
//        value = c0.value;
//        children.clear();
//        for (unsigned i=0; i<c0.children.size(); i++)
//          children.push_back(c0.children[i]);
//      }
//      ORdone: ;
//    }
//    break;
//  case LITERAL:
//  case TRUE:
//  case FALSE:
//  case NIL:
//  default:
//    /* Nothing */ ;
//  }
//}

bool TrivialExpression::set(signed l)
{
  if (type==LITERAL)
  {
    if (value==l)
    {
      type = TRUE;
      return true;
    }
    else if (value==-l)
    {
      type = FALSE;
      return true;
    }

    return false;
  }
  else
    assert(false);
  /*{
    bool did_something=false;
    for (unsigned i=0; i<children.size(); i++)
      if (children[i].set(l))
        did_something=true;
    return did_something;
  }*/
  return false;
}

TrivialExpression TrivialExpression::duplicate()
{
  assert(false);
  return NilTE();
}
