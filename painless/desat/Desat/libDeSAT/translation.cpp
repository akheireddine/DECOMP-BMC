// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <cassert>

#include "translation.h"


Translation::Translation(void)
{
}

Translation::~Translation(void)
{
}

void Translation::setMaxVariable(unsigned n)
{
  assert(_to.size()==_from.size());
}

void Translation::setLanguages(unsigned n)
{
  if (n>_to.size())
  {
    _to.resize(n);
    _from.resize(n);
  }
}

bool Translation::insert(unsigned who, signed x, signed who_x)
{ 
  unsigned x_v = (x<0) ? -x : x;
  unsigned wx_v = (who_x<0) ? -who_x : who_x;

  if (x_v >= _to[who].size())
    _to[who].resize(x_v+1, 0);
  if (wx_v >= _from[who].size())
    _from[who].resize(wx_v+1,0);

  if (_to[who][x_v]!=0)
    return false;

  _to[who][x_v] = wx_v;
  _from[who][wx_v] = x_v;

  return true;
}

void Translation::to(const std::vector<signed> &in, unsigned who, std::vector<signed> &out) const
{
  out.resize(in.size());

  for (unsigned i=0; i<in.size(); i++)
  {
    signed t=to(in[i], who);
    assert(t != 0);
    out[i] = t;
  }
}

signed Translation::to(signed l, unsigned who) const
{
  bool sgn=(l<0);
  unsigned v=(sgn)?-l:l;

  assert(v<_to[who].size());
  assert(_to[who][v] != 0);  

  signed t = _to[who][v];
  return (sgn) ? -t : t;
}

signed Translation::from(signed l, unsigned who) const
{  
  bool sgn=(l<0);
  unsigned v=(sgn)?-l:l;
  if (v>=_from[who].size()) return 0;
  signed t = _from[who][v];
  return (sgn) ? -t : t;
}

signed Translation::translate(signed l, unsigned f, unsigned t) const
{  
  signed g = from(l, f);
  if (g==0) return 0;
  return to(g, t);
}

bool Translation::hasVariable(signed l, unsigned who) const
{
  unsigned v = (l<0) ? -l : l;
  if (v>=_to[who].size()) return false;
  return _to[who][v] != 0;
}

bool Translation::isShared(signed l, unsigned x, unsigned y) const
{
  unsigned v = (l<0) ? -l : l;
  if (v>=_to[x].size() || v>=_to[y].size()) return false;
  return _to[x][v]!=0 && _to[y][v]!=0;
}

bool Translation::isExclusive(signed l, unsigned x) const
{  
  if (!hasVariable(l, x))
    return false;

  for (unsigned i=0; i<x; i++)
  {
    if (hasVariable(l, i))
      return false;
  }
  for (unsigned i=x+1; i<_to.size(); i++)
  {
    if (hasVariable(l, i))
      return false;
  }

  return true;
}


bool DirectTranslation::isShared(signed l) const
{  
  signed g = t->from(l, from);
  if (g==0) return false;
  return t->isShared(g, from, to);
}

signed DirectTranslation::operator() (signed l) const
{
  if (from==0 && to==0) return l;
  signed g = t->from(l, from);
  assert(g!=0);  
  return t->to(g, to);
}