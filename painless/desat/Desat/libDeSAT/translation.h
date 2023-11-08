// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _TRANSLATION_H_
#define _TRANSLATION_H_

#include <vector>

#include <ExpressionManager.h>
#include <expression.h>

class Translation
{
protected:
  class TranslationMap : public std::vector<signed> {};
  typedef std::vector<TranslationMap> Translations;
  Translations _to, _from;

public:  
  Translation(void);
  ~Translation(void);

  void setMaxVariable(unsigned n);
  void setLanguages(unsigned n);

  unsigned numLanguages(void) const { return (unsigned) _to.size(); }

  bool insert(unsigned who, signed x, signed who_x);
  
  void to(const std::vector<signed> &in, unsigned who, std::vector<signed> &out) const;
  virtual signed to(signed l, unsigned who) const;
  virtual signed from(signed l, unsigned who) const;
  signed translate(signed l, unsigned f, unsigned t) const;
  bool hasVariable(signed l, unsigned who) const;

  virtual bool isShared(signed l, unsigned x, unsigned y) const;
  virtual bool isExclusive(signed l, unsigned x) const;
};


class DirectTranslation
{
public:
  DirectTranslation(void) : t(NULL),from(0),to(0) {}
  DirectTranslation(const Translation &t, unsigned x, unsigned y) : t(&t), from(x), to(y) {}
  ~DirectTranslation(void) {}

  bool isShared(signed l) const;
  signed operator() (signed l) const;
  bool hasVariable(signed l) const;

protected:
  const Translation *t;
  unsigned from, to;
};

#endif
