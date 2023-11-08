// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _MODEL_H_
#define _MODEL_H_

typedef enum { M_FALSE, M_TRUE, M_UNDEF, M_CONFLICT } ModelValue;

class Model
{
public:
  Model(void) {}
  ~Model(void) {}

  ModelValue get(signed l)
  {
    bool sgn = (l<0);
    unsigned v = sgn ? -l : l;
    ModelValue vv = get(v);
    
    switch(vv)
    {
    case M_TRUE:
      vv = (sgn) ? M_FALSE : M_TRUE; break;
    case M_FALSE:
      vv = (sgn) ? M_TRUE : M_FALSE; break;
    default:
      vv = M_UNDEF;
    }

    return vv;
  }

  virtual ModelValue get(unsigned v) = 0;
};

#endif