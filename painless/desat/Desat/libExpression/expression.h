// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _EXPRESSION_H_
#define _EXPRESSION_H_

//#define _USE_TRIVIAL_EXPRESSIONS
#define _USE_RBC_EXPRESSIONS

#if defined(_USE_TRIVIAL_EXPRESSIONS)

#include "TE.h"
typedef TrivialExpression* Expression;
typedef const TrivialExpression* CExpression;

#elif defined(_USE_RBC_EXPRESSIONS)

typedef signed Expression;
typedef const signed CExpression;


#else
#error UNKNOWN EXPRESSION TYPE SELECTED
#endif

#endif
