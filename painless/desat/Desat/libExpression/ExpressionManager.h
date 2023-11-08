// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _EXPRESSIONMANAGER_H_
#define _EXPRESSIONMANAGER_H_

#include "expression.h"

#if defined(_USE_TRIVIAL_EXPRESSIONS)

#include "TEManager.h"
typedef TrivialExpressionManager ExpressionManager;

#elif defined(_USE_RBC_EXPRESSIONS)

// #include <RBCManager.h>
#include "../libRBC/RBCManager.h"
typedef RBCManager ExpressionManager;

#else
#error UNKNOWN EXPRESSION TYPE SELECTED
#endif

#endif
