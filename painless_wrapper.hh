
#pragma once

#include "painless/painless-src/sharing/Sharer.h"
#include "painless/painless-src/working/WorkingStrategy.h"

/// Is it the end of the search
extern atomic<bool> globalEnding;

/// Final result
extern SatResult finalResult;

/// Model for SAT instances
extern vector<int> finalModel;

/// Array of sharers
extern Sharer **sharers;

/// Size of the array of sharers
extern int nSharers;

// -------------------------------------------
// Declaration of global variables
// -------------------------------------------
atomic<bool> globalEnding(false);

SatResult finalResult = UNKNOWN;

vector<int> finalModel;

Sharer **sharers = NULL;

int nSharers = 0;

// -------------------------------------------
