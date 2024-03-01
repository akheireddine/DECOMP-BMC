// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <iostream>
#include <sstream>

#include "interpolator_m.h"
#include "interpolator_p.h"
#include "interpolator_im.h"

#include "partition.h"

using namespace Desat;


Partition::Partition(ExpressionManager &em, SharedVariables &sharedVariables, unsigned id, unsigned verbosity, bool proof) :
  m(em),
  SATSolver(em),
  sharedVariables(sharedVariables),
  id(id),
  interpolator(NULL)
{
  solver = new Cadical(em, proof);
	sharedVariables.add(&v);
}


Partition::~Partition(void)
{
  if (solver) delete solver;
  if (interpolator)
    delete interpolator;
}

void Partition::setInterpolationMode(InterpolationMode interpolationMode)
{
  if (interpolator)
    delete interpolator;

  if (interpolationMode==NONE)
    interpolator = NULL;
  if (interpolationMode==MCMILLAN)
    interpolator = new InterpolatorM(m, sharedVariables);
  else if (interpolationMode==INVERSE_MCMILLAN)
    interpolator  = new InterpolatorIM(m, sharedVariables);
  else if (interpolationMode==PUDLAK)
    interpolator  = new InterpolatorP(m, sharedVariables);
  else
    throw std::runtime_error("Unknown interpolation mode.");

  solver->setInterpolator(interpolator);
}

void Partition::setInterpolationMode(InterpolationMode interpolationMode, SharedVariables &masterSharedVariables)
{
  if (interpolator)
    delete interpolator;

  if (interpolationMode==NONE)
    interpolator = NULL;
  if (interpolationMode==MCMILLAN)
    interpolator = new InterpolatorM(m, masterSharedVariables);
  else if (interpolationMode==INVERSE_MCMILLAN)
    interpolator  = new InterpolatorIM(m, masterSharedVariables);
  else if (interpolationMode==PUDLAK)
    interpolator  = new InterpolatorP(m, masterSharedVariables);
  else
    throw std::runtime_error("Unknown interpolation mode.");

	  solver->setInterpolator(interpolator);
}

void Partition::setInterpolationMode(ExpressionManager &em, InterpolationMode interpolationMode, SharedVariables &masterSharedVariables)
{
  if (interpolator)
    delete interpolator;

  if (interpolationMode==NONE)
    interpolator = NULL;
  if (interpolationMode==MCMILLAN)
    interpolator = new InterpolatorM(em, masterSharedVariables);
  else if (interpolationMode==INVERSE_MCMILLAN)
    interpolator  = new InterpolatorIM(em, masterSharedVariables);
  else if (interpolationMode==PUDLAK)
    interpolator  = new InterpolatorP(em, masterSharedVariables);
  else
    throw std::runtime_error("Unknown interpolation mode.");

	 solver->setInterpolator(interpolator);
}


void Partition::setInterrupt()
{
  solver->setInterrupt();
}

void Partition::unsetInterrupt()
{
  solver->unsetInterrupt();
}

bool Partition::addClause(const std::vector<signed> &literals)
{
  for (unsigned i=0; i<literals.size(); i++)
  {
    unsigned var = (literals[i]<0) ? -literals[i] : literals[i];
    v.setOccurs(var);
    while (var>=numVars()) addVar();
  }

  return solver->addClause(literals);
}

void Partition::setPhase(const int var, const bool phase){
  solver->setPhase(var, phase);
}


Expression Partition::getInterpolant(const std::vector<signed> &assumptions)
{
  Expression res;

  if (interpolator==NULL)
  {
    bool r=solver->solve(assumptions);

    std::vector<int> cf_tmp;
    ((Cadical *)solver)->getConflict(cf_tmp);

    if (r)
      res = m.mkTrue();
    else
    {
      if (cf_tmp.size()==0)
        res = m.mkFalse();
      else if (cf_tmp.size()==1)
        //res = m.mkLiteral(-cf_tmp[0]);
		res = m.mkLiteral(cf_tmp[0]);
      else
      {
		//res = m.mkOr(m.mkLiteral(-cf_tmp[0]), m.mkLiteral(-cf_tmp[1]));
		res = m.mkOr(m.mkLiteral(cf_tmp[0]), m.mkLiteral(cf_tmp[1]));
        for (unsigned i=2; i<cf_tmp.size(); i++)
		{
          //res = m.mkOr(res, m.mkLiteral(-cf_tmp[i]));
		  res = m.mkOr(res, m.mkLiteral(cf_tmp[i]));
		}

	  }
	}
  }
  else
  {
	  res = solver->getInterpolant(assumptions);
  }

  if (verbosity>1)
  {
    std::string t = m.toString(res);
	std::cout << id << ": ITP: " << t << std::endl;
    print("%02d: ITP: %s\n", id, t.c_str());
    if (m.isTrue(res))
    {
      Expression q = solver->getModel();
      std::cout << "Model #" << id << ": " << m.toString(q) << std::endl;
	  //solutions_imported++;
    }
  }

#ifdef _DEBUG
  if (m.isTrue(res))
  {
    for (unsigned i=0; i<assumptions.size(); i++)
    {
      if (solver->get(assumptions[i]) == M_FALSE)
      {
        std::stringstream ss;
        ss << "Assumption " << assumptions[i] << " not satisfied";
        throw std::runtime_error(ss.str().c_str());
      }
    }
  }
#endif

  return res;
}
