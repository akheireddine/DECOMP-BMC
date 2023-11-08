// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _SHARED_VARIABLES_H_
#define _SHARED_VARIABLES_H_

#include <vector>

class VariableOccurrence : public std::vector<bool>
{
public:
  bool occurs(signed x) const 
  {
    unsigned v = (x<0)?-x:x;
    if (v>=size()) return false;
    return (*this)[ v ];
  }

  void setOccurs(signed x) 
  { 
    unsigned v = (x<0)?-x:x;
    if (v>=size()) resize(v+1, false);
    (*this)[ v ] = true; 
  }
};

class SharedVariables
{
public:
  bool isShared(signed x) const 
  {
    unsigned v = (x<0) ? -x : x ;
    if (v>=_isShared.size()) return false;
    return _isShared[ v ];
  }
  
  bool isShared(signed x, unsigned w1, unsigned w2) const 
  {
    return occurrences[w1]->occurs(x) && occurrences[w2]->occurs(x);
  }

  void update(void)
  {
    _isShared.clear();    

    for (unsigned i=0; i<occurrences.size(); i++)
      for (unsigned j=i+1; j<occurrences.size(); j++)
      {
        for (unsigned v = 0; v<occurrences[i]->size(); v++)
          if (isShared(v, i, j))
          {
            if (v>=_isShared.size())_isShared.resize(v+1, false);
            _isShared[v] = true;
          }
      }
  }

  bool occurs(signed x, unsigned i)  { return occurrences[i]->occurs(x); }

  void add(const VariableOccurrence* o) { occurrences.push_back(o); }

  void update(int * vec){
	 unsigned size = vec[0];
	 _isShared.clear();
	 _isShared.resize(size+1,false);
	 for(unsigned i = 1; i < size; i++){
		 if(vec[i]) _isShared[i] = true;
	 }
  }

  void setShared(signed x, bool value)
  {
	unsigned v = (x < 0) ? -x : x;
	if (v>=_isShared.size()) _isShared.resize(v+1, false);

	_isShared[v] = value;
  }

  void cleanShared(void){ _isShared.clear(); }

  void cleanOccurrences(void) { occurrences.clear(); }

protected:
  std::vector<bool> _isShared;
  std::vector<const VariableOccurrence*> occurrences;
};

#endif