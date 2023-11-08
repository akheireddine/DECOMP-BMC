// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#ifndef _RBC_MANAGER_H_
#define _RBC_MANAGER_H_

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <iostream>

#include "RBC.h"

class RBCManager
{
public:
	RBCManager(void);
	~RBCManager(void);

	//void not_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(1);}
	//void and_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(2);}
	//void or_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(3);}
	//void end_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(4);}
	//void nil_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(5);}
	//void true_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(6);}
	//void false_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(7);}
	////void literal_f(std::vector<signed> &vec){ vec.push_back(0); vec.push_back(8);}
	//void expression_f(std::vector<signed> &vec, signed q){ vec.push_back(0); vec.push_back(9); vec.push_back(q); } 

	signed mkNil(void) { return 0; }
	signed mkTrue(void) { return 1; }
	signed mkFalse(void) { return 2; }

	signed mkLiteral(signed l);
	signed mkAnd(signed a, signed b, bool swap = true);
	signed mkOr(signed a, signed b, bool optimize=true);
	signed mkEq(signed a, signed b);


	signed mkOr(const std::vector<signed> &l)
	{
		signed res;

		if (l.size()==0)
			res = mkTrue();
		else if (l.size()==1)
			res = l[0];
		else
		{
			res = mkOr(l[0], l[1]);
			for (unsigned i=2; i<l.size(); i++)
				res = mkOr(res, l[i]);
		}

		return res;
	}

	signed mkAnd(const std::vector<signed> &l)
	{
		signed res;

		if (l.size()==0)
			res = mkFalse();
		else if (l.size()==1)
			res = l[0];
		else
		{
			res = mkAnd(l[0], l[1]);
			for (unsigned i=2; i<l.size(); i++)
				res = mkAnd(res, l[i]);
		}

		return res;
	}

	signed mkNeg(signed x) { return -x; }

	bool s(signed x) const { return (x<0); }
	unsigned v(signed x) const { return s(x) ? -x : x; }  

	const RBC &get(signed x) const { return rbcs[v(x)]; }

	bool isNil(signed x) const { return get(x).op==R_NIL; }
	bool isTrue(signed x) const { return get(x).op==(s(x) ? R_FALSE : R_TRUE); }
	bool isFalse(signed x) const { return get(x).op==(s(x) ? R_TRUE : R_FALSE); }
	bool isLiteral(signed x) const { return get(x).op==R_LITERAL; }
	bool isAnd(signed x) const { return get(x).op==R_AND; }
	bool isOr(signed x) const { return get(x).op==R_OR; }
	bool isEq(signed x) const { return get(x).op==R_EQ; }
	bool isNegative(signed x) const { return x<0; }

	unsigned getSize(signed x) const 
	{
		seen.clear();
		cnfStack.clear();
		cnfStack.push_back(x);
		signed lvl = 1;
		unsigned count = 0;

		while(!cnfStack.empty())
		{
			signed q = cnfStack.back();
			cnfStack.pop_back();

			unsigned v = (q<0)?-q:q;

			//if (v<seen.size() && seen[v])
			if (seen.find(v)!=seen.end())
				continue;


			if(isAnd(q) || isOr(q) || isEq(q))
			{   
				cnfStack.push_back(get(q).left);
				cnfStack.push_back(get(q).right);
				//std::cout << "LEFT " << toString(get(q).left) << " RIGHT " << toString(get(q).right) << std::endl;
			}
			else if (isLiteral(q))
				count++;

			//if (v>=seen.size())
			//  seen.resize(v+1, false);
			//seen[v] = true;
			seen.insert(v);
		}

		return count;
	}

	bool getLiterals(signed x, std::vector<signed> &lits) const
	{
		seen.clear();
		cnfStack.clear();
		cnfStack.push_back(x);
		lits.clear();
		signed lvl = 1;

		while(!cnfStack.empty())
		{
			signed q = cnfStack.back();
			cnfStack.pop_back();

			unsigned v = (q<0)?-q:q;

			//if (v<seen.size() && seen[v])
			if (seen.find(v)!=seen.end())
				continue;


			if(isAnd(q) || isOr(q) || isEq(q))
			{   
				cnfStack.push_back(get(q).left);
				cnfStack.push_back(get(q).right);
				//std::cout << "LEFT " << toString(get(q).left) << " RIGHT " << toString(get(q).right) << std::endl;
			}
			else if (isLiteral(q))
				lits.push_back(getLiteral(q));

			//if (v>=seen.size())
			//  seen.resize(v+1, false);
			//seen[v] = true;
			seen.insert(v);
		}

		return true;
	}

	bool isCNF(signed x) const
	{
		cnfStack.clear();
		cnfStack.push_back(x);    

		while(!cnfStack.empty())
		{
			signed q = cnfStack.back();
			cnfStack.pop_back();

			if(q>0 && isAnd(q))
			{
				cnfStack.push_back(get(q).left);
				cnfStack.push_back(get(q).right);
			}
			else if (q<0 || !isClause(q))
				return false;
		}

		return true;
	}

	bool isClause(signed x) const
	{
		return x>0 && get(x).isClause;
	}

	bool isClauseAndGet(signed x, std::vector<signed> &lits) const
	{
		if (!isClause(x)) 
			return false;

		getLiterals(x, lits);
		return true;
	}

	bool isCube(signed x) const
	{
		return x>0 && get(x).isCube;
	}

	bool isCubeAndGet(signed x, std::vector<signed> &lits) const
	{
		tStack.clear();
		tStack.push_back(x);

		while(!tStack.empty())
		{
			signed q = tStack.back();
			tStack.pop_back();

			if(q>0 && isAnd(q))
			{
				tStack.push_back(get(q).left);
				tStack.push_back(get(q).right);
			}
			else if (!isLiteral(q))
				return false;
			else 
			{
				signed l = getLiteral(q);
				if (q<0) l = -l;
				if (std::find(lits.begin(), lits.end(), l)==lits.end())
					lits.push_back(l);
			}
		}

		return true;
	}

	signed getLiteral(signed x) const
	{
		assert(get(x).op==R_LITERAL);
		return s(x) ? -get(x).left : get(x).left;
	}

	unsigned nChildren(signed x) const { return (isLiteral(x) || isTrue(x) || isFalse(x) || isNil(x)) ? 0 : 2; }
	signed getChild(signed x, unsigned i) const { assert(i<2); if (i==0) return get(x).left; else return get(x).right; } 

	signed simplify(signed x) { return x; }

	std::string toString(signed x) const;

	bool sizeRBC(signed x, unsigned limit);

	signed duplicate(signed x, const RBCManager &other);

	void clear() { ands.clear(); ors.clear(); eqs.clear(); rbcs.resize(3); }

	void toDot(signed x, std::string filename);

	void toMPI(std::vector<signed> &buf, signed x);	

	void fromITPtoBuffer(std::vector<signed> &buf, signed x);

protected:
	unsigned nextId;  

	class RBCSet : public std::vector<RBC> {};
	class Index1Set : public std::vector<unsigned> {};  
	mutable std::vector<signed> tStack;
	mutable std::vector<signed> cnfStack;
	mutable std::set<signed> seen;

	RBCSet rbcs;
	Index1Set literals, nliterals;  

	signed newRBC(RBCOperator op, signed left, signed right);
	signed newRBC(RBCOperator op, signed l);



	class Index2Set : public std::map<signed, std::map<signed, signed> > 
	{
	public:
		void put(signed l, signed r, signed x)
		{
			(*this)[l][r] = x;
		}

		signed get(signed l, signed r)
		{
			return (*this)[l][r]; // implicitly 0 (i.e., nil) for nonexisting entries.
		}
	};

	Index2Set ands, ors, eqs;

	typedef std::unordered_map<int, int> Cache;
	Cache cache;

	std::vector<signed> t1, t2, t3;
};

#endif
