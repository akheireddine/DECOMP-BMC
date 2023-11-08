// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011

#include <iostream>
#include <sstream>
#include <fstream>

#include "RBCManager.h"

RBCManager::RBCManager(void)
{
	rbcs.push_back(RBC(R_NIL, 0, 0, false, false));
	rbcs.push_back(RBC(R_TRUE, 0, 0, true, true));
	rbcs.push_back(RBC(R_FALSE, 0, 0, true, true));
}

RBCManager::~RBCManager(void)
{
}

signed RBCManager::mkLiteral(signed l)
{
	bool sgn = (l<0);
	unsigned v = (sgn) ? -l : l;

	if (sgn)
	{
		if (nliterals.size() <= v)
		{
			size_t old_size = nliterals.size();
			nliterals.resize(v+1);

			for (unsigned i=old_size; i<nliterals.size(); i++)
				nliterals[i] = newRBC(R_LITERAL, -(signed)i);
		}
	}
	else
	{
		if (literals.size() <= v)
		{
			size_t old_size = literals.size();
			literals.resize(v+1);

			for (unsigned i=old_size; i<literals.size(); i++)
				literals[i] = newRBC(R_LITERAL, i);
		}
	}

	return (sgn) ? nliterals[v] : literals[v];
}

signed RBCManager::mkAnd(signed a, signed b, bool swap)
{
	if (a==b) return a;
	else if(a==-b) return mkFalse();
	else if( isTrue(a) ) return b;
	else if( isTrue(b) ) return a;
	else if( isFalse(a) ) return mkFalse();
	else if( isFalse(b) ) return mkFalse();
	else if ( isLiteral(a) && isLiteral(b) && getLiteral(a)==-getLiteral(b)) return mkFalse();
	else
	{
		signed res;
		signed l = a, r = b;
		if (swap && b<a) std::swap(l, r);

		signed q = ands.get(l, r);
		if (q!=0)
			res = q;
		else
		{
			signed x = newRBC(R_AND, l, r);
			ands.put(l, r, x);
			res = x;
		}

		//std::cout << "RES = " << res << " (" << get(res).left << "," << get(res).right << ")" << std::endl;
		return res;
	}
}

signed RBCManager::mkOr(signed a, signed b, bool optimize)
{
	if (a==b) return a;
	else if(a==-b) return mkTrue();
	else if( isTrue(a) ) return mkTrue();
	else if( isTrue(b) ) return mkTrue();
	else if( isFalse(a) ) return b;
	else if( isFalse(b) ) return a;
	else if ( isLiteral(a) && isLiteral(b) && getLiteral(a)==-getLiteral(b)) return mkTrue();
	else
	{
		signed res;

		if (optimize && isClause(a) && isClause(b)) {
			getLiterals(a, t1);
			getLiterals(b, t2);
			t3.clear();
			for (unsigned i=0; i<t1.size(); i++)
				t3.push_back(mkLiteral(t1[i]));
			for (unsigned i=0; i<t2.size(); i++) {
				signed q = mkLiteral(t2[i]);
				if (std::find(t3.begin(), t3.end(), q) == t3.end())
					t3.push_back(q);
			}
			assert(t3.size() > 0);
			res = t3[0];
			for (unsigned i=1; i<t3.size(); i++)
				res = mkOr(res, t3[i], false);
		}
		else {



			signed l=a, r=b;
			if (b<a) std::swap(l, r);

			// CMW: Unsure whether this is helpful.

			if ( isOr(l) && isOr(r) )
			{
				if ( getChild(l, 0) == getChild(r, 0) ) {
					r = getChild(r, 1);
					if (r<l) std::swap(l, r);
				}
				else if ( getChild(l, 1) == getChild(r, 1) ) {
					r = getChild(r, 0);
					if (r<l) std::swap(l, r);
				}
				else if ( getChild(l, 0) == getChild(r, 1) ) {
					r = getChild(r, 0);
					if (r<l) std::swap(l, r);
				}
				else if ( getChild(l, 1) == getChild(r, 0) ) {
					r = getChild(r, 1);
					if (r<l) std::swap(l, r);
				}
			}

			signed q = ors.get(l, r);
			if (q!=0)
				res = q;
			else
			{
				signed x = newRBC(R_OR, l, r);
				ors.put(l, r, x);
				res = x;
			}
		}
		//std::cout << "RES = " << res << " (" << get(res).left << "," << get(res).right << ")" << std::endl;
		return res;
	}
}

signed RBCManager::mkEq(signed a, signed b)
{
	throw std::runtime_error("NYI: RBC-EQ");
}

signed RBCManager::newRBC(RBCOperator op, signed left, signed right)
{
	//std::cout << "NEW: " << op << " " << left << " " << right << std::endl;

	rbcs.push_back(RBC());
	rbcs.back().op = op;
	rbcs.back().left = left;
	rbcs.back().right = right;

	rbcs.back().isClause =
		left>0 && right>0 && op==R_OR &&
		(rbcs[left].isClause && rbcs[right].isClause);

	rbcs.back().isCube =
		left>0 && right>0 && op==R_AND &&
		(rbcs[left].isCube && rbcs[right].isCube);

	return rbcs.size() - 1;
}

signed RBCManager::newRBC(RBCOperator op, signed l)
{
	rbcs.push_back(RBC());
	rbcs.back().op = op;
	rbcs.back().left =  l;
	// rbcs.back().right = NULL;
	rbcs.back().isClause = true;
	rbcs.back().isCube = true;
	return rbcs.size() - 1;
}

std::string RBCManager::toString(signed x) const
{
	std::string res;
	if (s(x)) res += "(NOT ";

	if (isNil(x))
		res += "NIL";
	else if (isTrue(x))
		res += "TRUE";
	else if (isFalse(x))
		res += "FALSE";
	else if (isLiteral(x))
	{
		std::stringstream s("");
		s << getLiteral(x);
		res += s.str();
	}
	else if (isClause(x))
	{

		res += "(OR";

		std::vector<signed> lits;
		getLiterals(x, lits);

		for (unsigned i=0; i<lits.size(); i++)
		{
			std::stringstream s("");
			s << lits[i];
			res += " " + s.str();
		}

		res += ")";

	}
	else if (isCube(x))
	{
		res += "(AND";

		std::vector<signed> lits;
		getLiterals(x, lits);

		for (unsigned i=0; i<lits.size(); i++)
		{
			std::stringstream s("");
			s << lits[i];
			res += " " + s.str();
		}

		res += ")";
	}
	else
	{
		switch (get(x).op)
		{
		case R_NIL: res += "NIL"; break;
		case R_TRUE: res += "TRUE"; break;
		case R_FALSE: res += "FALSE"; break;
		case R_AND:
			res += "(AND";
			res += " " + toString(get(x).left);
			res += " " + toString(get(x).right);
			res += ")";
			break;
		case R_OR:
			res += "(OR";
			res += " " + toString(get(x).left);
			res += " " + toString(get(x).right);
			res += ")";
			break;
		case R_EQ:
			res += "(=";
			res += " " + toString(get(x).left);
			res += " " + toString(get(x).right);
			res += ")";
			break;
		case R_LITERAL:
			{
				std::stringstream s("");
				s << getLiteral(x);
				res+=s.str();
				break;
			}
		default:
			throw std::runtime_error("Unexpected RBC type.");
		}
	}

	if (s(x)) res += ")";

	return res;
}

signed RBCManager::duplicate(signed x, const RBCManager &other)
{
	tStack.clear();
	tStack.push_back(x);
	cache.clear();


	while(!tStack.empty())
	{
		int q = tStack.back();

		//std::cout << "STACK SIZE: " << tStack.size() << std::endl;
		//std::cout << "NOW: " << q << std::endl;

		Cache::const_iterator it = cache.find(v(q));
		if (it!=cache.end())
			tStack.pop_back();
		else
		{
			if (other.isNil(q))
			{
				cache[v(q)] = mkNil();
				tStack.pop_back();
			}
			else if (other.isTrue(q))
			{
				cache[v(q)] = mkTrue();
				tStack.pop_back();
			}
			else if (other.isFalse(q))
			{
				cache[v(q)] = mkFalse();
				tStack.pop_back();
			}
			else if (other.isLiteral(q))
			{
				cache[v(q)] = mkLiteral(other.getLiteral(q));
				tStack.pop_back();
			}
			else if (other.isAnd(q) || other.isOr(q) || other.isEq(q))
			{
				bool children_done = true;
				int lChild = other.getChild(q, 0);
				int rChild = other.getChild(q, 1);

				//std::cout << "CHILDREN: " << lChild << ", " << rChild << std::endl;

				if (lChild == q || rChild == q)
					throw std::runtime_error("self-referential expression");

				if (cache.find(v(lChild))==cache.end())
				{
					tStack.push_back(lChild);
					children_done = false;
				}
				if (cache.find(v(rChild))==cache.end())
				{
					tStack.push_back(rChild);
					children_done = false;
				}

				if (children_done)
				{
					int c1 = cache[v(lChild)];
					int c2 = cache[v(rChild)];

					if (s(lChild)) c1 = mkNeg(c1);
					if (s(rChild)) c2 = mkNeg(c2);

					//std::cout << "CHILDREN: " << toString(c1) << ", " << toString(c2) << std::endl;

					if (other.isAnd(q)){
						cache[v(q)] = mkAnd(c1, c2);
					}
					else if (other.isOr(q)){
						cache[v(q)] = mkOr(c1, c2);
					}
					else if (other.isEq(q))
						cache[v(q)] = mkEq(c1, c2);

					tStack.pop_back();
				}
			}
			else
				throw std::runtime_error("unexpected expression type.");
		}
	}

	signed r = cache[x];
	return (s(x)) ? mkNeg(r) : r;
}


void RBCManager::fromITPtoBuffer(std::vector<signed> &buf, signed x)
{
	typedef enum { OP_UNDEF=0, OP_AND, OP_OR, OP_NOT, OP_LIT, OP_NIL, OP_FALSE, OP_TRUE } operations;

	std::vector<signed> nodeId;
	std::vector< std::vector< signed> > childrenId;
	std::vector<signed> nodeType;

	seen.clear();

	std::vector<signed> mpiStack;

	assert(!isNil(x));

	if(isTrue(x)) buf.push_back(OP_TRUE);
	else if(isFalse(x)) buf.push_back(OP_FALSE);
	else mpiStack.push_back(x);

	while (!mpiStack.empty()) {
		signed cur = mpiStack.back();
		mpiStack.pop_back();

		unsigned v = (cur<0)?-cur:cur;
		if (seen.find(v)!=seen.end())
			continue;

		if (cur < 0) {

			if(std::find(nodeId.begin(),nodeId.end(),-cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(-cur);
				nodeType.push_back(OP_NOT);
				childrenId.push_back(children);

			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), -cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_NOT;

			}

			mpiStack.push_back(-cur);
			continue;
		}
		else if (isLiteral(cur)) {
			signed l = getLiteral(cur);


			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;
				children.push_back(l);

				nodeId.push_back(cur);
				nodeType.push_back(OP_LIT);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF){
					nodeType[pos] = OP_LIT;

					childrenId[pos].push_back(l);
				}

			}

		}
		else if (isOr(cur)) {

			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(cur);
				nodeType.push_back(OP_OR);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_OR;

			}

		}
		else if (isAnd(cur)) {

			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(cur);
				nodeType.push_back(OP_AND);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_AND;

			}
		}


		std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
		unsigned pos_cur = std::distance(nodeId.begin(), it);

		for (unsigned i = 0 ; i < nChildren(cur); i++) {
			signed c = getChild(cur, i);
			mpiStack.push_back(c);
			if (c < 0){

				if(std::find(nodeId.begin(),nodeId.end(),-c)==nodeId.end())
				{
					std::vector<signed> children;

					nodeId.push_back(-c);
					nodeType.push_back(OP_UNDEF);
					childrenId.push_back(children);

				}

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), -c);
				unsigned pos = std::distance(nodeId.begin(), it);

				childrenId[pos_cur].push_back(pos);

			} else {

				if(std::find(nodeId.begin(),nodeId.end(),c)==nodeId.end())
				{
					std::vector<signed> children;

					nodeId.push_back(c);
					nodeType.push_back(OP_UNDEF);
					childrenId.push_back(children);

				}

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), c);
				unsigned pos = std::distance(nodeId.begin(), it);
				childrenId[pos_cur].push_back(pos);
			}


		}

		seen.insert(v);
	}

	for(unsigned t=0; t < nodeId.size(); t++)
	{
		buf.push_back(nodeType[t]);
		assert(nodeType[t] != OP_UNDEF);
		if(nodeType[t] == OP_LIT){
			assert( childrenId[t].size() == 1);
			buf.push_back(childrenId[t][0]);
		} else {

		assert(childrenId[t].size() > 0);
		for(unsigned w=0; w < childrenId[t].size(); w++)
		{
			signed z = childrenId[t][w];
			buf.push_back(z);
		}
		}
	}

}

void RBCManager::toMPI(std::vector<signed> &buf, signed x)
{

	int OP_AND = 1;
	int OP_OR = 2;
	int OP_NOT = 3;
	int OP_LIT = 4;
	int OP_NIL = 5;
	int OP_FALSE = 6;
	int OP_TRUE = 7;
	int OP_UNDEF = 0;

	std::vector<signed> nodeId;
	std::vector< std::vector< signed> > childrenId;
	std::vector<signed> nodeType;

	seen.clear();

	std::vector<signed> mpiStack;

	assert(!isNil(x));

	if(isTrue(x)) buf.push_back(OP_TRUE);
	else if(isFalse(x)) buf.push_back(OP_FALSE);
	else mpiStack.push_back(x);

	while (!mpiStack.empty()) {
		signed cur = mpiStack.back();
		mpiStack.pop_back();

		unsigned v = (cur<0)?-cur:cur;
		if (seen.find(v)!=seen.end())
			continue;

		if (cur < 0) {

			if(std::find(nodeId.begin(),nodeId.end(),-cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(-cur);
				nodeType.push_back(OP_NOT);
				childrenId.push_back(children);

			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), -cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_NOT;

			}

			mpiStack.push_back(-cur);
			continue;
		}
		else if (isLiteral(cur)) {
			signed l = getLiteral(cur);


			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;
				children.push_back(l);

				nodeId.push_back(cur);
				nodeType.push_back(OP_LIT);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF){
					nodeType[pos] = OP_LIT;

					childrenId[pos].push_back(l);
				}

			}

		}
		else if (isOr(cur)) {

			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(cur);
				nodeType.push_back(OP_OR);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_OR;

			}

		}
		else if (isAnd(cur)) {

			if(std::find(nodeId.begin(),nodeId.end(),cur)==nodeId.end())
			{
				std::vector<signed> children;

				nodeId.push_back(cur);
				nodeType.push_back(OP_AND);
				childrenId.push_back(children);


			} else {

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
				unsigned pos = std::distance(nodeId.begin(), it);
				if(nodeType[pos] == OP_UNDEF)
					nodeType[pos] = OP_AND;

			}
		}


		std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), cur);
		unsigned pos_cur = std::distance(nodeId.begin(), it);

		for (unsigned i = 0 ; i < nChildren(cur); i++) {
			signed c = getChild(cur, i);
			mpiStack.push_back(c);
			if (c < 0){

				if(std::find(nodeId.begin(),nodeId.end(),-c)==nodeId.end())
				{
					std::vector<signed> children;

					nodeId.push_back(-c);
					nodeType.push_back(OP_UNDEF);
					childrenId.push_back(children);

				}

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), -c);
				unsigned pos = std::distance(nodeId.begin(), it);

				childrenId[pos_cur].push_back(pos);

			} else {

				if(std::find(nodeId.begin(),nodeId.end(),c)==nodeId.end())
				{
					std::vector<signed> children;

					nodeId.push_back(c);
					nodeType.push_back(OP_UNDEF);
					childrenId.push_back(children);

				}

				std::vector<signed>::iterator it = std::find(nodeId.begin(), nodeId.end(), c);
				unsigned pos = std::distance(nodeId.begin(), it);
				childrenId[pos_cur].push_back(pos);
			}


		}

		seen.insert(v);
	}

	for(unsigned t=0; t < nodeId.size(); t++)
	{
		buf.push_back(nodeType[t]);
		assert(nodeType[t] != OP_UNDEF);
		if(nodeType[t] == OP_LIT){
			assert( childrenId[t].size() == 1);
			buf.push_back(childrenId[t][0]);
		} else {

		assert(childrenId[t].size() > 0);
		for(unsigned w=0; w < childrenId[t].size(); w++)
		{
			signed z = childrenId[t][w];
			buf.push_back(z);
		}
		}
	}

}

bool RBCManager::sizeRBC(signed x, unsigned limit)
{
	bool overflow_problem = false;
	seen.clear();

	std::vector<signed> dotStack;
	dotStack.push_back(x);

	unsigned num_edges =0;
	unsigned num_vertex = 0;

	while (!dotStack.empty()) {
		signed cur = dotStack.back();
		dotStack.pop_back();

		unsigned v = (cur<0)?-cur:cur;
		if (seen.find(v)!=seen.end())
			continue;

		if (cur < 0) {
			num_edges++;
			dotStack.push_back(-cur);
			num_vertex++;
			continue;
		}
		else if (isLiteral(cur)) {
			signed l = getLiteral(cur);
			num_vertex++;
		}
		else if (isOr(cur)) {
			num_vertex++;
		}
		else if (isAnd(cur)) {
			num_vertex++;
		}

		for (unsigned i = 0 ; i < nChildren(cur); i++) {
			signed c = getChild(cur, i);
			dotStack.push_back(c);
			if (c < 0)
				num_edges++;
			else
				num_edges++;
		}

		seen.insert(v);
	}

	//std::cout << " num_vertex " << num_vertex << " num_edges " << num_edges << std::endl;

	if(num_vertex > limit)
		overflow_problem = true;

	return overflow_problem;

}

void RBCManager::toDot(signed x, std::string filename)
{
	std::cout << "filename : " << filename << std::endl;
	std::ofstream f(filename.c_str(), std::ios::trunc);

	f << "digraph E {" << std::endl;

	seen.clear();

	std::vector<signed> dotStack;
	dotStack.push_back(x);

	while (!dotStack.empty()) {
		signed cur = dotStack.back();
		dotStack.pop_back();

		unsigned v = (cur<0)?-cur:cur;
		if (seen.find(v)!=seen.end())
			continue;

		if (cur < 0) {
			//std::cout << "nm" << -cur << " [label=\"NOT\"];" << std::endl;
			f << "nm" << -cur << " [label=\"NOT\"];" << std::endl;
			dotStack.push_back(-cur);
			f << "nm" << -cur << " -> n" << -cur << std::endl;
			continue;
		}
		else if (isLiteral(cur)) {
			signed l = getLiteral(cur);
			//std::cout << "n" << cur << " [label=\"L " << l << "\"];" << std::endl;
			f << "n" << cur << " [label=\"L " << l << "\"];" << std::endl;
		}
		else if (isOr(cur)) {
			//std::cout << "n" << cur << " [label=\"OR\"];" << std::endl;
			f << "n" << cur << " [label=\"OR\"];" << std::endl;
		}
		else if (isAnd(cur)) {
			//std::cout << "n" << cur << " [label=\"AND\"];" << std::endl;
			f << "n" << cur << " [label=\"AND\"];" << std::endl;
		}

		for (unsigned i = 0 ; i < nChildren(cur); i++) {
			signed c = getChild(cur, i);
			dotStack.push_back(c);
			if (c < 0)
				//std::cout << "n" << cur << " -> nm" << -c << std::endl;
				f << "n" << cur << " -> nm" << -c << " [label=\"" << i << "\"];" << std::endl;
			else
				//std::cout << "n" << cur << " -> n" << c << std::endl;
				f << "n" << cur << " -> n" << c << " [label=\"" << i << "\"];" << std::endl;

		}

		seen.insert(v);
	}


	f << "} " << std::endl;

	f.close();

}
