// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2011
// RM Martins, 2012
/*
 * This file contains common code to master and the worker of the static topologies (mpdesat_master and mpdesat_worker)
 * WARNING: due to some changes in the code it may have some memory problems
 * IN USE BUT IT IS BEING REFACTORED TO STATIC_WORKER.C		
 */

#include <mpi.h>

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <minisat1p.h>
#include <decomposition_batch.h>

#include "mpi_commands.h"
#include "mpdesat_node.h"
#include "returnvalue.h"

#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib") //added

MPDeSATNode::MPDeSATNode(ExpressionManager &m, unsigned rank, int master,  std::vector<int> worker) :
SATSolver(m),
	rank(rank),
	mpi_loadTime(0),
	mpi_solvingTime(0),
	mpi_idleTime(0),
	mpi_totalTime(0),
	rounds(0),
	interpolants_imported(0),
	solutions_imported(0),
	all_sat_iterations(0),
	mpi_verbose(0),
	mpi_stats(0),
	verbose(0),
	clauses_interpolants(0),
	variables_interpolants(0),
	print_rounds(1),
	maxSizeInterpolant(0),
	sizeInterpolant(0),
	maxVar(0),
	maxClauses(0),
	masterId(master),
	interpolants_exported(0),
	solutions_exported(0),
    interpolant_unit_clause(0),
	interpolant_binary_clause(0),
	interpolant_ternary_clause(0),
	interpolant_other_clause(0),
	interpolant_clauses(0),
	interpolant_not_clause(0),
	interpolant_size_clauses(0),
	memory_used(0),
	memory_flag(true)
{
	while(!worker.empty())
	{
		workersId.push_back(worker.back());
		worker.pop_back();
	}

	n_workers = workersId.size();
}

MPDeSATNode::~MPDeSATNode(void)
{
}

//stats methods
void MPDeSATNode::printClause(const std::vector<signed> &literals)
{
	std::cout << "rank " << rank << " (" << literals.size() << "): ";
	for(unsigned i = 0; i < literals.size(); i++)
		std::cout << literals[i] << " ";
	std::cout << std::endl;
}

void MPDeSATNode::saveMemoryUsed(double mem)
{
	if((mem / 1048576) > memory_used) memory_used = (mem / 1048576);
	if(memory_used > 1024 && memory_flag){
		std::cout << "Node " << rank << " warning memory! ( " << memory_used << " )" << std::endl;
		memory_flag = false;
	}

}

void MPDeSATNode::showTime(clock_t totalTime, const char *name, clock_t time)
{
	std::cout << "Node " << rank << " print_round: " << print_rounds << " : " << name << " Time: " << std::setprecision(3) <<
		time/(double)CLOCKS_PER_SEC << " sec" << 
		" (" << std::setw(3) << std::setfill(' ') <<
		100 * time/(double)totalTime << " %)" << std::endl;
}

void MPDeSATNode::statsSharedVariables(void)
{
	unsigned topShared = 0;
	unsigned botShared = 0;

	for (int v=1; v<=(signed)maxVar; v++)
	{
		if (sharedVariables.isShared(v))
			botShared++;

		if (rank != 0 && masterSharedVariables.isShared(v))
			topShared++;
	}

	std::cout << std::setprecision(4); 
	std::cout << "SharedVariables : Node " << rank << " [ " << 100*topShared/(double)maxVar << " | " << 100*botShared/(double)maxVar << " ] " << std::endl;
}

void MPDeSATNode::printStats(void)
{

	showTime(mpi_totalTime,"load",mpi_loadTime);
	showTime(mpi_totalTime,"idle",mpi_idleTime);
	showTime(mpi_totalTime,"solving",mpi_solvingTime);
	if(rank == 0) showTime(mpi_totalTime,"total",mpi_totalTime);

	//std::cout << "Node " << rank << " print_round " << print_rounds << ": #rounds " << rounds << std::endl;
	std::cout << "Node " << rank << " print_round " << print_rounds << ": #all_sat_iterations " << all_sat_iterations << std::endl;
	std::cout << "Node " << rank << " print_round " << print_rounds << ": #interpolants_exported " << interpolants_exported << std::endl;
	std::cout << "Node " << rank << " print_round " << print_rounds << ": #solutions_exported " << solutions_exported << std::endl;

	std::cout << "Node " << rank << " print_round " << print_rounds << ": maxBufferInterpolant " << maxSizeInterpolant << std::endl;
	if(interpolants_exported != 0) std::cout << "Node " << rank << " print_round " << print_rounds << ": avgBufferInterpolant "  << std::setprecision(3) << (double)sizeInterpolant/interpolants_exported << std::endl;

	std::cout << "Node " << rank << " interpolant_unit_clause " << interpolant_unit_clause << std::endl;
	std::cout << "Node " << rank << " interpolant_binary_clause " << interpolant_binary_clause << std::endl;
	std::cout << "Node " << rank << " interpolant_ternary_clause " << interpolant_ternary_clause << std::endl;
	std::cout << "Node " << rank << " interpolant_clauses " << interpolant_clauses << std::endl;
	if(interpolant_clauses !=0) std::cout << "Node " << rank << " interpolant_avg_size "  << std::setprecision(3) << (double)interpolant_size_clauses/interpolant_clauses << std::endl;
	std::cout << "Node " << rank << " interpolant_not_clause " << interpolant_not_clause << std::endl;

	//std::cout << "Node " << rank << " print_round " << print_rounds << ": clauses_interpolants " << clauses_interpolants << std::endl;
	//std::cout << "Node " << rank << " print_round " << print_rounds << " : variables_interpolants " << variables_interpolants << std::endl;

	std::cout << "Node " << rank << " memory_used " << memory_used << " MB " << std::endl;

}

void MPDeSATNode::interpolantStats(Expression t)
{
		
	if(m.isLiteral(t))
		{
			interpolant_unit_clause++;
			interpolant_size_clauses++;
			interpolant_clauses++;
		} else if(m.isClause(t))
		{
			unsigned s = m.getSize(t);
			if(s == 2){
				interpolant_binary_clause++;
				interpolant_size_clauses+=2;
				interpolant_clauses++;
			} else if(s == 3){
				interpolant_ternary_clause++;
				interpolant_size_clauses+=3;
				interpolant_clauses++;
			} else {
				interpolant_other_clause++;
				interpolant_size_clauses+=s;
				interpolant_clauses++;
			}
		} else interpolant_not_clause++;

}

//send methods
void MPDeSATNode::sendSharedVariables(void)
{
	sendBuffer.clear();

	sendBuffer.push_back(0);
	for(unsigned j = 1; j <= maxVar; j++)
	{
		if(sharedVariables.isShared(j)) sendBuffer.push_back(1);
		else sendBuffer.push_back(0);
	}

	sendBuffer[0] = sendBuffer.size();

	for(unsigned id = 0; id < workersId.size(); id++)
		MPI_safe_send(sendBuffer,workersId[id],MPI_DATA_SHARED_VARIABLES);

	sendBuffer.clear();
}

//receive methods
Expression MPDeSATNode::receiveInterpolant(int * buf)
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
	std::vector<Expression> nodeExpressions;

	unsigned size = buf[0];
	unsigned id = 0;

	for(unsigned t=1; t < size; t++)
	{

		std::vector<signed> children;

		switch(buf[t])
		{
		case 1:
			nodeId.push_back(id);
			nodeType.push_back(OP_AND);
			children.push_back(buf[++t]);
			children.push_back(buf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case 2:
			nodeId.push_back(id);
			nodeType.push_back(OP_OR);
			children.push_back(buf[++t]);
			children.push_back(buf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case 3:
			nodeId.push_back(id);
			nodeType.push_back(OP_NOT);
			children.push_back(buf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case 4:
			nodeId.push_back(id);
			nodeType.push_back(OP_LIT);
			children.push_back(buf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case 6:
			nodeExpressions.push_back(m.mkFalse());
			break;
		case 7:
			nodeExpressions.push_back(m.mkTrue());
			break;
		default:
			throw std::exception("RBC command not defined");
		}

		id++;

	}

	std::vector<bool> transformed;
	for(unsigned t=0; t < nodeType.size(); t++)
		transformed.push_back(false);

	while(m.isNil(nodeExpressions[0])){

		for(signed t=(signed)nodeType.size()-1; t >= 0; t--){

			if(transformed[t]) continue;

			switch(nodeType[t]){

			case 1:
				if(!m.isNil(nodeExpressions[childrenId[t][0]]) && !m.isNil(nodeExpressions[childrenId[t][1]])){
					nodeExpressions[t] = m.mkAnd(nodeExpressions[childrenId[t][0]],nodeExpressions[childrenId[t][1]]);
					transformed[t] = true;
				}
				break;
			case 2 :
				if(!m.isNil(nodeExpressions[childrenId[t][0]]) && !m.isNil(nodeExpressions[childrenId[t][1]])){
					nodeExpressions[t] = m.mkOr(nodeExpressions[childrenId[t][0]],nodeExpressions[childrenId[t][1]]);
					transformed[t] = true;
				}
				break;
			case 3:
				if(!m.isNil(nodeExpressions[childrenId[t][0]])){
					nodeExpressions[t] = m.mkNeg(nodeExpressions[childrenId[t][0]]);
					transformed[t] = true;
				}
				break;
			case 4:
				nodeExpressions[t] = m.mkLiteral(childrenId[t][0]);
				break;
			}

		}

	}


	for(unsigned t=1; t < nodeExpressions.size(); t++)
		nodeExpressions[t] = m.mkNil();

	return nodeExpressions[0];

}

//solving methods
bool MPDeSATNode::findDisagreement(void)
{

	bool res = false;

	std::vector<signed> new_trail;
	std::vector<unsigned> global_model(maxVar+1,M_UNDEF);
	std::vector<unsigned> preferencesTrue(maxVar+1,0);
	std::vector<unsigned> preferencesFalse(maxVar+1,0);


	for (unsigned i=0; i< workersId.size(); i++)
	{
		for(unsigned t=0; t < models[i].size(); t++)
		{
			signed mvi = models[i][t];

			if (sharedVariables.occurs(mvi, i)) //it is required to check if it is shared?
			{
				if(mvi < 0){
					preferencesFalse[-mvi]++;
					if(global_model[-mvi] == M_UNDEF || global_model[-mvi] == M_FALSE){
						global_model[-mvi] = M_FALSE;
					} else {
						//conflict in the set of models
						global_model[-mvi] = 4; //conflict

					}

				} else {
					preferencesTrue[mvi]++;
					if(global_model[mvi] == M_UNDEF || global_model[mvi] == M_TRUE){
						global_model[mvi] = M_TRUE;
					} else {
						//conflict in the set of models
						global_model[mvi] = 4; //conflict

					}

				}

			}
		}
	}

	for(signed t=1; t < (signed)global_model.size(); t++){

		if(global_model[t] == 4){

			signed v = t;

			if(preferencesTrue[t] > preferencesFalse[t])
				v = -t;

			if (std::find(trail.begin(), trail.end(), v) == trail.end())
				new_trail.push_back(v);

		}

	}


	for (unsigned i=0; i<new_trail.size(); i++)
	{
		trail.push_back(new_trail[i]);
		reasons.push_back(false);
	}

	global_model.clear();
	preferencesFalse.clear();
	preferencesTrue.clear();

	return new_trail.size()>0;

}