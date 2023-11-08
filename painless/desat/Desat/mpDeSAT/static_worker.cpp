// Copyright (C) 2011 Microsoft Research
// RM Martins, 2012
/*
 * This file contains the new refactored version of mpDeSAT-flat, mpDeSAT-tree and mpDeSAT-tree-flat using a state machine
 * WARNING: It is still under development and should not be used
 */

#include <decomposition_batch.h>
#include "static_worker.h"

StaticWorker::StaticWorker(ExpressionManager &em, const char *fn, int initial_nodes, InterpolationMode interpolation, int masterID, std::vector<int> workersID) :
    SATSolver(em),
    state(IDLE),
    d(NULL),
	overflow_limit(100000),
	numModels(0),
	numConflicts(0)
{

	if (initial_nodes < 2)
		throw std::exception("initial nodes must be larger than 2");

	MPI_Comm_rank(MPI_COMM_WORLD, & rank);
	MPI_Comm_size(MPI_COMM_WORLD, & size);
	filename = fn;
	active_nodes = initial_nodes;
	state = INITIALIZING;
	result = false;

	imode = interpolation;

	master = masterID;
    numWorkers = workersID.size();
	
	for (unsigned i = 0; i < numWorkers; i++)
        workers.push_back(workersID[i]);

	if (isMaster())
		globalSolver = new MiniSAT_1p(em, false, MiniSAT_1p::LIFTING);
	else 
	{
		if(imode == NONE) 
			solver = new Partition(em, sharedVariables, rank, verbosity, false);
		else
		{
			solver = new Partition(em, sharedVariables, rank, verbosity, true);
			solver->setInterpolationMode(imode, sharedVariables);
		}
	}

    for (unsigned i = 0; i<numWorkers; i++)
		occurrences.push_back(new VariableOccurrence());

    interpolants.resize(numWorkers, m.mkNil());
    models.resize(numWorkers);

	if (!isLeaf()) 
	{
		d = new BatchDecomposition();    
        d->setPartitions(numWorkers);
	}
}

StaticWorker::~StaticWorker()
{
    if (d) delete d;
}

bool StaticWorker::addClause(const std::vector<signed> &literals)
{
    assert(!isLeaf());
	unsigned w = d->where(literals);
    unsigned receiving_subordinate = workers[w];
    say(2, "sending clause to %d", receiving_subordinate);

	   sendbuf.clear();
        for (unsigned i = 0; i < literals.size(); i++)
		{
            sendbuf.push_back(literals[i]);
            occurrences[position(workers, receiving_subordinate)]->setOccurs(literals[i]);
		}
        safe_send(receiving_subordinate, CLAUSE_DATA);

	return true;
	
}

void StaticWorker::setVariableMax(unsigned n)
{
    maxVar = n;
    if (d) d->setVariableMax(n);

	if(!isMaster()) solver->setVariableMax(n);
    
    if (isMaster())
    {   
		globalSolver->setVariableMax(n);

        sendbuf.push_back(n);
        for (int i = 1; i < size; i++)
            safe_send(i, INIT_DATA);
    }
}

void StaticWorker::setClauseMax(unsigned n)
{
    maxClauses = n;
    if (d) d->setClauseMax(n);
}

unsigned StaticWorker::numClauses(void) const
{
    throw std::exception("NYI: numClauses");
}

unsigned StaticWorker::numVars(void) const
{
    throw std::exception("NYI: numVars");
}

signed StaticWorker::addVar(void)
{
    throw std::exception("NYI: addVar");
}

void StaticWorker::solveMaster(void)
{    
    say(0,"solving master (%d variables, %d clauses)", globalSolver->numVars(), globalSolver->numClauses());
	saybuffer(0,trail,"Trail: ");

	bool r=globalSolver->solve(trail);  

	while (!r && trail.size()>0)
	{ 
		tmp_clause.clear();
		((MiniSAT_1p*)globalSolver)->getConflict(tmp_clause);

		signed last = 0;
		bool reason = true;
		bool inConflict = false;

		inConflict = std::find(tmp_clause.begin(), tmp_clause.end(), -last)!=tmp_clause.end();    
		while (trail.size()>0 && !inConflict)
		{
			last = trail.back();
			trail.pop_back();
			reason = reasons.back();
			reasons.pop_back();      
			inConflict = std::find(tmp_clause.begin(), tmp_clause.end(), -last)!=tmp_clause.end();
		}

		if (inConflict && !reason)
		{      
			trail.push_back(-last);
			reasons.push_back(true);
		}

		if (trail.size()==0 && reason)
			r = false;
		else
			r=globalSolver->solve(trail);
	}

	if (!r)
	{
		result = false;
		state = DEAD;
	}
	else 
	{
		sendbuf.clear();
		for (int i = 1; i <= (signed)maxVar; i++)
		{
			switch(globalSolver->get(i))
			{
			case M_TRUE: sendbuf.push_back(i); break;
			case M_FALSE: sendbuf.push_back(-i); break;
			default: ;
			}
		}

		for (int i = 0; i < active_nodes; i++)
			if (i != rank)
				safe_send(i, TRAIL_DATA);

		state = IDLE;
	}
}

bool StaticWorker::findDisagreement(void)
{
	std::vector<signed> new_trail;
	std::vector<unsigned> global_model(maxVar+1,M_UNDEF);
	std::vector<unsigned> preferencesTrue(maxVar+1,0);
	std::vector<unsigned> preferencesFalse(maxVar+1,0);
	
	for (unsigned i=0; i < numWorkers; i++)
	{
		for (unsigned t=0; t < models[i].size(); t++)
		{
			signed mvi = models[i][t];

			if (sharedVariables.occurs(mvi, i))
			{
				if (mvi < 0){
				
					preferencesFalse[-mvi]++;
					if (global_model[-mvi] == M_UNDEF || global_model[-mvi] == M_FALSE)
						global_model[-mvi] = M_FALSE;
					else
						global_model[-mvi] = M_CONFLICT;
				
				} else {

					preferencesTrue[mvi]++;
					if (global_model[mvi] == M_UNDEF || global_model[mvi] == M_TRUE)
						global_model[mvi] = M_TRUE;
					else
						global_model[mvi] = M_CONFLICT;

				}

			}

		}
	}

	for (signed t=1; t < (signed)global_model.size(); t++){

		if (global_model[t] == M_CONFLICT){

			signed v = t;
			if (preferencesTrue[t] > preferencesFalse[t]) v = -t;

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

	saybuffer(0,trail,"Trail: ");

	return new_trail.size()>0;
}

void StaticWorker::resolveMaster(void)
{
	assert(isMaster());

	assert(numConflicts + numModels == numWorkers);
	bool had_disagreement = false;

	if (numConflicts == 0)
	{
		had_disagreement = findDisagreement();
		if (had_disagreement)
		{
			sendbuf.clear();
			for (unsigned i = 0; i < trail.size(); i++)
				sendbuf.push_back(trail[i]);

			for (unsigned i = 0; i < numWorkers; i++)
                safe_send(workers[i], TRAIL_DATA);

			state = IDLE;
		}
		else
		{
			result = true;
			state = DEAD;
		}
	}
	else 
	{
		
		for (unsigned i=0; i < numWorkers; i++)
		{
			CExpression itp = interpolants[i];
			if (m.isNil(itp)) continue;

			if (!m.isTrue(itp))
			{

				std::cout << "Node " << rank << " ITP: " << m.toString(itp) << std::endl;
	
				((MiniSAT_1p*)globalSolver)->clearNewClauses();
				globalSolver->addConstraint(itp);      
				((MiniSAT_1p*)globalSolver)->addNewClauses();
			}

			interpolants[i] = m.mkNil();
		}

		state = SOLVING_MASTER;
	}

	numModels = 0;
	numConflicts = 0;

	models.resize(numWorkers);
	interpolants.resize(numWorkers,m.mkNil());
}

void StaticWorker::resolveWorker(void)
{    
	assert(!isMaster());
	say(0,"solving worker (%d variables, %d clauses)", solver->numVars(), solver->numClauses());
	saybuffer(0,trail,"Trail: ");

	workerInterpolant = m.mkNil();
	Expression t = solver->getInterpolant(trail);

	bool r=solver->em().isTrue(t);

	while (!r && trail.size()>0)
	{ 
		tmp_clause.clear();
		((MiniSAT_1p*)solver->getSolver())->getConflict(tmp_clause);

		signed last = 0;
		bool reason = true;
		bool inConflict= false;

		inConflict = std::find(tmp_clause.begin(), tmp_clause.end(), -last)!=tmp_clause.end();    
		while (trail.size()>0 && !inConflict)
		{
			last = trail.back();
			trail.pop_back();
			reason = reasons.back();
			reasons.pop_back();      
			inConflict = std::find(tmp_clause.begin(), tmp_clause.end(), -last)!=tmp_clause.end();
		}

		if (inConflict && !reason)
		{      
			trail.push_back(-last);
			reasons.push_back(true);
		}

		if (trail.size()==0 && reason)
			r = false;
		else
			r=solver->solve(trail);

	}
	
	//t is already false when the unsatisfiability does not depend on the assumptions!
	
	if (!solver->em().isTrue(t))
	{
		workerInterpolant = m.mkNil();
		workerInterpolant = m.duplicate(t, solver->em());
		t = m.mkNil();

		bool overflow_problem = m.sizeRBC(workerInterpolant, overflow_limit);

		if (overflow_problem)
        {
			workerInterpolant = m.mkNil();

			vec<Lit> &cfs = ((MiniSAT_1p*)solver)->conflict;
			std::vector<signed> cf_tmp;
			for (int i=0; i<cfs.size(); i++)        
				cf_tmp.push_back(litToLiteral(cfs[i]));


			if (cf_tmp.size()==0)
				workerInterpolant = m.mkFalse();
			else if (cf_tmp.size()==1)
				workerInterpolant = m.mkLiteral(cf_tmp[0]);
			else
			{
				workerInterpolant = m.mkOr(m.mkLiteral(cf_tmp[0]), m.mkLiteral(cf_tmp[1]));
				for (unsigned i=2; i<cf_tmp.size(); i++)
				{
					workerInterpolant = m.mkOr(workerInterpolant, m.mkLiteral(cf_tmp[i]));
				}

			}
		}

		std::cout << "Node " << rank << " ITP: " << m.toString(workerInterpolant) << std::endl;
		m.fromITPtoBuffer(sendbuf, workerInterpolant);

		interpolants_exported++;
		safe_send(master, CONFLICT_DATA);

	}
	else 
	{

		sendbuf.clear();
		if (isLeaf())
		{

		for (int i = 1; i <= (signed)maxVar; i++)
		{
			switch(solver->get(i))
			{
			case M_TRUE: sendbuf.push_back(i); if (verbosity > 1) output << i << " "; break;
			case M_FALSE: sendbuf.push_back(-i); if (verbosity > 1) output << -i << " "; break;
			default: ;
			}
		}

		saybuffer(0, sendbuf, "Model: ");
		
		safe_send(master, MODEL_DATA);
		}
		else 
		{

			for (unsigned i = 0; i < trail.size(); i++)
				sendbuf.push_back(trail[i]);

			for (unsigned i = 0; i < numWorkers; i++)
                safe_send(workers[i], MODEL_DATA);
		}
	}

	sendbuf.clear();
	workerInterpolant = m.mkNil();
	t = m.mkNil();

	state = IDLE;
}

Expression StaticWorker::fromBufferToITP(void)
{

	typedef enum { OP_UNDEF=0, OP_AND, OP_OR, OP_NOT, OP_LIT, OP_NIL, OP_FALSE, OP_TRUE } operations;

	std::vector<signed> nodeId;
	std::vector< std::vector< signed> > childrenId;
	std::vector<signed> nodeType;
	std::vector<Expression> nodeExpressions;

	unsigned id = 0;

	for(unsigned t=0; t < recvbuf.size(); t++)
	{

		std::vector<signed> children;

		switch(recvbuf[t])
		{
		case OP_AND:
			nodeId.push_back(id);
			nodeType.push_back(OP_AND);
			children.push_back(recvbuf[++t]);
			children.push_back(recvbuf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case OP_OR:
			nodeId.push_back(id);
			nodeType.push_back(OP_OR);
			children.push_back(recvbuf[++t]);
			children.push_back(recvbuf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case OP_NOT:
			nodeId.push_back(id);
			nodeType.push_back(OP_NOT);
			children.push_back(recvbuf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case OP_LIT:
			nodeId.push_back(id);
			nodeType.push_back(OP_LIT);
			children.push_back(recvbuf[++t]);
			childrenId.push_back(children);
			nodeExpressions.push_back(m.mkNil());
			break;
		case OP_FALSE:
			nodeExpressions.push_back(m.mkFalse());
			break;
		case OP_TRUE:
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

			case OP_AND:
				if(!m.isNil(nodeExpressions[childrenId[t][0]]) && !m.isNil(nodeExpressions[childrenId[t][1]])){
					nodeExpressions[t] = m.mkAnd(nodeExpressions[childrenId[t][0]],nodeExpressions[childrenId[t][1]]);
					transformed[t] = true;
				}
				break;
			case OP_OR :
				if(!m.isNil(nodeExpressions[childrenId[t][0]]) && !m.isNil(nodeExpressions[childrenId[t][1]])){
					nodeExpressions[t] = m.mkOr(nodeExpressions[childrenId[t][0]],nodeExpressions[childrenId[t][1]]);
					transformed[t] = true;
				}
				break;
			case OP_NOT:
				if(!m.isNil(nodeExpressions[childrenId[t][0]])){
					nodeExpressions[t] = m.mkNeg(nodeExpressions[childrenId[t][0]]);
					transformed[t] = true;
				}
				break;
			case OP_LIT:
				nodeExpressions[t] = m.mkLiteral(childrenId[t][0]);
				break;
			}

		}

	}


	for(unsigned t=1; t < nodeExpressions.size(); t++)
		nodeExpressions[t] = m.mkNil();

	return nodeExpressions[0];

}

void StaticWorker::buildSharedVariables(void)
{
		for (unsigned i = 0; i < numWorkers; i++)
			sharedVariables.add(occurrences[i]);

		sharedVariables.update();

		sendbuf.clear();
		for (unsigned i = 1; i <= maxVar; i++)
		{
			if (sharedVariables.isShared(i))
				sendbuf.push_back(i);
		}

		for (unsigned i = 0; i < numWorkers; i++)
            safe_send(workers[i], SHARED_VARIABLES_DATA);

		state = SOLVING_MASTER;
}

void StaticWorker::importSharedVariables(void)
{
	
	for (unsigned i = 0; i < recvbuf.size(); i++)
		sharedVariables.setShared(recvbuf[i],true);
}

bool StaticWorker::solve(void)
{    
	say(0,"solving");
    bool done = false;

    while (!done)
    {
        switch (state)
        {
        case INITIALIZING:
            if (isMaster())
            {
                say(0,"reading file");
                if (!readDimacsFile(filename))
                    state = DEAD;
                else
                    state = SHARED_VARIABLES;
            }
            else
                getMessage();
            break;
        case IDLE:            
            getMessage();
            break;
		case SHARED_VARIABLES:
			buildSharedVariables();
			break;
        case SOLVING_MASTER:            
            solveMaster();
            break;
        case SOLVING_WORKER:            
            resolveWorker();
            break;
		case DEAD_LEAF:
			assert(isMaster());
			result = false;
			state = DEAD;
			break;
        case DEAD:
            if (isMaster()) 
            {
                sendbuf.clear();
                sendbuf.push_back(result ? 1 : 0);
                for (int i = 0; i < size; i++)
                    if (i != rank) safe_send(i, TERMINATE_DATA);
            }
            
            done = true;
            break;
        default:
            throw std::exception("unknown state");
        }
    }
    
    return result;
}

void StaticWorker::getMessage()
{
    Tag tag = receive();
    switch (tag)
    {    
    case INIT_DATA: 
        say(0,"received init message");
        setVariableMax(recvbuf[0]);
        state = IDLE;
        break;
    case CLAUSE_DATA: 
        say(0,"received clause message");
		if (isLeaf())
		{
			if (!solver->addClause(recvbuf))
			{
				sendbuf.clear();
				safe_send(0, DEAD_LEAF);
				state = IDLE;
			}
		} 
		else 
			addClause(recvbuf);
		break;
    case TRAIL_DATA:
		say(0,"received trail data message");
        //assumptions.clear();
		trail.clear();
		reasons.clear();
        for (unsigned i = 0; i < recvbuf.size(); i++)
		{
            trail.push_back(recvbuf[i]);
			reasons.push_back(true);
		}
		state = SOLVING_WORKER;
        break;
    case MODEL_DATA:
		say(0,"received model data message");
        assert(!isLeaf());
        models[position(workers, status.MPI_SOURCE)] = recvbuf;
		numModels++; 
        if (numModels + numConflicts >= numWorkers)
            resolveMaster();
        break;
    case CONFLICT_DATA:
		say(0,"received conflict data message");
        assert(!isLeaf());
        interpolants[position(workers, status.MPI_SOURCE)] = fromBufferToITP();
		numConflicts++;
        
		if (numModels + numConflicts >= numWorkers)
            resolveMaster();
        break;
	case SHARED_VARIABLES_DATA:
		say(0,"shared variables data message");
		for (unsigned i = 0; i < recvbuf.size(); i++)
			sharedVariables.setShared(recvbuf[i],true);
		state = IDLE;
		break;
    case TERMINATE_DATA:
		say(0,"terminate data message");
        result = (recvbuf[0] != 0);
        state = DEAD;
        break;
    default:
        std::exception("unknown mpi tag");
    }
}

void StaticWorker::safe_send(int to, int tag)
{  
    if (MPI_Send((void*) sendbuf.data(), sendbuf.size(), MPI_INT, to, tag, MPI_COMM_WORLD)!=MPI_SUCCESS)
		throw std::exception("MPI send failure");    
}

StaticWorker::Tag StaticWorker::receive()
{
    say(2, "waiting for message");

    Tag res = (Tag) -1;
    	
    while (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status)!=MPI_SUCCESS)
		say(0, "Probe failure while waiting for command.");
    
    int c = 0;
    MPI_Get_count(&status, MPI_INT, &c); 
    recvbuf.resize(c);
    
    res = (Tag) status.MPI_TAG;

    while (MPI_Recv(recvbuf.data(), status.count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status)!=MPI_SUCCESS)
		say(0, "Receive failure while waiting for command.");

    return res;
}