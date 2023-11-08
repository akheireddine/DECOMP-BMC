// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger & RM Martins, 2012
/*
 * This file contains the simple worker:
 * (i) it does not use interpolation and it was the start of the new framework
 * WARNING: Deprecated file! Should be used only if we want to extend it from scratch into a different direction.
 */

#include <mpi.h>
#define MPI_DEFAULT_MASTER 0;

#include <decomposition_batch.h>

#include "simple_worker.h"

SimpleWorker::SimpleWorker(const char *fn, int initial_nodes) :
    SATSolver(em),
    solver(em, false, MiniSAT_1p::LIFTING),
    state(IDLE),
    d(NULL)
{
    if (initial_nodes < 1)
        throw std::exception("initial nodes must be larger than 1");

    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    MPI_Comm_size(MPI_COMM_WORLD, & size);
    master = MPI_DEFAULT_MASTER;
    filename = fn;
    active_nodes = initial_nodes;
    state = INITIALIZING;
    result = false;

    if (isMaster()) 
    {
        d = new BatchDecomposition();    
        d->setPartitions(active_nodes);
    }
}

SimpleWorker::~SimpleWorker()
{
    if (d) delete d;
}

bool SimpleWorker::addClause(const std::vector<signed> &literals)
{
    unsigned w = d->where(literals);
    
    if (w == 0)
    {
        return solver.addClause(literals);
    }
    else
    {
        say("sending clause to %d", w);
        sendbuf.clear();
        for (unsigned i = 0; i < literals.size(); i++)
            sendbuf.push_back(literals[i]);
        safe_send(w, CLAUSE_DATA);
        return true;
    }
}

void SimpleWorker::setVariableMax(unsigned n)
{
    maxVar = n;
    if (d) d->setVariableMax(n);
    solver.setVariableMax(n);
    level.resize(n+1, 0);

    if (isMaster())
    {        
        sendbuf.push_back(n);
        for (int i = 1; i < size; i++)
            safe_send(i, INIT_DATA);
    }
}

void SimpleWorker::setClauseMax(unsigned n)
{
    maxClauses = n;
    if (d) d->setClauseMax(n);
}

unsigned SimpleWorker::numClauses(void) const
{
    throw std::exception("NYI: numClauses");
}

unsigned SimpleWorker::numVars(void) const
{
    throw std::exception("NYI: numVars");
}

signed SimpleWorker::addVar(void)
{
    throw std::exception("NYI: addVar");
}

void SimpleWorker::solveMaster(void)
{    
    say("solving master (%d variables, %d clauses)", solver.numVars(), solver.numClauses());

    if (solver.solve(trail))
    {
        sendbuf.clear();
        for (int i = 1; i <= (signed)maxVar; i++)
        {
            switch(solver.get(i))
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
    else
    {
        result = false;
        state = DEAD;
    }
}

void SimpleWorker::solveWorker(void)
{
    say("solving worker (%d variables, %d clauses)", solver.numVars(), solver.numClauses());

    sendbuf.clear();
    
    if (solver.solve(assumptions))
    {
        if (verbosity > 0)
        {
            output.str("");
            output << "Model: ";
        }
        
        say("sending solution");
        for (int i = 1; i <= (signed)maxVar; i++)
        {
            switch(solver.get(i))
            {
            case M_TRUE: sendbuf.push_back(i); if (verbosity > 0) output << i << " "; break;
            case M_FALSE: sendbuf.push_back(-i); if (verbosity > 0) output << -i << " "; break;
            default: ;
            }
        }
        safe_send(master, MODEL_DATA);

        if (verbosity > 0) 
            say(output.str().c_str());
    }
    else
    {
        say("sending conflict");
        solver.getConflict(sendbuf);
        safe_send(master, CONFLICT_DATA);
    }

    state = IDLE;
}

void SimpleWorker::resolve(void)
{
    for (unsigned i = 0; i < conflicts.size(); i++) 
    {
        if (!solver.addClause(conflicts[i]))
        {
            result = false;
            state = DEAD;
        }
    }

    if (state != DEAD) 
        state = SOLVING_MASTER;    
}

void SimpleWorker::conflictAnalysis(void)
{
    if (conflicts.size() == 0)
    {
        // only models.
        bool had_disagreement = false;
        assumptions.clear();
        tmp_model.resize(maxVar+1, M_UNDEF);
        
        for (unsigned i = 0; i < tmp_model.size(); i++)
            tmp_model[i] = M_UNDEF;

        for (unsigned i = 0; i < models.size(); i++)
        {
            std::vector<signed> & model_i = models[i];
            for (unsigned j = 0; j < model_i.size(); j++)
            {
                signed l = model_i[j];
                if (l == 0) continue;
                signed v = (l<0) ? -l : l;
                ModelValue nm = (l < 0) ? M_FALSE : M_TRUE;

                if (tmp_model[v] == M_UNDEF) 
                {
                    tmp_model[v] = nm;
                    assumptions.push_back(l);
                }
                else if (tmp_model[v] != nm)
                {
                    level[v] = trail.size();
                    trail.push_back(v);
                    reasons.push_back(std::vector<signed>());
                    had_disagreement = true;
                }
            }
        }
        
        say("disagreements: %d", had_disagreement ? 1 : 0);

        if (!had_disagreement) 
        {        
            bool r = solver.solve(assumptions);
            say("model extended: %d", r ? 1 : 0);

            if (r)
            {                
                result = true;
                state = DEAD;
            }
            else
            {
                // conflict.
                conflicts.push_back(std::vector<signed>());
                solver.getConflict(conflicts.back());
                resolve();
            }
        }
    }
    else
    {
        resolve();
    }

    models.clear();
    conflicts.clear();
}

bool SimpleWorker::solve(void)
{    
    bool done = false;

    while (!done)
    {
        switch (state)
        {
        case INITIALIZING:
            if (isMaster())
            {
                say("reading file");
                if (!readDimacsFile(filename))
                    state = DEAD;
                else
                    state = SOLVING_MASTER;
            }
            else
                getMessage();
            break;
        case IDLE:            
            getMessage();
            break;
        case SOLVING_MASTER:            
            solveMaster();
            break;
        case SOLVING_WORKER:            
            solveWorker();
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

void SimpleWorker::getMessage()
{
    Tag tag = receive();
    switch (tag)
    {    
    case INIT_DATA: 
        if (verbosity > 1) say("received init message");
        setVariableMax(recvbuf[0]);
        state = IDLE;
        break;
    case CLAUSE_DATA: 
        if (verbosity > 1) say("received clause message");        
        if (!solver.addClause(recvbuf))
            throw std::exception("NYI");
        break;
    case TRAIL_DATA:
        assumptions.clear();
        for (unsigned i = 0; i < recvbuf.size(); i++)
            assumptions.push_back(recvbuf[i]);
        state = SOLVING_WORKER;
        break;
    case MODEL_DATA:
        assert(isMaster());
        models.push_back(recvbuf);
        if ((models.size() + conflicts.size()) >= (unsigned)active_nodes-1)
            conflictAnalysis();
        break;
    case CONFLICT_DATA:
        assert(isMaster());        
        conflicts.push_back(recvbuf);
        if ((models.size() + conflicts.size()) >= (unsigned)active_nodes-1)
            conflictAnalysis();
        break;    
    case TERMINATE_DATA:
        result = (recvbuf[0] != 0);
        state = DEAD;
        break;
    default:
        std::exception("unknown mpi tag");
    }
}

void SimpleWorker::safe_send(int to, int tag)
{  
    if (MPI_Send((void*) sendbuf.data(), sendbuf.size(), MPI_INT, to, tag, MPI_COMM_WORLD)!=MPI_SUCCESS)
		throw std::exception("MPI send failure");    
}

SimpleWorker::Tag SimpleWorker::receive()
{
    say("waiting for message");

    Tag res = (Tag) -1;
    	
    while (MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status)!=MPI_SUCCESS)
	{
		std::cout << "Node " << rank << ": Probe failure while waiting for command." << std::endl;
	}
    
    int c = 0;
    MPI_Get_count(&status, MPI_INT, &c); 
    recvbuf.resize(c);
    
    res = (Tag) status.MPI_TAG;

    while (MPI_Recv(recvbuf.data(), status.count, MPI_INT, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status)!=MPI_SUCCESS)
	{		
		std::cout << "Node " << rank << ": Receive failure while waiting for command." << std::endl;
	}

    return res;
}