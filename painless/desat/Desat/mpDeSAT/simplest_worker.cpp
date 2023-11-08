// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger & RM Martins, 2012
/*
 * This file contains the latest version of our solver (mpDeSAT-memory):
 * (i)   it does not use interpolation;
 * (ii)  it uses a flat structure;
 * (iii) it can start the solving process with "w" nodes and will split the formula automatically if a cutoff is reached
 * (iv)  further testing to the memory limit of MiniSAT should be done (including the way we split the database)
 * (v) ... still a lot of work needs to be done so it is competitive with MiniSAT1.14p
 * WARNING: Do not use HALFLEARNT clause since it makes the SAT solver performs much worse
 */

#include <mpi.h>
#include <iostream>
#include <iomanip>

#include <decomposition_batch.h>
#include "returnvalue.h"
#include "simplest_worker.h"

#define MPI_DEFAULT_MASTER 0;
//#define LEARN_FROM_ASSUMPTIONS 1; //Should be the target of further testing - not sure if it is working properly
#define FORCE_MASTER_CHANGE 1;
//#define CLEAN_LEARNED_CLAUSES 1;

SimplestWorker::SimplestWorker(const char *fn, int initial_nodes, int working_nodes, int mlimit) :
    SATSolver(em),
    solver(em, false, MiniSAT_1p::LIFTING),
    state(IDLE),
    d(NULL),
	master_changed(0),
	rounds(0),
    numConflicts(0),
    numModels(0),
	memory_limit(mlimit),
	memory_exhausted(false),
	startTime(0)
{
    if (initial_nodes < 1)
        throw std::exception("initial nodes must be larger than 1");

    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
    MPI_Comm_size(MPI_COMM_WORLD, & size);
    master = MPI_DEFAULT_MASTER;
    filename = fn;
    active_nodes = working_nodes;
	total_nodes = initial_nodes;
    state = INITIALIZING;
    result = true;
	active_status = false;

    conflicts.resize(active_nodes);
    models.resize(active_nodes);

    if (isMaster()) 
    {
		active_status = true;
        d = new BatchDecomposition();    
        d->setPartitions(active_nodes);
    }

	solver.setMemoryLimit(memory_limit);
}

SimplestWorker::~SimplestWorker()
{
    if (d) delete d;
}

bool SimplestWorker::safe_addClause(const std::vector<signed> &literals, int type)
{
	try
	{
			tmp_clause.clear();
			for (unsigned i = 0; i < literals.size(); i++)
			{
				unsigned var = literals[i] < 0 ? -literals[i] : literals[i];
					while (var > solver.numVars())
						solver.addVar();
				tmp_clause.push_back(literals[i]);

#ifdef LEARN_FROM_ASSUMPTIONS 
				if(!variableOccurs[var]) variableOccurs[var] = true;
#endif
			}

		switch (type)
		{
			case ORIGINAL_CLAUSE:
				return solver.addClause(tmp_clause);
			case LEARNT_CLAUSE:
				return solver.addLearntClause(tmp_clause);
			case HALFLEARNT_CLAUSE:
				return solver.addHalfLearntClause(tmp_clause);
			default:
				throw std::exception("unknown clause type");
		}

	}
	catch (std::bad_alloc &)
	{
		
		splitMaster();

		switch (type)
		{
		case ORIGINAL_CLAUSE:
			return solver.addClause(tmp_clause);
		case LEARNT_CLAUSE:
			return solver.addLearntClause(tmp_clause);
		case HALFLEARNT_CLAUSE:
			return solver.addHalfLearntClause(tmp_clause);
		default:
			throw std::exception("unknown clause type");
		}
	}

}

void SimplestWorker::splitWorker(void)
{
	say(2,"splitWorker init: variables = %d\t clauses = %d\t learnts = %d \t halflearnts = %d",solver.nVars(), solver.nClauses(), solver.nLearnts(), solver.nHalfLearnts());

	unsigned num_clauses = solver.nClauses();
	unsigned num_learnts = solver.nLearnts();
	unsigned num_halflearnts = solver.nHalfLearnts();

	master = recvbuf[0];
	active_nodes = recvbuf[1];

	unsigned next_worker = active_nodes-1;

	assumptions.clear();
	for (unsigned i = 2; i < recvbuf.size(); i++)
		assumptions.push_back(recvbuf[i]);

	for (unsigned i = 0; i < ceil((double)num_clauses/2); i++)
	{
		say(2, "sending clause to %d", next_worker);
		sendbuf.clear();
		//solver.getAndRemoveLast(sendbuf, ORIGINAL_CLAUSE);
		int position = rand() % num_clauses;
		solver.getAndRemovePosition(sendbuf, position, ORIGINAL_CLAUSE);
		safe_send(next_worker, CLAUSE_DATA);
	}

#ifdef CLEAN_LEARNED_CLAUSES
	solver.cleanLearnts();
#else
	for (unsigned i = 0; i < ceil((double)num_learnts/2); i++)
	{
		say(2, "sending learnt clause to %d", next_worker);
		sendbuf.clear();
		//solver.getAndRemoveLast(sendbuf, LEARNT_CLAUSE);
		int position = rand() % num_learnts;
		solver.getAndRemovePosition(sendbuf, position, LEARNT_CLAUSE);
		if (sendbuf.empty()){
			say(2,"no learnt clauses are available to be sent to %d",next_worker);
			break;
		}
		safe_send(next_worker, LEARNT_DATA);
	}
#endif

	for (unsigned i = 0; i < ceil((double)num_halflearnts/2); i++)
	{
		say(2, "sending half learnt clause to %d", next_worker);
		sendbuf.clear();
		//solver.getAndRemoveLast(sendbuf, HALFLEARNT_CLAUSE);
		int position = rand() % num_halflearnts;
		solver.getAndRemovePosition(sendbuf, position, HALFLEARNT_CLAUSE);
		safe_send(next_worker, HALFLEARNT_DATA);	
	}

	say(2,"splitMaster end: variables = %d\t clauses = %d\t learnts = %d \t halflearnts = %d",solver.nVars(), solver.nClauses(), solver.nLearnts(), solver.nHalfLearnts());
	
	conflicts.resize(active_nodes);
	models.resize(active_nodes);

#ifdef LEARN_FROM_ASSUMPTIONS 
	variableOccurs.resize(maxVar+1,false);
	solver.updateFormulaVariables(variableOccurs);
#endif

	sendbuf.clear();
	sendbuf.push_back(master);
	sendbuf.push_back(active_nodes);
	for (unsigned i=0; i < assumptions.size(); i++)
		sendbuf.push_back(assumptions[i]);

	safe_send(next_worker, TRAIL_DATA);

	solveWorker();
}

void SimplestWorker::splitMaster(void)
{
	say(2,"splitMaster init: variables = %d\t clauses = %d\t learnts = %d \t halflearnts = %d",solver.nVars(), solver.nClauses(), solver.nLearnts(), solver.nHalfLearnts());

	if (active_nodes == total_nodes)
	{
		say(0, "solver memory = %d\t active nodes = %d",solver.getUsedMemory(), active_nodes);
		say(0,"MEMORY EXHAUSTED");
		MPI_Abort(MPI_COMM_WORLD, S_ERROR);
		
	}
	else
	{
		unsigned num_clauses = solver.nClauses();
		unsigned num_learnts = solver.nLearnts();
		unsigned num_halflearnts = solver.nHalfLearnts();

		unsigned next_worker = active_nodes;
		active_nodes++;

		for (unsigned i = 0; i < ceil((double)num_clauses/2); i++)
		{
			say(2, "sending clause to %d", next_worker);
			sendbuf.clear();
			//solver.getAndRemoveLast(sendbuf, ORIGINAL_CLAUSE);
			int position = rand() % num_clauses;
			solver.getAndRemovePosition(sendbuf, position, ORIGINAL_CLAUSE);
			safe_send(next_worker, CLAUSE_DATA);
		}

#ifdef CLEAN_LEARNED_CLAUSES
		solver.cleanLearnts();
#else
		for (unsigned i = 0; i < ceil((double)num_learnts/2); i++)
		{
			say(2, "sending learnt clause to %d", next_worker);
			sendbuf.clear();
			//solver.getAndRemoveLast(sendbuf, LEARNT_CLAUSE);
			int position = rand() % num_learnts;
			solver.getAndRemovePosition(sendbuf, position, LEARNT_CLAUSE);
			if(sendbuf.empty()) break;
			safe_send(next_worker, LEARNT_DATA);
		}
#endif

		for (unsigned i = 0; i < ceil((double)num_halflearnts/2); i++)
		{
			say(2, "sending half learnt clause to %d", next_worker);
			sendbuf.clear();
			//solver.getAndRemoveLast(sendbuf, HALFLEARNT_CLAUSE);
			int position = rand() % num_halflearnts;
			solver.getAndRemovePosition(sendbuf, position, HALFLEARNT_CLAUSE);
			safe_send(next_worker, HALFLEARNT_DATA);	
		}

#ifdef LEARN_FROM_ASSUMPTIONS 
		variableOccurs.resize(maxVar+1,false);
		solver.updateFormulaVariables(variableOccurs);
#endif

		conflicts.resize(active_nodes);
		models.resize(active_nodes);

		say(2,"splitMaster end: variables = %d\t clauses = %d\t learnts = %d \t halflearnts = %d",solver.nVars(), solver.nClauses(), solver.nLearnts(), solver.nHalfLearnts());
		//solver.printDB(rank);
		
	}

}

bool SimplestWorker::safe_solve(const std::vector<signed> &trail)
{
	try
	{
		//solver.setMemoryLimit(memory_limit);

		if (trail.empty()) return solver.solve();
		else return solver.solve(trail);
	}
	catch (std::bad_alloc &)
	{
		if (isMaster())
		{
			splitMaster();
			return safe_solve(tmp_trail);
		}
		else 
		{
			memory_exhausted = true;
			return false;
		}
	}
}

bool SimplestWorker::addClause(const std::vector<signed> &literals)
{
    unsigned w = d->where(literals);
    
	if (w == 0)
	{
		return safe_addClause(literals, ORIGINAL_CLAUSE);
	}
    else
    {
        say(2, "sending clause from file to %d", w);
        sendbuf.clear();
        for (unsigned i = 0; i < literals.size(); i++)
            sendbuf.push_back(literals[i]);
        safe_send(w, CLAUSE_DATA);
        return true;
    }
}

void SimplestWorker::setVariableMax(unsigned n)
{
    maxVar = n;
    if (d) d->setVariableMax(n);
	solver.setVariableMax(n);

#ifdef LEARN_FROM_ASSUMPTIONS 
	variableOccurs.resize(maxVar+1, false);
#endif

    if (isMaster())
    {        
        sendbuf.push_back(n);
        for (int i = 1; i < size; i++)
            safe_send(i, INIT_DATA);
    }
}

void SimplestWorker::setClauseMax(unsigned n)
{
    maxClauses = n;
    if (d) d->setClauseMax(n);
}

unsigned SimplestWorker::numClauses(void) const
{
    throw std::exception("NYI: numClauses");
}

unsigned SimplestWorker::numVars(void) const
{
    throw std::exception("NYI: numVars");
}

signed SimplestWorker::addVar(void)
{
    throw std::exception("NYI: addVar");
}

void SimplestWorker::getStats(void)
{
	std::cout << "Active nodes: " << active_nodes << std::endl;
	std::cout << "Total nodes: " << total_nodes << std::endl;
	std::cout << "Rounds: " << rounds << std::endl;
	std::cout << "Master changed: " << master_changed << std::endl;
	std::cout << "Time: " << std::setprecision(3) << (clock() - startTime)/(double)CLOCKS_PER_SEC << " sec" << std::endl;
}

void SimplestWorker::solveMaster(void)
{    
    say(2, "solving master (%d variables, %d clauses)", solver.numVars(), solver.numClauses());    

    numConflicts = numModels = 0;
	
	if (safe_solve(tmp_trail))
	{

		if(active_nodes == 1) 
		{
			//only the master is active

			result = true;
			state = DEAD;
		}
		else 
		{
			sendbuf.clear();
			sendbuf.push_back(rank);
			sendbuf.push_back(active_nodes);
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
	}
	else
	{
		result = false;
		state = DEAD;
	}

}

void SimplestWorker::solveWorker(void)
{
    say(1, "solving worker (%d variables, %d clauses)", solver.numVars(), solver.numClauses());

    sendbuf.clear();

    if (!solver.okay())
	{
		// empty conflict
		say(2, "sending empty conflict to %d", master);
		safe_send(master, CONFLICT_DATA);
	}
	else {

		if (safe_solve(assumptions))
		{
			if (verbosity > 1)
			{
				output.str("");
				output << "Model: ";
			}

			say(2, "sending solution");
			for (int i = 1; i <= (signed)maxVar; i++)
			{
				switch(solver.get(i))
				{
				case M_TRUE: sendbuf.push_back(i); if (verbosity > 1) output << i << " "; break;
				case M_FALSE: sendbuf.push_back(-i); if (verbosity > 1) output << -i << " "; break;
				default: ;
				}
			}
			safe_send(master, MODEL_DATA);

			say(2, output.str().c_str());
		}
		else
		{
			if (memory_exhausted)
			{
				say(2, "sending out of memory to %d", master);
				sendbuf.clear();
				safe_send(master, MEMORY_DATA);
				memory_exhausted = false;
			}
			else
			{
				say(2, "sending conflict to %d", master);

				solver.getConflict(sendbuf);

				saybuffer(1, sendbuf, "conflict: ");
				safe_send(master, CONFLICT_DATA);

#ifdef LEARN_FROM_ASSUMPTIONS
				worker_clause.clear();
				for (unsigned i = 0; i < assumptions.size(); i++)
				{
					unsigned var = assumptions[i] < 0 ? -assumptions[i] : assumptions[i];
					if(variableOccurs[var])
						worker_clause.push_back(-assumptions[i]);
				}
				//safe_addClause(worker_clause, HALFLEARNT_CLAUSE);
				safe_addClause(worker_clause, LEARNT_CLAUSE);

				saybuffer(2, assumptions, "assumptions:");
				saybuffer(2, worker_clause, "learning conflict from the assumptions:");
#endif
			}
		}

	}

    state = IDLE;
}

void SimplestWorker::resolve(void)
{
    if (result == false)
    {
        say(2, "Resolve: already false.");
        state = DEAD;
    }
    else
    {
        for (unsigned i = 0; i < conflicts.size(); i++) 
        {
            if (conflicts[i].size() != 0 &&
				!safe_addClause(conflicts[i], ORIGINAL_CLAUSE))
				//!safe_addClause(conflicts[i], HALFLEARNT_CLAUSE))
            {
                result = false;
                state = DEAD;
            }
        }
    }

    if (state != DEAD) 
        state = SOLVING_MASTER;    
}

void SimplestWorker::conflictAnalysis(void)
{
    assert(numConflicts + numModels == (unsigned)(active_nodes-1));

    bool had_disagreement = false;
    unsigned new_master = rank;

    say(1, "Conflict analysis: %d conflicts, %d models", numConflicts, numModels);

    if (numConflicts == 0)
    {
        // every worker sent a model.
        assumptions.clear();
        tmp_model.resize(maxVar+1, M_UNDEF);		        
        tmp_from.resize(maxVar+1, -1);

        for (unsigned i = 0; i < tmp_model.size(); i++)
            tmp_model[i] = M_UNDEF;

        for (unsigned i = 0; i < models.size(); i++)
        {
            const std::vector<signed> & model_i = models[i];
            for (unsigned j = 0; j < model_i.size(); j++)
            {
                signed l = model_i[j];
                if (l == 0) continue;
                signed v = (l<0) ? -l : l;
                ModelValue nm = (l < 0) ? M_FALSE : M_TRUE;

                if (tmp_model[v] == M_UNDEF) 
                {
                    tmp_model[v] = nm;
                    tmp_from[v] = i;
                    assumptions.push_back(l);
                }
                else if (tmp_model[v] != nm)
                {
                    had_disagreement = true;                    
					new_master = (rand() >= (RAND_MAX/2)) ? tmp_from[v] : i;
					break;
                }
            }
        }
    }

    say(2, "disagreements: %s (new master: %d)", had_disagreement ? "yes" : "no", new_master);

    if (!had_disagreement && numConflicts == 0) 
    {
        // No disagreements + no conflicts.

		if (safe_solve(assumptions))
		{                
			result = true;
			state = DEAD;
		}
		else
		{
			// workers ok, but local conflict.            
			solver.getConflict(conflicts[rank]);
			numConflicts++;
		}		
	}
    
    if (numConflicts > 0) 
        resolve(); // resolve all conflicts.
    
    if (state != DEAD)
    {
        // Should we change the master?    
        // If we had a disagreement, we're forced to change;
        // otherwise we can chose to give up our master status.
	    if (had_disagreement) 
	    {			
		    assumptions.clear();

		    sendbuf.clear();
		    sendbuf.push_back(new_master);
			sendbuf.push_back(active_nodes);
		    for (unsigned i = 0; i < models[new_master].size() ; i++)
		    {
			    assumptions.push_back(models[new_master][i]);
			    sendbuf.push_back(models[new_master][i]);
		    }		    
	    }
#ifdef FORCE_MASTER_CHANGE // enable unforced master changes?
        else
        {
            // do we want to give up the master status? Only if master_changed < log(rounds).
            unsigned status = master_changed;
            unsigned target = (unsigned) (log((double)rounds));
            // unsigned target = rounds / (2 * active_nodes);
            say(2, "rounds = %d, target = %d, status = %d", rounds, target, status);
            if (status < target)
            {
                new_master = (rank + 1) % active_nodes;
                say(1, "voluntarily giving master to %d after %d rounds (status=%d < target=%d)", new_master, rounds, status, target);
            
		        sendbuf.clear();
		        sendbuf.push_back(new_master);
				sendbuf.push_back(active_nodes);
		        for (unsigned i = 0; i < assumptions.size() ; i++)
                    sendbuf.push_back(assumptions[i]);
            }
        }
#endif

        if (new_master != rank)
        {
            if (verbosity > 1) {
                output.str("");
                output << "Transferring current assignment: ";
                for (unsigned i = 0; i < sendbuf.size(); i++)
                    output << " " << sendbuf[i];
                say(2, output.str().c_str());
            }

            for (int i = 0; i < active_nodes ; i++)
			    if (i != rank) safe_send(i,TRAIL_DATA);

            master_changed++;
		    master = new_master;
		    state = SOLVING_WORKER;
        }
    }

    // Clean up the conflict/models structures
    numConflicts = numModels = 0;
    for (unsigned i = 0; i < conflicts.size(); i++)
        conflicts[i].clear();
    for (unsigned i = 0; i < models.size(); i++)
        models[i].clear();
}

bool SimplestWorker::solve(void)
{    
    bool done = false;
	clock_t startTime = clock();

    while (!done)
    {        
        switch (state)
        {
        case INITIALIZING:
            if (isMaster())
            {
                say(2, "reading file");
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
			rounds++;
			solveMaster();
            break;
        case SOLVING_WORKER:
            rounds++;
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

void SimplestWorker::getMessage()
{
    Tag tag = receive();
    switch (tag)
    {    
    case INIT_DATA: 
        say(2, "received init message");
        setVariableMax(recvbuf[0]);
        state = IDLE;
        break;
	case LEARNT_DATA:
		say(2,"received learnt clause message");
		if (solver.okay() && !safe_addClause(recvbuf, LEARNT_CLAUSE))
		{
			result = false;
			state = DEAD;
		}
		break;
    case CLAUSE_DATA:
        say(2, "received clause message");
		if (solver.okay() && !safe_addClause(recvbuf, ORIGINAL_CLAUSE))
		{
			result = false;
			state = DEAD;
		}
		break;
    case TRAIL_DATA:		
        assumptions.clear();

        if (recvbuf[0] != master) 
        {
            master_changed++;
            say(1, "Node %d is now the master", recvbuf[0]);
        }
		master = recvbuf[0];

		if (recvbuf[1] != active_nodes)
		{
			say(1,"update on the number of active nodes from %d to %d",active_nodes,recvbuf[1]);
			active_nodes = recvbuf[1];

			conflicts.resize(active_nodes);
			models.resize(active_nodes);
		}
		
	    for (unsigned i = 2; i < recvbuf.size(); i++)
	        assumptions.push_back(recvbuf[i]);
	
		saybuffer(2,assumptions,"trail: ");

        if (isMaster())
        {            
            rounds++; // count the round, even though the master does nothing.
            state = IDLE;
        }
        else
            state = SOLVING_WORKER;
		if (!active_status)
			say(0,"active at %.3f seconds.",(double)((clock() - startTime)/(double)CLOCKS_PER_SEC));
		active_status = true;
        break;
    case MODEL_DATA:
        say(2, "receiving model from %d, me = %d, master = %d", status.MPI_SOURCE, rank, master);
        assert(isMaster());
        models[status.MPI_SOURCE] = recvbuf;
        numModels++;
        if ((numModels+numConflicts) >= (unsigned)active_nodes-1)
            conflictAnalysis();
        break;
    case CONFLICT_DATA:
        say(2, "receiving conflict from %d, me = %d, master = %d", status.MPI_SOURCE, rank, master);
        assert(isMaster());
        conflicts[status.MPI_SOURCE] = recvbuf;
        if (conflicts[status.MPI_SOURCE].size() == 0)
            result = false;
        numConflicts++;
        if ((numModels+numConflicts) >= (unsigned)active_nodes-1)
            conflictAnalysis();        
        break;
	case MEMORY_DATA:
		say(1,"receiving out of memory from %d", status.MPI_SOURCE);
		assert(isMaster());

		if (active_nodes == total_nodes)
		{
			say(0,"MEMORY EXHAUSTED");
			MPI_Abort(MPI_COMM_WORLD, S_ERROR);
		} 
		else 
		{
			
			active_nodes++;

			sendbuf.clear();
			sendbuf.push_back(rank);
			sendbuf.push_back(active_nodes);
			for (unsigned i = 0; i < assumptions.size() ; i++)
				sendbuf.push_back(assumptions[i]);
			safe_send(status.MPI_SOURCE, SPLIT_DATA);
			
			conflicts.resize(active_nodes);
			models.resize(active_nodes);
		}

		break;
	case SPLIT_DATA:
		say(2,"receiving split message from %d to send half of the clauses to %d", status.MPI_SOURCE, recvbuf[1]-1);
		assert(!isMaster());
		splitWorker();
		break;
    case TERMINATE_DATA:
        result = (recvbuf[0] != 0);
        state = DEAD;
        break;
    default:
        std::exception("unknown mpi tag");
    }
	
}

void SimplestWorker::safe_send(int to, int tag)
{  
    if (MPI_Send((void*) sendbuf.data(), sendbuf.size(), MPI_INT, to, tag, MPI_COMM_WORLD)!=MPI_SUCCESS)
		throw std::exception("MPI send failure");    
}

SimplestWorker::Tag SimplestWorker::receive()
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