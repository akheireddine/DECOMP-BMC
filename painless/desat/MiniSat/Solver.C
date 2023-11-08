/****************************************************************************************[Solver.C]
MiniSat -- Copyright (c) 2003-2005, Niklas Een, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/
#include <iostream>

#include "Solver.h"
#include "Sort.h"
#include <cmath>
#include <iostream>

#define CC_MINIMIZATION

using namespace Desat;

//=================================================================================================
// Helper functions:

void removeWatch(vec<Clause *> &ws, Clause *elem)
{
    if (ws.size() == 0)
        return; // (skip lists that are already cleared)
    int j = 0;
    for (; ws[j] != elem; j++)
    {
        assert(j < ws.size());
    }
    for (; j < ws.size() - 1; j++)
        ws[j] = ws[j + 1];
    ws.pop();
}

//=================================================================================================
// Operations on clauses:

/*_________________________________________________________________________________________________
|
|  newClause : (ps : const vec<Lit>&) (learnt : bool) (id : ClauseId)  ->  [void]
|
|  Description:
|    Allocate and add a new clause to the SAT solvers clause database. If a conflict is detected,
|    the 'ok' flag is cleared and the solver is in an unusable state (must be disposed).
|
|  Input:
|    ps     - The new clause as a vector of literals.
|    learnt - Is the clause a learnt clause? For learnt clauses, 'ps[0]' is assumed to be the
|             asserting literal. An appropriate 'enqueue()' operation will be performed on this
|             literal. One of the watches will always be on this literal, the other will be set to
|             the literal with the highest decision level.
|    id     - If logging proof, learnt clauses should be given an ID by caller.
|
|  Effect:
|    Activity heuristics are updated.
|________________________________________________________________________________________________@*/
void Solver::newClause(const vec<Lit> &ps_, bool learnt, ClauseId id, bool enqueue_asserting, bool halflearnt)
{
    assert(learnt || id == ClauseId_NULL);
    if (!ok)
        return;

    vec<Lit> qs;
    if (!learnt)
    {
        assert(decisionLevel() == 0);
        ps_.copyTo(qs); // Make a copy of the input vector.

        // Remove duplicates:
        sortUnique(qs);

        // Check if clause is satisfied:
        for (int i = 0; i < qs.size() - 1; i++)
        {
            if (qs[i] == ~qs[i + 1])
                return;
        }
        for (int i = 0; i < qs.size(); i++)
        {
            if (value(qs[i]) == D_True)
                return;
        }

        // Remove false literals:
        int i, j;
        if (proof != NULL)
            proof->beginChain(proof->addRoot(qs));
        for (i = j = 0; i < qs.size(); i++)
            if (value(qs[i]) != D_False)
                qs[j++] = qs[i];
            else
            {
                //std::cout << "removing " << (sign(qs[i])?"-":"") << var(qs[i]) << std::endl;
                // CMW: Need phase of pivot.
                //if (proof != NULL) proof->resolve(unit_id[var(qs[i])], var(qs[i]));
                if (proof != NULL)
                    proof->resolve(unit_id[var(qs[i])], qs[i]);
            }
        qs.shrink(i - j);
        if (proof != NULL)
            id = proof->endChain();
    }
    const vec<Lit> &ps = learnt ? ps_ : qs; // 'ps' is now the (possibly) reduced vector of literals.

    if (ps.size() == 0)
    {
        ok = false;
    }
    else if (ps.size() == 1)
    {
        //std::cout << "LEARNED UNIT: " << (sign(ps[0])?"-":"") << var(ps[0]) << std::endl;
        // NOTE: If enqueue takes place at root level, the assignment will be lost in incremental use (it doesn't seem to hurt much though).
        if (id != ClauseId_NULL)
            unit_id[var(ps[0])] = id;
        if (!enqueue(ps[0]))
            ok = false;
    }
    else
    {
        // Allocate clause:

#ifdef _MEMORY_ALLOCATION
        Clause *c = Clause_new(learnt, ps, memory_counter, memory_limit, id);
#else
        Clause *c = Clause_new(learnt, ps, id);
#endif
        if (learnt)
        {
            // Put the second watch on the literal with highest decision level:
            int max_i = 1;
            int max = level[var(ps[1])];
            for (int i = 2; i < ps.size(); i++)
                if (level[var(ps[i])] > max)
                    max = level[var(ps[i])],
                    max_i = i;
            (*c)[1] = ps[max_i];
            (*c)[max_i] = ps[1];

            // Bumping:
            claBumpActivity(c); // (newly learnt clauses should be considered active)

            // Enqueue asserting literal:
            if (enqueue_asserting)
                check(enqueue((*c)[0], c));

            // Store clause:
            watches[index(~(*c)[0])].push(c);
            watches[index(~(*c)[1])].push(c);

            if (halflearnt)
                halflearnts.push(c);
            else
                learnts.push(c);
            stats.learnts_literals += c->size();
        }
        else
        {
            // Store clause:
            watches[index(~(*c)[0])].push(c);
            watches[index(~(*c)[1])].push(c);
            clauses.push(c);
            stats.clauses_literals += c->size();
        }
    }
}

void Solver::importClauses()
{
    assert(decisionLevel() == 0);
    if (importClauseCallback == NULL)
        return;

    std::vector<int> importedClause;
    int sender = 0;
    while (importClauseCallback(issuer, importedClause, sender))
    {
        vec<Lit> cls_lit;
        for (unsigned i = 0; i < importedClause.size(); i++)
        {
            bool sign = importedClause[i] < 0;
            while(abs(importedClause[i]) >= nVars())
                newVar();
            Lit l = sign ? Lit(-importedClause[i], true) : Lit(importedClause[i], false);
            cls_lit.push(l);
        }
        addHalfLearntClause(cls_lit);
        importedClause.clear();
    }
}

// Disposes a clauses and removes it from watcher lists. NOTE! Low-level; does NOT change the 'clauses' and 'learnts' vector.
//
void Solver::remove(Clause *c, bool just_dealloc)
{

#ifdef _MEMORY_ALLOCATION

    uint decrease_memory = sizeof(Clause) + sizeof(uint) * (c->size() + (int)(c->learnt()));
    if (proof != NULL)
        decrease_memory++;
    assert(memory_counter >= decrease_memory);
    memory_counter -= decrease_memory;

#endif

    if (!just_dealloc)
    {
        removeWatch(watches[index(~(*c)[0])], c),
            removeWatch(watches[index(~(*c)[1])], c);

        if (c->learnt())
            stats.learnts_literals -= c->size();
        else
            stats.clauses_literals -= c->size();

        if (proof != NULL)
            proof->deleted(c->id());
    }

    xfree(c);
}

// Can assume everything has been propagated! (esp. the first two literals are != D_False, unless
// the clause is binary and satisfied, in which case the first literal is true)
// Returns True if clause is satisfied (will be removed), False otherwise.
//
bool Solver::simplify(Clause *c) const
{
    assert(decisionLevel() == 0);
    for (int i = 0; i < c->size(); i++)
    {
        if (value((*c)[i]) == D_True)
            return true;
    }
    return false;
}

//=================================================================================================
// Minor methods:

// Creates a new SAT variable in the solver. If 'decision_var' is cleared, variable will not be
// used as a decision variable (NOTE! This has effects on the meaning of a SATISFIABLE result).
//
Var Solver::newVar()
{
    int index;
    index = nVars();
    watches.push(); // (list for positive literal)
    watches.push(); // (list for negative literal)
    reason.push(NULL);
    assigns.push(toInt(D_Undef));
    varDecision.push(-1);
    level.push(-1);
    trail_pos.push(-1);
    activity.push(0);
    order.newVar();

    analyze_seen.push(0);
    seen2.push(0);

    polarity.push(true); //RM: polarity controls the phase of the branching variable
    if (proof != NULL)
        unit_id.push(ClauseId_NULL);
    return index;
}

// Returns FALSE if immediate conflict.
bool Solver::assume(Lit p)
{
    trail_lim.push(trail.size());
    return enqueue(p);
}

// Revert to the state at given level.
void Solver::cancelUntil(int level)
{
    if (decisionLevel() > level)
    {
        for (int c = trail.size() - 1; c >= trail_lim[level]; c--)
        {
            Var x = var(trail[c]);
            polarity[x] = sign(trail[c]); //RM: progress saving heuristic
            assigns[x] = toInt(D_Undef);
            reason[x] = NULL;
            order.undo(x);
        }
        trail.shrink(trail.size() - trail_lim[level]);
        trail_lim.shrink(trail_lim.size() - level);
        qhead = trail.size();
    }
}

//=================================================================================================
// Major methods:

/*_________________________________________________________________________________________________
|
|  analyze : (confl : Clause*) (out_learnt : vec<Lit>&) (out_btlevel : int&)  ->  [void]
|
|  Description:
|    Analyze conflict and produce a reason clause ('out_learnt') and a backtracking level
|    ('out_btlevel').
|
|    Pre-conditions:
|      * 'out_learnt' is assumed to be cleared.
|      * Current decision level must be greater than root level.
|
|    Post-conditions:
|      * 'out_learnt[0]' is the asserting literal at level 'out_btlevel'.
|      * If performing proof-logging, the last derived clause in the proof is the reason clause.
|________________________________________________________________________________________________@*/

class lastToFirst_lt
{ // Helper class to 'analyze' -- order literals from last to first occurance in 'trail[]'.
    const vec<int> &trail_pos;

public:
    lastToFirst_lt(const vec<int> &t) : trail_pos(t) {}
    bool operator()(Lit p, Lit q) { return trail_pos[var(p)] > trail_pos[var(q)]; }
};

void Solver::analyze(Clause *confl, vec<Lit> &out_learnt, int &out_btlevel, ClauseId &out_id, int &out_lbd)
{
    vec<char> &seen = analyze_seen;
    int pathC = 0;
    Lit p = lit_Undef;
    out_id = ClauseId_NULL;

    // Generate conflict clause:
    //
    if (proof != NULL)
        proof->beginChain(confl->id());
    out_learnt.push(); // (leave room for the asserting literal)
    out_btlevel = 0;
    int index = trail.size() - 1;
    for (;;)
    {
        assert(confl != NULL); // (otherwise should be UIP)

        Clause &c = *confl;
        if (c.learnt())
            claBumpActivity(&c);

        for (int j = (p == lit_Undef) ? 0 : 1; j < c.size(); j++)
        {
            Lit q = c[j];
            if (!seen[var(q)])
            {
                if (level[var(q)] > 0)
                {
                    varBumpActivity(q);
                    seen[var(q)] = 1;
                    if (level[var(q)] == decisionLevel())
                        pathC++;
                    else
                    {
                        out_learnt.push(q);
                        out_btlevel = max_minisat(out_btlevel, level[var(q)]);
                    }
                }
                else
                    // CMW: Need phase of pivot.
                    //if (proof != NULL) proof->resolve(unit_id[var(q)], var(q));
                    if (proof != NULL)
                    proof->resolve(unit_id[var(q)], q);
            }
        }

        // Select next clause to look at:
        while (!seen[var(trail[index--])])
            ;
        p = trail[index + 1];
        confl = reason[var(p)];
        seen[var(p)] = 0;
        pathC--;
        if (pathC == 0)
            break;

        // CMW: Need phase of pivot.
        //if (proof != NULL) proof->resolve(confl->id(), var(p));
        if (proof != NULL)
            proof->resolve(confl->id(), ~p);
    }
    out_learnt[0] = ~p;

    // Conflict clause minimization:
    //
#ifdef CC_MINIMIZATION
    int i, j;
    if (expensive_ccmin)
    {
        // Simplify conflict clause (a lot):
        //
        uint min_level = 0;
        for (i = 1; i < out_learnt.size(); i++)
            min_level |= 1 << (level[var(out_learnt[i])] & 31); // (maintain an abstraction of levels involved in conflict)

        analyze_toclear.clear();
        for (i = j = 1; i < out_learnt.size(); i++)
            if (reason[var(out_learnt[i])] == NULL || !analyze_removable(out_learnt[i], min_level))
                out_learnt[j++] = out_learnt[i];
    }
    else
    {
        // Simplify conflict clause (a little):
        //
        analyze_toclear.clear();
        for (i = j = 1; i < out_learnt.size(); i++)
        {
            Clause *r = reason[var(out_learnt[i])];
            if (r == NULL)
                out_learnt[j++] = out_learnt[i];
            else
            {
                Clause &c = *r;
                for (int k = 1; k < c.size(); k++)
                    if (!seen[var(c[k])] && level[var(c[k])] != 0)
                    {
                        out_learnt[j++] = out_learnt[i];
                        goto Keep;
                    }
                analyze_toclear.push(out_learnt[i]);
            Keep:;
            }
        }
    }
#else
    int i = 0, j = 0;
#endif

    // Finilize proof logging with conflict clause minimization steps:
    //
    if (proof != NULL)
    {
        sort(analyze_toclear, lastToFirst_lt(trail_pos));
        for (int k = 0; k < analyze_toclear.size(); k++)
        {
            Var v = var(analyze_toclear[k]);
            assert(level[v] > 0);
            Clause &c = *reason[v];
            // CMW: Need phase of pivot
            // proof->resolve(c.id(), v);
            proof->resolve(c.id(), analyze_toclear[k]);
            for (int l = 1; l < c.size(); l++)
                if (level[var(c[l])] == 0)
                    // CMW: Need phase of pivot
                    // proof->resolve(unit_id[var(c[k])], var(c[k]));
                    proof->resolve(unit_id[var(c[l])], c[l]);
        }
        out_id = proof->endChain();
    }
    // Clean up:
    //
    for (int j = 0; j < out_learnt.size(); j++)
        seen[var(out_learnt[j])] = 0;
    for (int j = 0; j < analyze_toclear.size(); j++)
        seen[var(analyze_toclear[j])] = 0; // ('seen[]' is now cleared)

    stats.max_literals += out_learnt.size();
    out_learnt.shrink(i - j);
    stats.tot_literals += out_learnt.size();

    out_lbd = computeLBD(out_learnt);
}

// Check if 'p' can be removed. 'min_level' is used to abort early if visiting literals at a level that cannot be removed.
//
bool Solver::analyze_removable(Lit p, uint min_level)
{
    assert(reason[var(p)] != NULL);
    analyze_stack.clear();
    analyze_stack.push(p);
    int top = analyze_toclear.size();
    while (analyze_stack.size() > 0)
    {
        assert(reason[var(analyze_stack.last())] != NULL);
        Clause &c = *reason[var(analyze_stack.last())];
        analyze_stack.pop();
        for (int i = 1; i < c.size(); i++)
        {
            Lit p = c[i];
            if (!analyze_seen[var(p)] && level[var(p)] != 0)
            {
                if (reason[var(p)] != NULL && ((1 << (level[var(p)] & 31)) & min_level) != 0)
                {
                    analyze_seen[var(p)] = 1;
                    analyze_stack.push(p);
                    analyze_toclear.push(p);
                }
                else
                {
                    for (int j = top; j < analyze_toclear.size(); j++)
                        analyze_seen[var(analyze_toclear[j])] = 0;
                    analyze_toclear.shrink(analyze_toclear.size() - top);
                    return false;
                }
            }
        }
    }
    analyze_toclear.push(p);

    return true;
}

/*_________________________________________________________________________________________________
|
|  analyzeFinal : (confl : Clause*) (skip_first : bool)  ->  [void]
|
|  Description:
|    Specialized analysis procedure to express the final conflict in terms of assumptions.
|    'root_level' is allowed to point beyond end of trace (useful if called after conflict while
|    making assumptions). If 'skip_first' is TRUE, the first literal of 'confl' is  ignored (needed
|    if conflict arose before search even started).
|________________________________________________________________________________________________@*/
void Solver::analyzeFinal(Clause *confl, bool skip_first)
{
    // -- NOTE! This code is relatively untested. Please report bugs!
    conflict.clear();
    if (root_level == 0)
    {
        if (proof != NULL)
            conflict_id = proof->last();
        return;
    }

    vec<char> &seen = analyze_seen;
    if (proof != NULL)
        proof->beginChain(confl->id());
    for (int i = skip_first ? 1 : 0; i < confl->size(); i++)
    {
        Var x = var((*confl)[i]);
        if (level[x] > 0)
            seen[x] = 1;
        else
            // CMW: Need phase of pivot.
            //if (proof != NULL) proof->resolve(unit_id[x], x);
            if (proof != NULL)
            proof->resolve(unit_id[x], (*confl)[i]);
    }

    int start = (root_level >= trail_lim.size()) ? trail.size() - 1 : trail_lim[root_level];
    for (int i = start; i >= trail_lim[0]; i--)
    {
        Var x = var(trail[i]);
        if (seen[x])
        {
            Clause *r = reason[x];
            if (r == NULL)
            {
                assert(level[x] > 0);
                conflict.push(~trail[i]);
            }
            else
            {
                Clause &c = *r;
                // CMW: Need phase of pivot.
                //if (proof != NULL) proof->resolve(c.id(), x);
                if (proof != NULL)
                    proof->resolve(c.id(), ~trail[i]);
                for (int j = 1; j < c.size(); j++)
                    if (level[var(c[j])] > 0)
                        seen[var(c[j])] = 1;
                    else
                        // CMW: Need phase of pivot.
                        // if (proof != NULL) proof->resolve(unit_id[var(c[j])], var(c[j]));
                        if (proof != NULL)
                        proof->resolve(unit_id[var(c[j])], c[j]);
            }
            seen[x] = 0;
        }
    }
    if (proof != NULL)
        conflict_id = proof->endChain();
}

/*_________________________________________________________________________________________________
|
|  enqueue : (p : Lit) (from : Clause*)  ->  [bool]
|
|  Description:
|    Puts a new fact on the propagation queue as well as immediately updating the variable's value.
|    Should a conflict arise, FALSE is returned.
|
|  Input:
|    p    - The fact to enqueue
|    from - [Optional] Fact propagated from this (currently) unit clause. Stored in 'reason[]'.
|           Default value is NULL (no reason).
|
|  Output:
|    TRUE if fact was enqueued without conflict, FALSE otherwise.
|________________________________________________________________________________________________@*/
bool Solver::enqueue(Lit p, Clause *from)
{
    if (value(p) != D_Undef)
        return value(p) != D_False;
    else
    {
        Var x = var(p);
        assigns[x] = toInt(lbool(!sign(p)));
        level[x] = decisionLevel();
        trail_pos[x] = trail.size();
        reason[x] = from;
        trail.push(p);
        varDecision[x] = decisionLevel();
        return true;
    }
}

/*_________________________________________________________________________________________________
|
|  propagate : [void]  ->  [Clause*]
|
|  Description:
|    Propagates all enqueued facts. If a conflict arises, the conflicting clause is returned,
|    otherwise NULL. NOTE! This method has been optimized for speed rather than readability.
|
|    Post-conditions:
|      * The propagation queue is empty, even if there was a conflict.
|________________________________________________________________________________________________@*/
Clause *Solver::propagate()
{
    Clause *confl = NULL;
    while (qhead < trail.size())
    {
        stats.propagations++;
        simpDB_props--;

        Lit p = trail[qhead++]; // 'p' is enqueued fact to propagate.
        vec<Clause *> &ws = watches[index(p)];
        Clause **i, **j, **end;

        for (i = j = (Clause **)ws, end = i + ws.size(); i != end;)
        {
            Clause &c = **i;
            i++;
            // Make sure the false literal is data[1]:
            Lit false_lit = ~p;
            if (c[0] == false_lit)
                c[0] = c[1], c[1] = false_lit;

            assert(c[1] == false_lit);

            // If 0th watch is true, then clause is already satisfied.
            Lit first = c[0];
            lbool val = value(first);
            if (val == D_True)
            {
                *j++ = &c;
            }
            else
            {
                // Look for new watch:
                for (int k = 2; k < c.size(); k++)
                    if (value(c[k]) != D_False)
                    {
                        c[1] = c[k];
                        c[k] = false_lit;
                        watches[index(~c[1])].push(&c);
                        goto FoundWatch;
                    }

                // Did not find watch -- clause is unit under assignment:
                if (decisionLevel() == 0 && proof != NULL)
                {
                    // Log the production of this unit clause:
                    proof->beginChain(c.id());
                    for (int k = 1; k < c.size(); k++)
                    {
                        // CMW: Need phase of pivot.
                        // proof->resolve(unit_id[var(c[k])], var(c[k]));
                        Lit l = c[k];
                        int v = var(l);
                        ClauseId uid = unit_id[v];
                        proof->resolve(uid, l);
                    }
                    ClauseId id = proof->endChain();
                    assert(unit_id[var(first)] == ClauseId_NULL || value(first) == D_False); // (if variable already has 'id', it must be with the other polarity and we should have derived the empty clause here)
                    if (value(first) != D_False)
                        unit_id[var(first)] = id;
                    else
                    {
                        // Empty clause derived:
                        proof->beginChain(unit_id[var(first)]);
                        // CMW: Need phase of pivot.
                        //proof->resolve(id, var(first));
                        proof->resolve(id, ~first);
                        proof->endChain();
                    }
                }

                *j++ = &c;
                if (!enqueue(first, &c))
                {
                    if (decisionLevel() == 0)
                        ok = false;
                    confl = &c;
                    qhead = trail.size();
                    // Copy the remaining watches:
                    while (i < end)
                        *j++ = *i++;
                }
            FoundWatch:;
            }
        }
        ws.shrink((int)(i - j));
    }

    return confl;
}

/*_________________________________________________________________________________________________
|
|  reduceDB : ()  ->  [void]
|
|  Description:
|    Remove half of the learnt clauses, minus the clauses locked by the current assignment. Locked
|    clauses are clauses that are reason to some assignment. Binary clauses are never removed.
|________________________________________________________________________________________________@*/
struct reduceDB_lt
{
    bool operator()(Clause *x, Clause *y) { return x->size() > 2 && (y->size() == 2 || x->activity() < y->activity()); }
};
void Solver::reduceDB()
{
    int i, j;
    double extra_lim = cla_inc / learnts.size(); // Remove any clause below this activity

    sort(learnts, reduceDB_lt());
    for (i = j = 0; i < learnts.size() / 2; i++)
    {
        if (learnts[i]->size() > 2 && !locked(learnts[i]))
            remove(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    for (; i < learnts.size(); i++)
    {
        if (learnts[i]->size() > 2 && !locked(learnts[i]) && learnts[i]->activity() < extra_lim)
            remove(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    learnts.shrink(i - j);

    reduce_callback();
}

void Solver::splitDB(vec<vec<Lit>> &split_clauses)
{
    int j;

    if (clauses.size() % 2 == 0)
        j = clauses.size() / 2;
    else
        j = (clauses.size() - 1) / 2;

    int removed = 0;

    while (removed < j)
    {
        int last = clauses.size() - 1;

        vec<Lit> read_clause;
        for (int z = 0; z < clauses[last]->size(); z++)
        {
            Clause &c = *clauses[last];
            Lit p = c[z];
            //if(sign(p)) std::cout << -var(p) << " ";
            //else std::cout << var(p) << " ";
            read_clause.push(p);
        }
        //std::cout << std::endl;

        split_clauses.push(read_clause);
        //std::cout << "removing clause: "<< removed << " remaining: " << clauses.size() << " goal: " << j << std::endl;
        remove(clauses[last]);

        clauses.pop();

        removed++;
    }
}

void Solver::getAndRemovePosition(std::vector<signed> &last_clause, int position, int type)
{
    //type : ORIGINAL 0; LEARNT 1; HALFLEARNT 2
    cancelUntil(0);

    last_clause.clear();

    int last = position;
    int size = 0;

    switch (type)
    {
    case 0:
        if (clauses.size() > 0)
        {
            size = clauses[last]->size();
            Clause &c = *clauses[last];

            for (int i = 0; i < size; i++)
            {
                Lit p = c[i];
                if (sign(p))
                    last_clause.push_back(-var(p));
                else
                    last_clause.push_back(var(p));
            }
            remove(clauses[last]);

            for (int i = last; i < clauses.size() - 1; i++)
                clauses[i] = clauses[i + 1];

            clauses.pop();
        }
        break;

    case 1:
        if (learnts.size() > 0)
        {
            last = position;
            while (last >= 0)
            {
                if (locked(learnts[last]))
                    last--;
                else
                    break;
            }

            if (last >= 0)
            {
                size = learnts[last]->size();
                Clause &c = *learnts[last];
                assert(!locked(learnts[last]));
                for (int i = 0; i < size; i++)
                {
                    Lit p = c[i];
                    if (sign(p))
                        last_clause.push_back(-var(p));
                    else
                        last_clause.push_back(var(p));
                }
                remove(learnts[last]);

                for (int i = last; i < learnts.size() - 1; i++)
                    learnts[i] = learnts[i + 1];

                learnts.pop();
            }

        case 2:
            if (halflearnts.size() > 0)
            {
                size = halflearnts[last]->size();
                Clause &c = *halflearnts[last];

                for (int i = 0; i < size; i++)
                {
                    Lit p = c[i];
                    if (sign(p))
                        last_clause.push_back(-var(p));
                    else
                        last_clause.push_back(var(p));
                }
                remove(halflearnts[last]);

                for (int i = last; i < halflearnts.size() - 1; i++)
                    halflearnts[i] = halflearnts[i + 1];

                halflearnts.pop();
            }
            break;

        default:
            throw std::runtime_error("type of clause not defined");
        }
    }
}

void Solver::getAndRemoveLast(std::vector<signed> &last_clause, int type)
{

    //type : ORIGINAL 0; LEARNT 1; HALFLEARNT 2
    cancelUntil(0);

    last_clause.clear();

    int last = 0;
    int size = 0;

    switch (type)
    {
    case 0: //ORIGINAL
        if (clauses.size() > 0)
        {
            last = clauses.size() - 1;
            size = clauses[last]->size();
            Clause &c = *clauses[last];
            for (int i = 0; i < size; i++)
            {
                Lit p = c[i];
                if (sign(p))
                    last_clause.push_back(-var(p));
                else
                    last_clause.push_back(var(p));
            }
            remove(clauses[last]);
            clauses.pop();
        }
        break;
    case 1: //LEARNT
        if (learnts.size() > 0)
        {
            last = learnts.size() - 1;
            while (last >= 0)
            {
                if (locked(learnts[last]))
                    last--;
                else
                    break;
            }

            if (last >= 0)
            {
                size = learnts[last]->size();
                Clause &c = *learnts[last];
                assert(!locked(learnts[last]));
                for (int i = 0; i < size; i++)
                {
                    Lit p = c[i];
                    if (sign(p))
                        last_clause.push_back(-var(p));
                    else
                        last_clause.push_back(var(p));
                }
                remove(learnts[last]);

                if (last != learnts.size() - 1)
                {
                    int i, j;
                    for (i = j = last; i < learnts.size(); i++)
                        learnts[j++] = learnts[i];

                    learnts.shrink(i - j);
                }
                else
                    learnts.pop();
            }
        }
        break;
    case 2: //HALFLEARNT
        if (halflearnts.size() > 0)
        {
            last = halflearnts.size() - 1;
            size = halflearnts[last]->size();
            Clause &c = *halflearnts[last];
            for (int i = 0; i < size; i++)
            {
                Lit p = c[i];
                if (sign(p))
                    last_clause.push_back(-var(p));
                else
                    last_clause.push_back(var(p));
            }
            remove(halflearnts[last]);
            halflearnts.pop();
        }
        break;
    default:
        throw std::runtime_error("type of clause not defined");
    }
}

void Solver::removeLast(vec<vec<Lit>> &split_clauses)
{
    int last = clauses.size() - 1;

    vec<Lit> read_clause;
    for (int z = 0; z < clauses[last]->size(); z++)
    {
        Clause &c = *clauses[last];
        Lit p = c[z];
        //if(sign(p)) std::cout << -var(p) << " ";
        //else std::cout << var(p) << " ";
        read_clause.push(p);
    }
    //std::cout << std::endl;

    split_clauses.push(read_clause);
    //std::cout << "removing clause: "<< removed << " remaining: " << clauses.size() << " goal: " << j << std::endl;
    remove(clauses[last]);

    clauses.pop();
}

void Solver::printDB(int rank)
{
    for (int i = 0; i < clauses.size(); i++)
    {
        Clause &c = *clauses[i];
        std::cout << "[ " << rank << " ] ";
        for (int z = 0; z < clauses[i]->size(); z++)
        {
            Lit p = c[z];
            if (sign(p))
                std::cout << -var(p) << " ";
            else
                std::cout << var(p) << " ";
        }
        std::cout << "0" << std::endl;
    }

    cancelUntil(0);
    for (int i = 1; i < nVars(); i++)
    {
        if (value(i) != D_Undef)
        {
            if (value(i) == D_True)
                std::cout << "c " << i << " 0" << std::endl;
            else
                std::cout << "c " << -i << " 0" << std::endl;
        }
    }
}

void Solver::getUnitClauses(std::set<int> &vars, int maxVar)
{
    cancelUntil(0);
    for (int i = 1; i <= maxVar; i++)
    {
        if (value(i) != D_Undef)
        {
            if (value(i) == D_True)
                vars.insert(i);
            else
                vars.insert(-i);
        }
    }
}

void Solver::getUsedVariables(std::set<int> &vars, int maxVar)
{
    for (int i = 0; i < clauses.size(); i++)
    {
        Clause &c = *clauses[i];
        for (int z = 0; z < clauses[i]->size(); z++)
        {
            Lit p = c[z];
            if (var(p) <= maxVar)
                vars.insert(var(p));
        }
    }
}

void Solver::cleanLearnts()
{

    int i, j;

    for (i = j = 0; i < learnts.size(); i++)
    {
        if (learnts[i]->size() > 2 && !locked(learnts[i]))
            remove(learnts[i]);
        else
            learnts[j++] = learnts[i];
    }

    learnts.shrink(i - j);

    //for(int i = 0; i < learnts.size(); i++)
    //	remove(learnts[i]);

    //learnts.clear();
}

/*_________________________________________________________________________________________________
|
|  simplifyDB : [void]  ->  [bool]
|
|  Description:
|    Simplify the clause database according to the current top-level assigment. Currently, the only
|    thing done here is the removal of satisfied clauses, but more things can be put here.
|________________________________________________________________________________________________@*/
void Solver::simplifyDB()
{
    if (!ok)
        return; // GUARD (public method)
    assert(decisionLevel() == 0);

    if (propagate() != NULL)
    {
        ok = false;
        return;
    }

    if (nAssigns() == simpDB_assigns || simpDB_props > 0) // (nothing has changed or preformed a simplification too recently)
        return;

    // Clear watcher lists:
    // CMW: New clauses over the units may be added through interpolation!
    for (int i = simpDB_assigns; i < nAssigns(); i++)
    {
        //for (int i = 0; i < nAssigns(); i++){
        Lit p = trail[i];
        watches[index(p)].clear(true);
        watches[index(~p)].clear(true);
    }

    // Remove satisfied clauses:
    for (int type = 0; type < 2; type++)
    {
        vec<Clause *> &cs = type ? learnts : clauses;
        int j = 0;
        for (int i = 0; i < cs.size(); i++)
        {
            if (!locked(cs[i]) && simplify(cs[i]))
                remove(cs[i]);
            else
                cs[j++] = cs[i];
        }
        cs.shrink(cs.size() - j);
    }

    simpDB_assigns = nAssigns();
    simpDB_props = stats.clauses_literals + stats.learnts_literals; // (shouldn't depend on 'stats' really, but it will do for now)
}

/*_________________________________________________________________________________________________
|
|  search : (nof_conflicts : int) (nof_learnts : int) (params : const SearchParams&)  ->  [lbool]
|
|  Description:
|    Search for a model the specified number of conflicts, keeping the number of learnt clauses
|    below the provided limit. NOTE! Use negative value for 'nof_conflicts' or 'nof_learnts' to
|    indicate infinity.
|
|  Output:
|    'D_True' if a partial assigment that is consistent with respect to the clauseset is found. If
|    all variables are decision variables, this means that the clause set is satisfiable. 'D_False'
|    if the clause set is unsatisfiable. 'D_Undef' if the bound on number of conflicts is reached.
|________________________________________________________________________________________________@*/
lbool Solver::search(int nof_conflicts, int nof_learnts, const SearchParams &params)
{

    if (!ok)
        return D_False; // GUARD (public method)
    assert(root_level == decisionLevel());

    stats.starts++;
    int conflictC = 0;
    var_decay = 1 / params.var_decay;
    cla_decay = 1 / params.clause_decay;
    model.clear();

    for (;;)
    {
        if (decisionLevel() == 0)
        {
            importClauses();
        }

        if (resources_limit)
        {
            double solving_time = 0;
            if (resources_time > 0)
                solving_time = (clock() - start) / (double)CLOCKS_PER_SEC;

            if ((resources_conflicts > 0 && stats.conflicts_limit > resources_conflicts) ||
                (resources_clauses > 0 && learnts.size() > resources_clauses) ||
                (resources_time > 0 && solving_time > resources_time))
            {
                resources_exhausted = true;
                cancelUntil(0);
                return D_False;
            }
        }

        Clause *confl = propagate();
        if (confl != NULL)
        {
            // CONFLICT

            stats.conflicts++;
            conflictC++;
            stats.conflicts_limit++;
            //vec<Lit>    learnt_clause;
            //int         backtrack_level;
            if (decisionLevel() == root_level)
            {
                // Contradiction found:
                analyzeFinal(confl);
                return D_False;
            }

            //analyze(confl, learnt_clause, backtrack_level);
            //cancelUntil(max_minisat(backtrack_level, root_level));
            //newClause(learnt_clause, true, (proof != NULL) ? proof->last() : ClauseId_NULL);
            //if (learnt_clause.size() == 1) level[var(learnt_clause[0])] = 0;    // (this is ugly (but needed for 'analyzeFinal()') -- in future versions, we will backtrack past the 'root_level' and redo the assumptions)

            // CMW: overridable conflict analysis
            if (!resolve_conflict(confl))
                return D_False;

            varDecayActivity();
            claDecayActivity();
        }
        else
        {
            // NO CONFLICT

            if (nof_conflicts >= 0 && conflictC >= nof_conflicts)
            {
                // Reached bound on number of conflicts:
                progress_estimate = progressEstimate();
                cancelUntil(root_level);
                return D_Undef;
            }

            if (decisionLevel() == 0)
                // Simplify the set of problem clauses:
                simplifyDB(), assert(ok);

            if (nof_learnts >= 0 && learnts.size() - nAssigns() >= nof_learnts)
                // Reduce the set of learnt clauses:
                reduceDB();

            if (partial_models)
            {
                if (checkSatisfaction())
                {
                    // Model found:
                    model.growTo(nVars());
                    for (int i = 0; i < nVars(); i++)
                        model[i] = value(i);
                    cancelUntil(root_level);
                    return D_True;
                }

                // New variable decision:
                stats.decisions++;
                Var next = order.select(0.0); //params.random_var_freq);
                //RM 16/05/2012: phase saving heuristic when choosing a new variable decision
                (!phase_saving || polarity[next]) ? check(next != var_Undef && assume(~Lit(next))) : check(next != var_Undef && assume(Lit(next)));
            }
            else
            {
                // New variable decision:
                stats.decisions++;
                Var next = order.select(0.0); //params.random_var_freq);

                if (next == var_Undef)
                {
                    // Model found:
                    model.growTo(nVars());
                    for (int i = 0; i < nVars(); i++)
                        model[i] = value(i);
                    cancelUntil(root_level);
                    return D_True;
                }
                //RM: phase saving heuristic when choosing a new variable decision
                (!phase_saving || polarity[next]) ? check(assume(~Lit(next))) : check(assume(Lit(next)));
            }
        }
    }
}

// Return search-space coverage. Not extremely reliable.
//
double Solver::progressEstimate()
{
    double progress = 0;
    double F = 1.0 / nVars();
    for (int i = 0; i < nVars(); i++)
        if (value(i) != D_Undef)
            progress += pow(F, level[i]);
    return progress / nVars();
}

// Divide all variable activities by 1e100.
//
void Solver::varRescaleActivity()
{
    for (int i = 0; i < nVars(); i++)
        activity[i] *= 1e-100;
    var_inc *= 1e-100;
}

// Divide all constraint activities by 1e100.
//
void Solver::claRescaleActivity()
{
    for (int i = 0; i < learnts.size(); i++)
        learnts[i]->activity() *= (float)1e-20;
    cla_inc *= 1e-20;
}

/*_________________________________________________________________________________________________
|
|  solve : (assumps : const vec<Lit>&)  ->  [bool]
|
|  Description:
|    Top-level solve. If using assumptions (non-empty 'assumps' vector), you must call
|    'simplifyDB()' first to see that no top-level conflict is present (which would put the solver
|    in an undefined state).
|
|  Input:
|    A list of assumptions (unit clauses coded as literals). Pre-condition: The assumptions must
|    not contain both 'x' and '~x' for any variable 'x'.
|________________________________________________________________________________________________@*/
bool Solver::solve(const vec<Lit> &assumps)
{

    //std::cout << "Assumptions inside the solver: ";
    //for(int i=0; i < assumps.size(); i++){
    //	Lit l = assumps[i];
    //	int v = var(l);
    //	if(sign(l)) std::cout << -v << " ";
    //	else std::cout << v << " ";
    //}
    //std::cout << std::endl;
    if (asynch_interrupt)
        return false;

    start = clock();
    stats.conflicts_limit = 0;

    cancelUntil(0);

    simplifyDB();
    if (!ok)
        return false;

    SearchParams params(default_params);
    double nof_conflicts = 100;
    double nof_learnts = nClauses() / 3;
    lbool status = D_Undef;

    double inner = 100;
    double outer = 100;

    // Perform assumptions:
    root_level = assumps.size();
    for (int i = 0; i < assumps.size(); i++)
    {
        Lit p = assumps[i];
        assert(var(p) < nVars());
        if (!assume(p))
        {
            if (reason[var(p)] != NULL)
            {
                analyzeFinal(reason[var(p)], true);
                conflict.push(~p);
            }
            else
            {
                assert(proof == NULL || unit_id[var(p)] != ClauseId_NULL); // (this is the pre-condition above)
                conflict.clear();
                conflict.push(~p);
                if (proof != NULL)
                    conflict_id = unit_id[var(p)];
            }
            cancelUntil(0);
            return false;
        }
        Clause *confl = propagate();
        if (confl != NULL)
        {
            analyzeFinal(confl), assert(conflict.size() > 0);
            cancelUntil(0);
            return false;
        }
    }
    assert(root_level == decisionLevel());

    // Search:
    if (verbosity >= 1)
    {
        reportf("==================================[MINISAT]===================================\n");
        reportf("| Conflicts |     ORIGINAL     |              LEARNT              | Progress |\n");
        reportf("|           | Clauses Literals |   Limit Clauses Literals  Lit/Cl |          |\n");
        reportf("==============================================================================\n");
    }

    // for (int i=0; i<clauses.size(); i++)
    // {
    //   const Clause &c = *clauses[i];
    //	std::cout << "( ";
    //    for (int j=0; j<c.size(); j++){
    //	 const Lit &l = c[j];
    //		int v = var(l);
    //		if(sign(l)) std::cout << -v << " ";
    //		else std::cout << v << " ";
    //}
    // std::cout << ")" << std::endl;
    // }

    while (status == D_Undef && !asynch_interrupt)
    {
        if (verbosity >= 1)
        {
            reportf("| %9d | %7d %8d | %7d %7d %8d %7.1f | %6.3f %% |\n", (int)stats.conflicts, nClauses(), (int)stats.clauses_literals, (int)nof_learnts, nLearnts(), (int)stats.learnts_literals, (double)stats.learnts_literals / nLearnts(), progress_estimate * 100);
            fflush(stdout);
        }
        status = search((int)nof_conflicts, (int)nof_learnts, params);

        //default restart strategy
        if (!fast_restart)
        {
            nof_conflicts *= 1.5;
            nof_learnts *= 1.1;
        }
        else
        {

            //armin restart strategy
            if (inner >= outer)
            {
                outer *= 1.1;
                inner = 100;
            }
            else
            {
                inner *= 1.1;
            }
            nof_conflicts = ceil(inner);
        }
    }
    if (verbosity >= 1)
        reportf("==============================================================================\n");

    cancelUntil(0);

    if (model_lifting && status == D_True)
        liftModel(assumps);

    return status == D_True;
}

bool Solver::resolve_conflict(Clause *confl)
{
    vec<Lit> learnt_clause;
    int backtrack_level;
    int lbd;
    ClauseId c_id = ClauseId_NULL;
    analyze(confl, learnt_clause, backtrack_level, c_id, lbd);

    if (exportClauseCallback != NULL)
    {
        std::vector<int> learnt_clause_int;
        for (unsigned i = 0; i < learnt_clause.size(); i++)
        {
            int variable = sign(learnt_clause[i]) ? -var(learnt_clause[i]) : var(learnt_clause[i]);
            learnt_clause_int.emplace_back(variable);
        }
        exportClauseCallback(issuer, lbd, learnt_clause_int);
    }

    assert(!proof || c_id != ClauseId_NULL);
    //std::cout << "BTL: " << backtrack_level << std::endl;
    cancelUntil(max_minisat(backtrack_level, root_level));
    newClause(learnt_clause, true, c_id);
    //std::cout << "ASSERTING: " << var(learnt_clause[0]) << std::endl;
    if (learnt_clause.size() == 1)
        level[var(learnt_clause[0])] = 0; // (this is ugly (but needed for 'analyzeFinal()') -- in future versions, we will backtrack past the 'root_level' and redo the assumptions)
    return true;
}

void Solver::liftModel(const vec<Lit> &assumps)
{
    //std::cout << "Lifting model " << std::endl;
    //   unsigned lifted = 0;
    vec<bool> requiredVars(nVars(), false);
    bool progress = true;

    for (int i = 0; i < nVars(); i++)
        if (level[i] == 0)
            requiredVars[var(trail[i])] = true;

    while (progress)
    {
        //std::cout << "Lifting ... " << std::endl;
        progress = false;

        for (int i = 0; i < clauses.size(); i++)
        {
            const Clause &c = *clauses[i];
            unsigned trueCnt = 0;
            Lit lastTrue;

            for (int j = 0; j < c.size(); j++)
            {
                const Lit &l = c[j];
                unsigned v = var(l);
                bool sgn = sign(l);

                if ((!sgn && model[v] == D_True) ||
                    (sgn && model[v] == D_False))
                {
                    if (requiredVars[v])
                        goto nextClause;
                    else
                    {
                        trueCnt++;
                        lastTrue = l;
                    }
                }

                requiredVars[var(lastTrue)] = true;
                progress = true;
            }

        nextClause:;
        }
    }

    for (int i = 0; i < assumps.size(); i++)
        requiredVars[var(assumps[i])] = true;

    for (int i = 0; i < model.size(); i++)
        if (!requiredVars[i] && level[i] > root_level)
        {
            // std::cout << "Not required: " << i << std::endl;
            model[i] = D_Undef;
        }
}

bool Solver::checkSatisfaction()
{
    bool sat;

    for (int i = 0; i < clauses.size(); i++)
    {
        const Clause &c = *clauses[i];
        sat = false;
        for (int j = 0; j < c.size(); j++)
            if (value(c[j]) == D_True)
            {
                sat = true;
                break;
            }
        if (!sat)
            return false;
    }

    return true;
}

/*
void Solver::findWatch(const vec<Clause*>& ws, const Clause* elem) const
{
  if (ws.size()==0)
  {
    std::cout << "CLAUSE NOT FOUND" << std::endl;
    for (int i=0; i<elem->size(); i++)
    {
      lbool v = value((*elem)[i]);
      std::cout << " " << toDimacs((*elem)[i]) << "(" << ((v==D_True)?1:(v==D_Undef)?2:0) << ")";
    }
    std::cout << std::endl;
    throw std::runtime_error("Watchlist empty!");
  }

  int j = 0;
  for (; ws[j] != elem  ; j++)
  {
    if (j==ws.size())
    {
      std::cout << "CLAUSE NOT FOUND" << std::endl;
        for (int i=0; i<elem->size(); i++)
        {
          lbool v = value((*elem)[i]);
          std::cout << " " << toDimacs((*elem)[i]) << "(" << ((v==D_True)?1:(v==D_Undef)?2:0) << ")";
        }
        std::cout << std::endl;
      throw std::runtime_error("Clause not found");
    }
  }
}

void Solver::checkWatches(const vec<Clause *> &cl) const
{
  for (int i=0; i<cl.size(); i++)
  {
    const Clause &c = *cl[i];
    if (level[var(c[0])]!=0) findWatch(watches[index(~c[0])], &c);
    if (level[var(c[1])]!=0) findWatch(watches[index(~c[1])], &c);
  }
}
*/
