#pragma once

#include <vector>
#include <streambuf>
#include <iostream>
#include <zlib.h>
#include <set>
#include <string>



//-------------------------------------------------------------------------------------------------
// A simple buffered character stream class:

class StreamBuffer_ {
    gzFile        in;
    unsigned char buf[1048576];
    int           pos;
    int           size;

    void assureLookahead() {
        if (pos >= size) {
            pos  = 0;
            size = gzread(in, buf, sizeof(buf)); } }

public:
    explicit StreamBuffer_(gzFile i) : in(i), pos(0), size(0) { assureLookahead(); }

    int  operator *  () const { return (pos >= size) ? EOF : buf[pos]; }
    void operator ++ ()       { pos++; assureLookahead(); }
    int  position    () const { return pos; }
};


static bool eagerMatch_(StreamBuffer_& in, const char* str) {
    for (; *str != '\0'; ++str, ++in)
        if (*str != *in)
            return false;
    return true; }


static void skipWhitespace_(StreamBuffer_& in)
{
    while ((*in >= 9 && *in <= 13) || *in == 32)
        ++in; 
}


static void skipLine_(StreamBuffer_& in)
{
    for (;;)
    {
        if ( *in == EOF ) return;
        if (*in == '\n') 
            { ++in; return; }
        ++in; 
    } 
}


static int parseInt_(StreamBuffer_& in)
{
    int     val = 0;
    bool    neg = false;
    skipWhitespace_(in);
    if      (*in == '-') 
        neg = true, ++in;
    else if (*in == '+') 
        ++in;
    if (*in < '0' || *in > '9')
    {
        std::cerr << "PARSE ERROR! Unexpected char: "<< *in <<std::endl;
        exit(3);
    }
    while (*in >= '0' && *in <= '9')
        val = val*10 + (*in - '0'),
        ++in;
    return neg ? -val : val; 
}


static int parseVar_(StreamBuffer_& in)
{
    int parsed_val = -1, var = -1;
    for (;;)
    {
        skipWhitespace_(in);
        if( *in < '0' || *in > '9' )
        {
            ++in;
            continue;
        }
        parsed_val = parseInt_(in);
        if (parsed_val > 0)
        {
            var = parsed_val ;
            break;
        }
    }
    return var;
}

static void readClause_(StreamBuffer_& in, std::vector<int>& lits)
{
    int     parsed_lit;
    lits.clear();
    for (;;)
    {
        parsed_lit = parseInt_(in);
        if (parsed_lit == 0) 
            break;
        lits.emplace_back( parsed_lit );
    }
}


static void parse_DIMACS_main(const char* filename, std::vector<std::vector<int>> & clauses, int& max_vars)
{
    gzFile in_gz = gzopen(filename, "rb");
    StreamBuffer_ in(in_gz);

    std::vector<int> lits;
    int nb_clauses = 0;
    int cnt     = 0;
    for (;;)
    {
        skipWhitespace_(in);
        if (*in == EOF) 
            break;
        else if (*in == 'p')
        {
            if (eagerMatch_(in, "p cnf"))
            {
                max_vars    = parseInt_(in);
                nb_clauses = parseInt_(in);
            }
            else
            {
                std::cerr <<"PARSE ERROR! Unexpected char: "<<*in << std::endl;
                exit(3);
            }
        } 
        else if (*in == 'c' || *in == 'p')
            skipLine_(in);
        else
        {
            cnt++;
            readClause_(in, lits);
            clauses.emplace_back(lits); 
        }
    }
       gzclose(in_gz);

    if (cnt != nb_clauses)
        std::cerr << "PARSE ERROR! vars:"<<max_vars<<" DIMACS header mismatch: wrong number of clauses" << std::endl;
}



static std::vector<std::vector<int>> parse_NuSMV_DIMACS_main(std::string filename,
                                    std::vector<std::unordered_set<int>>& variable_bloc,
                                    int& max_vars)
{
    std::vector<std::vector<int>> clauses;
    gzFile in_gz = gzopen(filename.c_str(), "rb");
    StreamBuffer_ in(in_gz);

    std::vector<int> lits;
    int nb_clauses = 0;
    int cnt        = 0;
    int time       = -1;
    std::unordered_set<int> set_vars_time;

    // Get information about pure variables
    for(;;)
    {
        if (*in == EOF)
            break;

        if ( *in == 'c' )
        {
            ++in; ++in;
            // c @@@@@@ Time
            if ( *in == '@' )
            {
                if( !set_vars_time.empty() )
                    variable_bloc.emplace_back( set_vars_time );
                set_vars_time.clear();
                time++;
            }
            // c CNF variable
            else if ( *in == 'C' )
            {
                int var_time = parseVar_(in);
                set_vars_time.insert( var_time );
            }
            // c model
            else if ( *in == 'm' )
            {
                skipLine_(in);
                break;
            }
            skipLine_(in);
        }
        else
        {
            ++in;
        }
    }
    if( set_vars_time.empty() )
        throw std::runtime_error("dimacs: parse_NUMSV_DIMACS");
    variable_bloc.emplace_back( set_vars_time );


    // Do normal parsing
    for (;;)
    {
        skipWhitespace_(in);
        if (*in == EOF)
            break;
        else if (*in == 'p')
        {
            if (eagerMatch_(in, "p cnf"))
            {
                max_vars    = parseInt_(in);
                nb_clauses = parseInt_(in);
            }
            else
            {
                std::cerr <<"PARSE ERROR! Unexpected char: "<<*in << std::endl;
                exit(3);
            }
        }
        else if (*in == 'c' || *in == 'p')
            skipLine_(in);
        else
        {
            cnt++;
            readClause_(in, lits);
            clauses.emplace_back(lits); 
        }
    }
    gzclose(in_gz);

    if (cnt != nb_clauses)
        std::cerr << "PARSE ERROR! vars:"<<max_vars<<" DIMACS header mismatch: wrong number of clauses" << std::endl;
    return clauses;
}