#pragma once

#include <vector>
#include <streambuf>
#include <iostream>
#include <zlib.h>
#include <math.h> 
#include <set>
#include <unordered_set>
#include "EnvBMC.hh"



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


static bool eagerMatchStr(std::string in, const char* str) {
    int i = 0;
    for (; *str != '\0'; ++str, ++i)
        if (*str != in[i])
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

static void readClauseAddPartition_(StreamBuffer_& in, EnvBMC &env)
{
    int     parsed_lit;
    int     id_partition = -1;
    std::vector<int> lits;

    for (;;)
    {
        parsed_lit = parseInt_(in);
        if(abs(parsed_lit) > env.nb_variables)
            env.nb_variables = abs(parsed_lit);
        if(env.log)
            env.logFile<<parsed_lit<<" ";
        if (parsed_lit == 0){
            if(env.log)
                env.logFile<<"\n";
            break;
        }
        while(env.info_variables.size() <= abs(parsed_lit))
            env.info_variables.emplace_back(VB());
        // printf(" %d ",env.info_variables[abs(parsed_lit)].num_partition);
        if(id_partition == -1)
            id_partition = env.info_variables[abs(parsed_lit)].num_partition;
        else if (id_partition != env.info_variables[abs(parsed_lit)].num_partition && 
                 env.info_variables[abs(parsed_lit)].num_partition >= 0)
            id_partition = -2;
        lits.emplace_back( parsed_lit );
    }
    if(id_partition < 0 || env.nb_leafs==0)// goes to reconciliator G
        env.clauses_g.emplace_back(lits);
    else
        env.clauses_partition[id_partition].emplace_back(lits);
    // printf(" >> update partition %d\n",id_partition);
    // printf("G %d   Partitions:\n",env.clauses_g.size());
    // for(int i = 0; i < env.nb_leafs; i++)
        // printf("\t Partition %d : %d\n",i, env.clauses_partition[i].size());
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
        std::cerr << "c PARSE ERROR! vars:"<<max_vars<<" DIMACS header mismatch: wrong number of clauses" << std::endl;
}

static bool operator_ltlspec(char in){
    return (in == ')'  || in == '&' || in == '|'  || 
        in == '>'  || in == ' ' || in == '(' || in == '!' || in == 'G' || in == 'F' || in == 'U' || 
        in == 'X' || in == 'R' || in == 'W' || in == '&' || in == '|' || in == '<' || 
        in == '>' || in == ' ' || in == 'V' || in == '(' || in == '!' || in == '<');
}

static void readLTLSPEC(std::string filename, std::unordered_set<std::string> &ltl_varnames){
    std::string ltl="";
    std::string lastLine="";            
    std::ifstream in(filename.c_str(),std::ifstream::binary);
    if(!in) return;
    in.seekg(0, std::ifstream::end);
    const std::streamoff len = in.tellg();
    //empty file
    if(len == 0)
    {
        printf("c No LTL found, skip\n");
        return;
    }
    int buf_size = 128;
    std::vector<char> buf;
    while(in)
    {   
        if(buf_size > len)
        {
            buf_size = len;
        }
        buf.resize(buf_size);
        in.seekg(0 - buf_size, std::ifstream::end);
        in.read(&buf[0],buf_size);
        //all content is in the buffer or we already have the complete last line
        if(len == buf_size || std::count(buf.begin(), buf.end(), '\n') > 1)
        {
            break;
        }
        //try enlarge the buffer
        buf_size *= 2;
    }
    //find the second line seperator from the end if any
    auto i = std::find(++buf.rbegin(),buf.rend(), '\n');
    lastLine.assign(i == buf.rend() ?  buf.begin() : buf.begin() + std::distance(i, buf.rend()), buf.begin() + buf_size);

    if(!eagerMatchStr(lastLine, "c LTLSPEC "))
    {
        printf("c No LTL found, skip\n");
        return;
    }

    // Temporary string used to split the string.
    bool start = false, end= false;
    for (int i = 10; i < lastLine.length(); i++){
        char c = lastLine[i];
        if(c == '('){
            start = true;
            continue;
        }
        if (c == ')' || c == ' '){
            end = true;
            ltl_varnames.insert(ltl);
            ltl.clear();
        }
        if(!operator_ltlspec(c))
            ltl += c; 
    }
}


std::string parseNameVar(StreamBuffer_& in)
{
    std::string name = "";
    for (;;){
        skipWhitespace_(in);
        if(eagerMatch_(in,"Variable")){
            skipWhitespace_(in);
            while(*in != '\n'){
                name += *in; ++in;}
            break;}
        else{
            ++in; continue;}}
    return name;
}

int parseVarLoopStep(StreamBuffer_& in)
{
    int parsed_val = -1, var = -1;
    for (;;){
        skipWhitespace_(in);
        if( *in < '0' || *in > '9' ){
            ++in; continue;}
        parsed_val = parseInt_(in);
        if (parsed_val >= 0){
            var = parsed_val; break;}}
    return var;
}

static void parse_NuSMV_DIMACS_main(const char* filename, EnvBMC &env)
{
    gzFile in_gz = gzopen(filename, "rb");
    StreamBuffer_ in(in_gz);

    int cnt     = 0;
    int time        = -1;
    int id_group    = 0;
    env.idmax_model = 0;

    readLTLSPEC(filename, env.ltl_varnames);

    // Get information about pure variables
    for(;;)
    {
        skipWhitespace_(in);
        if (*in == EOF) 
            break;
        if ( *in == 'c' )
        {
            ++in; ++in;
            // c Time steps from...
            if ( *in == 'T'){
                if(eagerMatch_(in,"Time steps from 0 to ")){
                    env.k = parseVar_(in); // not adding the initial step k=0
                    if(env.nb_leafs != 0)
                        env.nb_grp_steps = ceil((env.k)/env.nb_leafs);
                    else if (env.nb_grp_steps != 0)
                        env.nb_leafs = ceil((env.k)/env.nb_grp_steps);
                    else{
                        env.nb_leafs = 0;
                        env.nb_grp_steps = 0;
                    }
                    env.clauses_partition.resize(env.nb_leafs);
                }
            }
            // c @@@@@@ Time
            if ( *in == '@' ){
                if(env.nb_grp_steps > 0 && time > 0 && time%env.nb_grp_steps==0 && (id_group+1) < env.nb_leafs)
                    id_group++;
                time++;
            }
            // c CNF variable
            else if ( *in == 'C' )
            {
                int var_time = parseVar_(in);
                std::string name = parseNameVar(in);
                env.idmax_model = var_time;
                while(env.info_variables.size() <= var_time)
                    env.info_variables.emplace_back(VB());
                env.info_variables[var_time].num_step = time;
                if(time == 0)
                    env.info_variables[var_time].num_partition = -1;
                else
                    env.info_variables[var_time].num_partition = id_group;
                // printf("time %d var %d  ==> %d\n",time, var_time,env.info_variables[var_time].num_partition);
                env.info_variables[var_time].in_property = false;
                if(env.ltl_varnames.find(name) != env.ltl_varnames.end())
                    env.info_variables[var_time].in_property = true;
            }
            // c LOOP variable
            else if ( *in == 'L' )
                env.rename_loop.emplace_back(parseVarLoopStep(in));
            // c STEP variable
            else if ( *in == 'S'){
                int v = parseVarLoopStep(in);
                while(env.info_variables.size() <= abs(v))
                    env.info_variables.emplace_back(VB());
                env.info_variables[v].num_step = parseVarLoopStep(in);
                env.rename_step.emplace_back(v);
            }
            // c model
            else if ( *in == 'm' )
                skipLine_(in);
            skipLine_(in);
        }
        else if (*in == 'p')
        {
            if (eagerMatch_(in, "p cnf"))
            {
                env.nb_variables    = parseInt_(in);
                env.nb_clauses      = parseInt_(in);
                if(env.log)
                    env.logFile<< "p cnf "<<env.nb_variables<<" "<<(env.nb_clauses+env.k)<<"\n";
            }
            else
            {
                std::cerr <<"PARSE ERROR! Unexpected char: "<<*in << std::endl;
                exit(3);
            }
        } 
        else
        {
            cnt++;
            readClauseAddPartition_(in, env);
        }
    }
    if(env.log)
        env.logFile.flush();
    gzclose(in_gz);
    env.nb_clauses = cnt;
    if (cnt != env.nb_clauses)
        std::cerr << "c PARSE ERROR! vars:"<<env.nb_variables<<" DIMACS header mismatch: wrong number of clauses " <<cnt<< std::endl;
}

