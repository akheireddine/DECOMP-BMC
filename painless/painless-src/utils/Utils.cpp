
#include "Utils.h"

#include <stdlib.h>

std::vector<bool> getBooleanVectorVar(std::vector< std::vector< int > >& vec, int len)
{
    std::vector<bool> bool_vec(len,false);

    for(auto& cl : vec)
        for(int v : cl){
            if( abs(v) <= len )
                bool_vec[ abs( v ) ] = true;
        }
    return bool_vec;
}