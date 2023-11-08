// -----------------------------------------------------------------------------
// Copyright (C) 2017  Ludovic LE FRIOUX
//
// This file is part of PaInleSS.
//
// PaInleSS is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See theprintf GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program.  If not, see <http://www.gnu.org/licenses/>.
// -----------------------------------------------------------------------------

#pragma once

#include <iostream>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <string>
#include <string.h>


extern std::map<std::string, std::string> params;

extern char * filename;

/// Class to manage the parameters
class Parameters
{
public:
   /// Init the parameters.
   static void init(int argc, char ** argv)
   {
      for (int i = 1; i < argc; i++) {
         char * arg = argv[i];

         if (arg[0] != '-' && filename == NULL) {
            filename = arg;
            continue;
         }

         char * eq = strchr(arg, '=');

         if (eq == NULL) {
            params[arg+1];
         } else {
            *eq = 0;

            char * left  = arg+1;
            char * right = eq+1;

            params[left] = right;
         }
      }
   }

   /// Get the input cnf filename.
   static char * getFilename()
   {
      static char * ret = filename;

      return ret;
   }

   static std::string getFilenameStr()
   {
      std::string ret = filename;

      return ret;
   }

   static void removeExtensionFromFilename(std::string extension)
   {
      std::string str_filename(filename);
      if( str_filename.find(extension) != std::string::npos )
         str_filename.erase(str_filename.find_last_of("."), std::string::npos);

      filename = strdup(str_filename.c_str());
   }

   /// Print all parameters.
   static void printParams(std::ofstream& opf)
   {
      opf << "c filename " << filename <<std::endl;;

      opf << "c ";

      for (std::map<std::string,std::string>::iterator it = params.begin();
           it != params.end(); it++)
      {
         if (it->second.empty()) {
            opf << it->first << ", ";
         } else {
            opf << it->first << "=" << it->second << ", ";
         }
      }

      opf << std::endl;;
   }

   /// Return true if the parameters is set otherwise false.
   static bool isSet(const std::string & name)
   {
      return params.find(name) != params.end();
   }

   /// Get the string value of a parameters with a default value.
   static const std::string getParam(const std::string & name,
                                const std::string & defaultValue)
   {
      if (isSet(name))
         return params[name];

      return defaultValue;
   }

   /// Get the string value of a parameter.
   static const std::string getParam(const std::string & name)
   {
      return getParam(name, "");
   }

   /// Get the int value of a parameter with a default value.
   static int getIntParam(const std::string & name, int defaultValue)
   {
      if (isSet(name))
         return atoi(params[name].c_str());

      return defaultValue;
   }

   static float getFloatParam(const std::string & name, float defaultValue)
   {
      if (isSet(name))
         return atof(params[name].c_str());
      return defaultValue;
   }

   static void setIntParam(const std::string & name, int new_val)
   {
      params[name] = std::to_string(new_val);
   }
};
