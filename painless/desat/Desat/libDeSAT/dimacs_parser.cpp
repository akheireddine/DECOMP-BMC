// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2010

#include <iostream>

#include "dimacs_parser.h"

DimacsParser::DimacsParser(void) : parser_inx(0),parser_size(0),file(0) {
}

DimacsParser::~DimacsParser(void) {
}

bool DimacsParser::readDimacsFile(const char *filename)
{
  file = fopen(filename, "r");
  if (file == 0)
      throw std::runtime_error("File cannot be opened.");

  buffer();

  if (!readHeader()) return false;
  if (!readClauses()) return false;

  fclose(file);
  file = 0;
  parser_inx = 0;
  parser_size = 0;

  return true;
}

bool DimacsParser::readDimacsFile(const char *filename, long long fraction, long long total)
{
  file = fopen(filename, "r");
  if (file == 0)
      throw std::runtime_error("File cannot be opened.");

  buffer();

  unsigned vmax, cmax;
  if (!readHeader(&vmax, &cmax)) return false;

  long long fsize = ftello(file);
  long long frac_size = fsize/total;

  if (fraction != 0)
  {
      fseek(file, fraction * frac_size, SEEK_SET);
      while (readInt() != 0 && !eof()) ; // Note: this could leaves us inside a comment
  }

  throw std::runtime_error("NYI"); // parser does not remember file position for data in the buffer.

  if (!readClauses((fraction+1) * frac_size)) return false;

  fclose(file);
  file = 0;
  parser_inx = 0;
  parser_size = 0;

  return true;
}

bool DimacsParser::readDimacsFile(std::string &filename)
{
  return readDimacsFile(filename.c_str());
}

bool DimacsParser::readHeader(unsigned * vmax, unsigned * cmax)
{
  std::string str;
  int i;

  readComments();

  if (!readString("p cnf"))
      throw std::runtime_error("File is not in DIMACS format.");

  i = readInt();
  if (i <= 0)
      throw std::runtime_error("Invalid number of variables in the file.");
  if (vmax) *vmax = i;
  setVariableMax(i);

  i = readInt();
  if (i <= 0)
      throw std::runtime_error("Invalid number of clauses in the file.");
  if (cmax) *cmax = i;
  setClauseMax(i);

  return true;
}

bool DimacsParser::readClauses(long long limit)
{
  std::vector<signed> temp;

  while (!eof()) // check limit.
  {
    readWhitespace();
    readComments();

	if (eof() || peek()=='%')
		return true; // old satlib format?

    readClause(temp);
    /*
    std::cout << "Adding: ";
    for (unsigned i = 0; i < temp.size(); i++)
        std::cout << " " << temp[i];
    std::cout << std::endl;
	*/
    if (!addClause(temp))
		return false;
  }

  return true;
}

bool DimacsParser::readComments()
{
  while (peek()=='c' && !eof())
  {
    do {
      fwd();
    } while (peek()!='\n' && !eof());
    fwd();
  }
  return true;
}

void DimacsParser::readWhitespace()
{
  while (!eof() && (peek()=='\n' || peek()=='\r' || peek()==' '))
    fwd();
}

void DimacsParser::readClause(std::vector<signed> & temp)
{
  int i;
  temp.clear();
  for (i = readInt(); i != 0; i = readInt())
    temp.push_back(i);
}

int DimacsParser::readInt() {
    bool neg = false;
    int val = 0;

    while ((peek() >= 9 && peek() <= 13) || peek() == 32)
        fwd();

    if (peek() == '-') {
        neg = true;
        fwd();
    }

    if (peek() < '0' || peek() > '9')
        throw new std::runtime_error("Unexpected character in input");

    while (peek() >= '0' && peek() <= '9') {
        val = val*10 + (peek() - '0'),
        fwd();
    }

    return neg ? -val : val;
}

bool DimacsParser::readString(const char *string) {
    size_t i = 0;
    while (peek() == string[i]) {
        i++;
        fwd();
        if (string[i] == '\0') return true;
    }
    return false;
}
