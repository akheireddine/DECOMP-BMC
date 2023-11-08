// Copyright (C) 2011 Microsoft Research
// CM Wintersteiger, 2012

#ifndef _DIMACS_PARSER_H_
#define _DIMACS_PARSER_H_

#include <vector>

class DimacsParser
{
public:
  DimacsParser(void);
  ~DimacsParser(void);

  bool readDimacsFile(const char *filename);
  bool readDimacsFile(std::string &filename);
  bool readDimacsFile(const char *filename, long long fraction, long long total);

  virtual bool addClause(const std::vector<signed> &literals) = 0;
  virtual void setVariableMax(unsigned n) = 0;
  virtual void setClauseMax(unsigned n) = 0;

protected:
  bool readHeader(unsigned * vmax = NULL, unsigned * cmax = NULL);  
  bool readClauses(long long limit = -1);
  bool readComments();
  void readWhitespace();
  void readClause(std::vector<signed> & temp);
  int  readInt();
  bool readString(const char *string);

private:
  static const unsigned buf_size = 1048576;
  unsigned char fbuffer[buf_size];
  size_t parser_inx, parser_size;
  FILE * file;
  
  inline unsigned char peek(void) const;
  inline void fwd(void);
  inline void buffer(void);
  inline bool eof(void) const;
};

// --- Inlines

inline unsigned char DimacsParser::peek(void) const {
    return fbuffer[parser_inx];
}

inline void DimacsParser::fwd(void) {
    parser_inx++;
    buffer();
}

inline void DimacsParser::buffer(void) {
    if (parser_inx >= parser_size) {
        parser_inx = 0;
        parser_size = fread(fbuffer, 1, sizeof(fbuffer), file);
    }
}

inline bool DimacsParser::eof(void) const {
    return parser_size == 0;
}

#endif