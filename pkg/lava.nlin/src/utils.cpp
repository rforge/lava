#include "utils.h"

/* get the list element named str, or return NULL */
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = Rf_getAttrib(list, R_NamesSymbol);
  int i;
  
  for (i = 0; i < Rf_length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  return elmt;
}

