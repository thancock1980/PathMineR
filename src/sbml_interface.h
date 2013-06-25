#ifndef __sbml_interface__h_
#define __sbml_interface__h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>   
using namespace std; 

#include <sbml/SBMLTypes.h>

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

extern "C" SEXP readsbmlfile(SEXP FILENAME);

SEXP getSpeciesFrame(Model* model);
SEXP getReactionList(Model *model);

#endif
