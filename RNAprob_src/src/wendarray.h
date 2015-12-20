/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef WENDARRAY_H
#define WENDARRAY_H

#undef BOUNDS

#if defined (BOUNDS)
	#include <iostream>
#endif

#include "defines.h"

//This class encapsulates the w3 and w5 arrays for Dynalign
class wendarray {
private:
	int N1,N2;
	//int *itN2dN1; // array to store i times N2 divided by N1 for
                  // values of i
	int *Lowlimit;

#if defined (BOUNDS)
	//Bounds checking is on.
	int *Highlimit;
#endif

	public:
		integersize **array;
		wendarray();
		wendarray(int n1, int n2, int *lowlimit, int *highlimit);
		void allocate(int n1, int n2, int *lowlimit, int *highlimit);
		~wendarray();
		integersize &f(int i, int j);
};

inline integersize &wendarray::f(int i, int j) {

  return array[i][j/*-itN2dN1[i] + M*/];
}

#endif
