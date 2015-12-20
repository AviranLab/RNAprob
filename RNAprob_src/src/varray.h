/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef VARRAY_H
#define VARRAY_H

#include "defines.h"


//Varray encapsulates the v array. i.e. the lowest free energy for fragments from i to j in seq1,
//k to l in seq 2, with i paired to j, k paired to l, i aligned to k, and j aligned to l.

//Note that highlimit[i] and lowlimit[i] are the spans of allowed alignments in seq2 for nucletide 
//i from sequence 1.  IMPORTANT:  These arrays must persist until after the call of ~varray.

class varray {
	private:
		//int M;//maxseparation parameter
		int N,N2,Ndiff;//length of sequence one and two
		bool optimalonly;
		integersize infinite;
		bool **tem;
		int *Lowlimit, *Highlimit;

	public:
		integersize ****array;

		//constructor: n1 is the length of sequence 1, n2 is the length of sequence 2, 
		varray(int n1, int n2, int *lowlimit, int *highlimit, bool **Tem, bool Optimalonly=false);
		varray();
		~varray();
		void allocate(int n1, int n2, int *lowlimit, int *highlimit, bool **Tem, bool Optimalonly=false);
		integersize &f(int i, int j,int k, int l);
};

inline integersize &varray::f(int i, int j,int k, int l) {
  if (i > N && j > N) {
    i -= N;
    j -= N;
    k -= N2;
    l -= N2;
  }
  
  if (j > N) {
    if (tem[i][j-N]) {
      return array[i][j][k][l];
    } else {
      return infinite;
    }
  } else if (tem[j][i]) {
    return array[i][j][k][l];
  } else {
    return infinite;
  }
}

#endif
