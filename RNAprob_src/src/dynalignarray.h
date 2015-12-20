/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005,2006
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNARRAY_H
#define DYNALIGNARRAY_H

#include "defines.h"

// This class encapsulates the large 4-dimensional arrays used
// by Dynalign

//Note that highlimit[i] and lowlimit[i] are the spans of allowed alignments in seq2 for nucletide 
//i from sequence 1.  IMPORTANT:  These arrays must persist until after the call of ~dynalignarray.
class dynalignarray {
	private:
  int *Lowlimit,*Highlimit;          
  int N,N2,Ndiff; // length of sequence one and two
  bool optimalonly;
  int infinite;

	public:
		integersize ****array;
		dynalignarray(int n1, int n2, int *lowlimit, int *highlimit, bool Optimalonly=false);
		dynalignarray();
		~dynalignarray();
		void allocate(int n1, int n2, int *lowlimit, int *highlimit, bool Optimalonly=false);
		integersize &f(int i, int j,int k, int l);
};

inline integersize &dynalignarray::f(int i, int j,int k, int l) {
  if (i>N&&j>N) {
    i -= N;
    j -= N;
    k -= N2;
    l -= N2;
  }
  return array[i][j][k][l];
}

#endif
