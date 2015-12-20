/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef DYNALIGNSTACKCLASS_H
#define DYNALIGNSTACKCLASS_H

#include "defines.h"

class dynalignstackclass {
	int **stack;
	int size, max;
	integersize *stackenergy;
	bool *openness;
        void allocate_stack();
#ifdef DYNALIGN_II
        bool *vmodness;
#else
#endif
public:
	dynalignstackclass(int stacksize = 50);
#ifdef DYNALIGN_II
    	bool pull(int *i,int *j, int *a, int *b, 
		  integersize *energy, bool *open,bool *if_vmod);
	void push(int i,int j, int a, 
		  int b, integersize energy,bool open=false, bool if_vmod=false);
#else
	bool pull(int *i,int *j, int *a, int *b, 
            integersize *energy, bool *open);
	void push(int i,int j, int a, 
            int b, integersize energy, bool open=false);
#endif	
	void delete_array();
	~dynalignstackclass();
};

#endif
