/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#ifndef STACKCLASS_H
#define STACKCLASS_H

#include "defines.h"

class stackclass {
private:
	void allocate_stack();

public:
	int size,**stack,maximum;
	integersize *stackenergy;

	stackclass(int stacksize = 50);
	~stackclass();

	bool pull(int *i,int *j, int *open, integersize *energy, int *pair);
	void push(int i,int j, int open, integersize energy, int pair);
	
	void delete_array();
};

#endif
