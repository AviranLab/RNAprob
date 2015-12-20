/*
 * Dynalign
 *  by David Mathews, copyright 2002, 2003,2004,2005
 *
 * Contributors: Chris Connett and Andrew Yohn, 2006
 */

#include "stackclass.h"

#include "defines.h"

void stackclass::allocate_stack() {
	int i;
  
	stackenergy =new integersize [maximum];
	stack=new int *[maximum];
	for (i=0;i<maximum;i++) stack[i] = new int [4];
}

stackclass::stackclass(int stacksize) {
	maximum = stacksize;
	size = 0;
	allocate_stack();
}

bool stackclass::pull(int *i,int *j, int *open, 
                      integersize *energy, int *pair) {
		
	if (size==0) return false;
	else {
		size--;
		*i = stack[size][0];
		*j = stack[size][1];
		*open = stack[size][2];
		*energy = stackenergy[size];
		*pair = stack[size][3];
		return true;
	}
}
	
void stackclass::push(int i,int j, int open, 
                      integersize energy, int pair){
	int k;

	if (size == maximum) {
		//allocate more space:
		stackclass *temp;
		temp = new stackclass(maximum);
		for (k=0;k<maximum;k++) {
			temp->push(stack[k][0],stack[k][1],stack[k][2],stackenergy[k],stack[k][3]);
		}
		delete_array();
		maximum = 2*maximum;

		allocate_stack();
		for (k=0;k<(maximum/2);k++) {
			temp->pull(&(stack[k][0]),&(stack[k][1]),&(stack[k][2]),&(stackenergy[k]),&(stack[k][3]));
		}

		delete temp;
	}
		
	stack[size][0] = i;
	stack[size][1] = j;
	stack[size][2] = open;
	stackenergy[size] = energy;
	stack[size][3] = pair;
	size++;
}
	
void stackclass::delete_array() {
	for (int i = 0; i < maximum; i++) {
    delete[] stack[i];
  }
	delete[] stack;
	delete[] stackenergy;
}

stackclass::~stackclass() {
	delete_array();
}
