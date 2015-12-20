#if !defined(ALLTRACE_H)
#define ALLTRACE_H

#ifdef _WINDOWS
#include "../RNAstructure_windows_interface/TProgressDialog.h"
#else

#ifdef _JAVA_GUI
#include "../RNAstructure_java_interface/SWIG/TProgressDialog.h"
#else
#include "TProgressDialog.h"
#endif // JAVA GUI

#endif //WINDOWS

#include "defines.h"
#include "structure.h"
#include "algorithm.h"
#include "rna_library.h"

////////////////////////////////////////////////////////////////////////
//arrayclass encapsulates the large 2-d arrays of w and v, used by the dynamic
//	algorithm
class atarrayclass {
   int Size;

   public:
   	
      int k;
      integersize **dg;
      integersize infinite;
	  bool allocated;
      

      //the constructor allocates the space needed by the arrays
   		atarrayclass(int size);
		atarrayclass();
		void allocate(int size);

      //the destructor deallocates the space used
      ~atarrayclass();

      //f is an integer function that references the correct element of the array
   	inline integersize &f(int i, int j);
};

//f is an integer function that references the correct element of the array
inline integersize &atarrayclass::f(int i, int j) {
      	
   if (i>j) {
        return infinite;
    }
   else if (i>Size) return f(i-Size,j-Size);//dg[i-Size][j-Size];
   else return dg[i][j-i];
         
}

//Function to generate all suboptimal structures for sequence in structure ct and store them in structure ct.
	//data is a pointer to datatable, where the thermodynamic parameters are stored.
	//percentdeta is the maximum percent difference from the lowest free energy structure.
	//absolutedelta is the maximum folding free energy change (in kcal/mol * conversionfactor, as defined in defines.h)
	//update is a TProgressDialog, used to track progress.
	//save is the name of a savefile, which generates save files that can be used by realltrace
	//NoMBLoop = whether multibranch loops are allowed, wehere true indicates NO multibranch loops
void alltrace(structure* ct,datatable* data, int percentdelta, int absolutedelta, TProgressDialog* update, char* save, bool NoMBLoop=false);
void readalltrace(char *filename, structure *ct, 
			 int *w5,  
			 atarrayclass *v, atarrayclass *w, atarrayclass *wmb, atarrayclass *wmbl, atarrayclass *wl, atarrayclass *wcoax,
			 atarrayclass *w2, atarrayclass *wmb2, forceclass *fce, bool *lfce, bool *mod, datatable *data);

void realltrace(char *savefilename, structure *ct, int percentdelta, int absolutedelta, char *ctname = NULL);


#define startingsize 500  //maximum number of structure fragments to start in alltracestructurestack (below)
//#define startingrefinementstacksize 25


//a stack to keep track of partially refined structures
	//contains a stackclass to keep track of where refinements need to occur 
class alltracestructurestack {

	public:
		int **basepairs;
		int maximumsize;
		int current; //current location in stack
		alltracestructurestack(int size, int sizeofstack=startingsize);
		stackclass *refinementstack;//a stack to keep track of where refinement is occuring
		~alltracestructurestack();
		void addpair(int i, int j, int index);
		void push();
		void push(integersize totalenergy, bool topair,int pairi, int pairj, bool tostack, int i, int j, int open, integersize energy, int pair, bool topair2,int pairi2, int pairj2, bool tostack2, int i2, int j2, int open2, integersize energy2, int pair2,bool tostack3=false, int i3=0, int j3=0, int open3=0, integersize energy3=0, int pair3=0);
		int numberofnucs; //length of the sequence
		void allocatearrays();
		void deletearrays();
		void pull();
		void pushtorefinement(int a, int b, int c, integersize d, int e);
		bool pullfromrefinement(int *a, int *b, int *ct, integersize *d, int *e);
		int *energy;
		integersize peekatenergy();
		void placeenergy(int energy);
		
		void flushbullpen();

		bool refined;
		bool bullpentopair,bullpentopair2;
		bool bullpentostack,bullpentostack2,bullpentostack3;
		int bullpeni,bullpenj,bullpenopen,bullpenpair,bullpenpairi,bullpenpairj;
		integersize bullpenenergy,bullpentotalenergy;
		int bullpeni2,bullpenj2,bullpenopen2,bullpenpair2,bullpenpairi2,bullpenpairj2;
		integersize bullpenenergy2;
		int bullpeni3,bullpenj3,bullpenopen3,bullpenpair3;
		integersize bullpenenergy3;
		int readpair(int i);
		void stackup(int index);

		//The following is infrastructure to keep track of stacked nucleotides:
		int stack1[2],stack2[2];
		void nstack(int i, int j, int k=0, int l=0);
		int **stacks;
		int readstacking(int i);

		#if defined (pfdebugmode) 
		ofstream out;

		#endif

};

#endif //!defined ALLTRACE_H
