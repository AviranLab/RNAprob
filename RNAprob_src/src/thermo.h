
#if !defined(THERMO_H)
#define THERMO_H


#include "defines.h"


struct thermo //this structure contains the DH and DS values for helical stacks

{
	char DH[maxfil],DS[maxfil],HELIX[maxfil];
 	int dh[5][5][5][5];
   int ds[5][5][5][5];
   int dhi,dsi,dss; //initiation dh and ds and the ds penalty for symmetry
   int dha,dsa; //dh and ds penalty for terminal AU pair
   int read();
   thermo(char *path); //path is a path to the datafiles

};

#endif