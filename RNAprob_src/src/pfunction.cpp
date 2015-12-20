/*Code for partition function calculation of base pair probabilities.
Copyright 2003,2004,2005,2006 by David H. Mathews

Contributions by Zhi John Lu, 2006.

*/


//NOTE: Two compiler flags exist to modify the behaviour of MB Loops:
//SIMPLEMBLOOP uses a model in which every helix end has a 5' and 3' dangle (Except at sequence ends)
//disablecoax is a flag that disables coaxial stacking, but leaves the rest of the model intact



#include "pfunction.h"
#include "boltzmann.h" //for boltzman
#include <math.h>
#include <cstdlib>

#ifdef SMP
	#include <omp.h>
#endif

using namespace std;

#undef pfdebugmode  //flag to indicate debugging
//#define pfdebugmode  //flag to indicate debugging

#undef equiout
//#define equiout //flag to indicate equilibria should be written to file

#undef timer

//#define timer //flag to indicate the code execution should be timed

#undef disablecoax
//#define disablecoax



//These have been moved to rna_library
//inline void write(ofstream *out,double *i) {

//	out->write((char *) i,sizeof(*i));
//}

//inline void read(ifstream *out,double *i) {

//	out->read((char *) i,sizeof(*i));
//}

#define maxinter 30   //maximum length of the internal loops
#define maxasym 30  //maximum asymetry in the internal loops


void calculatepfunction(structure* ct,pfdatatable* data, TProgressDialog* update, char* save, bool quickQ, PFPRECISION *Q,
	pfunctionclass *w, pfunctionclass *v, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl,
	pfunctionclass *wcoax, forceclass *fce,PFPRECISION *w5,PFPRECISION *w3,bool *mod, bool *lfce) {


int ip,jp,ii,jj,jpf,jf,bl,ll,dp;
register int i,j,h,d;
int k,p;
register int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
//int before,after;
//register PFPRECISION e;
register int number,maxj;
register PFPRECISION twoscaling,rarray;
PFPRECISION **curE,**prevE;
PFPRECISION **tempE,**wca;
//bool calculatev;


#ifdef timer
	#include <time.h>
	ofstream timeout;
	int seconds;
	char timerstring[100];
	char timelength[10];
	strcpy(timerstring,"time_pf_");
	sprintf(timelength,"%i",ct->GetSequenceLength());
	strcat(timerstring,timelength);
	strcat(timerstring,".out");

	timeout.open(timerstring);
	timeout<<time(NULL)<<"\n";
	seconds = time(NULL);
#endif

number = (ct->GetSequenceLength());//place the number of bases in a registered integer

//#ifndef SMP //These (for internal loop calculations) are only used in serial code

#ifndef SMP //These (for internal loop calculations) are only used in serial code

curE= new  PFPRECISION *[number+1];
prevE= new PFPRECISION *[number+1];
#else
curE = NULL;
prevE = NULL;
#endif
wca = new PFPRECISION *[number+1];




for (i=0;i<=number;i++) {
	#ifndef SMP //These (for internal loop calculations) are only used in serial code
		curE[i]= new  PFPRECISION [number+1];
   		prevE[i]= new PFPRECISION [number+1];
	#endif
	wca[i]=new PFPRECISION [number+1];
	for(j=0;j<=number;j++) {

		wca[i][j] = (PFPRECISION) 0;
		#ifndef SMP //These (for internal loop calculations) are only used in serial code
			curE[i][j]= (PFPRECISION) 0;
			prevE[i][j]= (PFPRECISION) 0;
		#endif
	}

}

w5[0] = (PFPRECISION) 1;//initialize the random coil contribution to the partition function
w3[number+1] = (PFPRECISION) 1;

force(ct,fce,lfce);

twoscaling = data->scaling*data->scaling;


//This is the fill routine:

if (quickQ) maxj = number;
else maxj = 2*number-1;


for (h=0;h<=( quickQ?(maxj-1):(maxj-1-minloop) );h++){

	d=(h<=(number-1))?h:(h-number+1);

	#ifndef SMP //These (for internal loop calculations) are only used in serial code
	if (h==number) {
		for(i=0;i<=number;i++) {
			for(j=0;j<=number;j++) {
				curE[i][j]= (PFPRECISION) 0;
				prevE[i][j]= (PFPRECISION) 0;
			}
		}
	}
//<<<<<<< pfunction.cpp
//	#endif

	
//=======
	#endif
//>>>>>>> 1.30

	if (((h%10)==0)&&update) update->update((100*h)/(2*ct->GetSequenceLength()));

//<<<<<<< pfunction.cpp
//	int start;
//	int end;
//	if (h<=(number-1)) {
//		start = 1;
//		end = number-h;
//	}
//	else {
//		start = 2*number-h;
//		end = number;
//	}
	
	
	//for (int locali=((h<=(number-1))?1:(2*number-h));locali<=((h<=(number-1))?(number-h):number);locali++){
//	#ifdef SMP
//		#pragma omp parallel for
//	#endif
//	for (int locali=start;locali<=end;locali++){
//		int localj=locali+d;
		
//=======

	int start;
	int end;
	if (h<=(number-1)) {
		start = 1;
		end = number-h;
	}
	else {
		start = 2*number-h;
		end = number;
	}
	
	
	//for (int locali=((h<=(number-1))?1:(2*number-h));locali<=((h<=(number-1))?(number-h):number);locali++){
	#ifdef SMP
		#pragma omp parallel for
	#endif
	for (int locali=start;locali<=end;locali++){
		int localj=locali+d;
		
//>>>>>>> 1.30
		//Test to make sure the fragment is large enough to form pairs:
		if (!(localj<=(number)&&((localj-locali)<=minloop))) {



			//First, calculate V(locali,localj): (The partition function from nucleotides locali to localj, given that locali and localj pair)



			//Start the value of V:
			double localrarray=0.0;
			double locale;


			//Now, test some conditions as to whether V should be evaluated:
//<<<<<<< pfunction.cpp
//			bool calculatev=true;
//			int localii,localjj,localp,localafter;
//			int localip, localjp;
//=======
			bool calculatev=true;
			int localii,localjj,localp,localafter;
//>>>>>>> 1.30

			if (ct->templated) {
				//Test whether these nucleotides are allowed to pair because there is a pairing mask in ct ("tem").
//<<<<<<< pfunction.cpp

//				if (locali>ct->GetSequenceLength()) localii = locali - ct->GetSequenceLength();
//				else localii = locali;
//				if (localj>ct->GetSequenceLength()) localjj = localj - ct->GetSequenceLength();
//				else localjj = localj;
//				if (localjj<localii) {
//					localp = localjj;
//      				localjj = localii;
//					localii = localp;
//				}
//   				if (!ct->tem[localjj][localii]) calculatev = false;

//=======
				if (locali>ct->GetSequenceLength()) localii = locali - ct->GetSequenceLength();
				else localii = locali;
				if (localj>ct->GetSequenceLength()) localjj = localj - ct->GetSequenceLength();
				else localjj = localj;
				if (localjj<localii) {
					localp = localjj;
      				localjj = localii;
					localii = localp;
				}
   				if (!ct->tem[localjj][localii]) calculatev = false;
//>>>>>>> 1.30
			}


			if (fce->f(locali,localj)&SINGLE) {
				//locali or localj is forced single-stranded
				//v->f(locali,localj) = 0;
				calculatev = false;
			}
			if (fce->f(locali,localj)&NOPAIR) {
				//locali or localj is forced into a pair elsewhere
   				//v->f(locali,localj)= 0;
				calculatev = false;
			}



			if (inc[ct->numseq[locali]][ct->numseq[localj]]==0) {
				//These are two nucleotides that cannot form a canonical pair
				//v->f(locali,localj)= 0;
				calculatev = false;
			}


			//force u's into gu pairs, if the user has specified these restraints
//<<<<<<< pfunction.cpp

			for (int localip=0;localip<ct->GetNumberofGU();localip++) {
				if (ct->GetGUpair(localip)==locali) {
					if (ct->numseq[localj]!=3) {

						
//=======
//			for (int localip=0;localip<ct->ngu;localip++) {
//				if (ct->gu[localip]==locali) {
//					if (ct->numseq[localj]!=3) {
//>>>>>>> 1.30

         				calculatev = false;
					}
				}

//<<<<<<< pfunction.cpp

				else if (ct->GetGUpair(localip)==localj) {
       				if (ct->numseq[locali]!=3) {

//=======
//				else if (ct->gu[localip]==localj) {
//       				if (ct->numseq[locali]!=3) {
//>>>>>>> 1.30

         				calculatev = false;
					}
				}
//<<<<<<< pfunction.cpp

				
				else if ((ct->GetGUpair(localip)+number)==localj) {
       				if (ct->numseq[locali]!=3) {

//=======
//				else if ((ct->gu[localip]+number)==localj) {
//       				if (ct->numseq[locali]!=3) {
//>>>>>>> 1.30

						calculatev = false;
					}
				}

			}



			//now check to make sure that this isn't an isolated pair:
			//	(consider a pair separated by a bulge as not! stacked)

			//before = 0 if a stacked pair cannot form 5' to locali
			int localbefore =0;
			if ((locali>1&&localj<(2*number)&&localj!=number)) {
				if ((localj>number&&((locali-localj+number)>minloop+2))||localj<number) {
					localbefore = inc[ct->numseq[locali-1]][ct->numseq[localj+1]];
				}
			}


			//after = 0 if a stacked pair cannot form 3' to locali
			if ((((localj-locali)>minloop+2)&&(localj<=number)||(localj>number+1))&&(locali!=number)) {
				localafter = inc[ct->numseq[locali+1]][ct->numseq[localj-1]];

			}
			else localafter = 0;

			//if there are no stackable pairs to locali.localj then don't allow a pair locali,localj
			if ((localbefore==0)&&(localafter==0)) {
				//v->f(locali,localj)= 0;
				calculatev = false;
			}


			//A large code block for filling V, if all the above conditions pass:
			if (calculatev) {


				//Test to make sure this isn't the end of the sequence
				if (!(locali==(number)||localj==((number)+1))) {





   					//Perhaps locali and localj close a hairpin:
					localrarray=erg3(locali,localj,ct,data,fce->f(locali,localj));

					if ((localj-locali-1)>=(minloop+2)||localj>(number)) {
      					//Perhaps locali,localj stacks over locali+1,localj-1
						if (!mod[locali]&&!mod[localj])  //make sure this is not a site of chemical modification
							localrarray+=erg1(locali,localj,locali+1,localj-1,ct,data)*v->f(locali+1,localj-1);
						else {
							//allow G-U to be modified or a pair next to a G-U to be modified
							if ((ct->numseq[locali]==3&&ct->numseq[localj]==4)||(ct->numseq[locali]==4&&ct->numseq[localj]==3)) {
								localrarray+=erg1(locali,localj,locali+1,localj-1,ct,data)*v->f(locali+1,localj-1);

							}
							else if ((ct->numseq[locali+1]==3&&ct->numseq[localj-1]==4)||(ct->numseq[locali+1]==4&&ct->numseq[localj-1]==3)) {

								localrarray+=erg1(locali,localj,locali+1,localj-1,ct,data)*v->f(locali+1,localj-1);

							}
							else if (locali-1>0&&localj+1<2*number) {
								if ((ct->numseq[locali-1]==3&&ct->numseq[localj+1]==4)||(ct->numseq[locali-1]==4&&ct->numseq[localj+1]==3)) {

									localrarray+=erg1(locali,localj,locali+1,localj-1,ct,data)*v->f(locali+1,localj-1);

								}

							}

						}
					}



//<<<<<<< pfunction.cpp
					//Perhaps locali,localj closes an interior or bulge loop, search all possibilities

//					#ifndef SMP //USe the O(N^3) code, below, only for serial calculations:
					
						/*fill the interior loops' energy localrarray first
						calculate the small loop (size<=5) first
						the larger loop is prefilled with curE[dp][locali] (after sub2)

						d= localj-locali, dp= jp-ip (interior loop)
						locali<ip<number<jp<localj or locali<ip<jp<localj<number
						*/
//						if ((d-1)>=(minloop+3)||localj>number)
//						for (dp=d-3;dp>=((localj>number)?1:minloop+1);dp--) {
//							ll=d-dp-2;

							//calculate every ip,jp when ll <=5: 0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
//							if(ll>=1&&ll<=5) {

//								for (ip=locali+1;ip<=localj-1-dp;ip++){
//	  								jp=ip+dp;
//									if (inc[ct->numseq[ip]][ct->numseq[jp]])
//	  								if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) ) {
   										//using jpf and jf and  instead of jp and localj when localj,jp>number
//										jpf=( (jp<=number)?jp:jp-number);
//										jf=((localj<=number)?localj:localj-number);
//										localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),
//              									fce->f(jpf,jf)) * v->f(ip,jp);
										//locali or localj is modified
//										if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
//		  									localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
//                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
//=======
					//Perhaps locali,localj closes an interior or bulge loop, search all possibilities

					#ifndef SMP //USe the O(N^3) code, below, only for serial calculations:
					
						/*fill the interior loops' energy localrarray first
						calculate the small loop (size<=5) first
						the larger loop is prefilled with curE[dp][locali] (after sub2)

						d= localj-locali, dp= jp-ip (interior loop)
						locali<ip<number<jp<localj or locali<ip<jp<localj<number
						*/
						if ((d-1)>=(minloop+3)||localj>number)
						for (dp=d-3;dp>=((localj>number)?1:minloop+1);dp--) {
							ll=d-dp-2;

							//calculate every ip,jp when ll <=5: 0x1,0x2,0x3,0x4,0x5,1x1,1x2,1x3,1x4,2x2,2x3
							if(ll>=1&&ll<=5) {

								for (ip=locali+1;ip<=localj-1-dp;ip++){
	  								jp=ip+dp;
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
	  								if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) ) {
   										//using jpf and jf and  instead of jp and localj when localj,jp>number
										jpf=( (jp<=number)?jp:jp-number);
										jf=((localj<=number)?localj:localj-number);
										localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),
              									fce->f(jpf,jf)) * v->f(ip,jp);
										//locali or localj is modified
										if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
		  									localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

      								}

								}
							}
							//when size >=6 and <=30;

							else if (ll>=6&&ll<=maxinter)
							{

								//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is stored in rna_library_inter.cpp
								localrarray+=(curE[dp][locali]) *erg2ex(locali,localj,ll,ct,data);

							//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
								for (bl=0;bl<=1;bl++) {
									ip=locali+1+bl;
									jp=ip+dp;
									jpf=( (jp<=number)?jp:jp-number);
									jf=((localj<=number)?localj:localj-number);
									if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) )
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//      								}
//=======
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
									if (abs(ip-locali+jp-localj)<=maxasym) {
		  								localrarray +=(erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf))) * (v->f(ip,jp));
										//locali or localj is modified
		   								if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
											localrarray +=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//								}
//							}
							//when size >=6 and <=30;

//							else if (ll>=6&&ll<=maxinter)
//							{
//=======
									}
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
								//interior energy prefilled (after sub2:) + exterior engergy (the function for erg2in and erg2ex is stored in rna_library_inter.cpp
//								localrarray+=(curE[dp][locali]) *erg2ex(locali,localj,ll,ct,data);

							//considering loop 1x(n-1) (bl=1)and 0 x n(bulge)(bl=0) as stacking bonus on 1x(n-1) and bulge is not allowed
//								for (bl=0;bl<=1;bl++) {
//									ip=locali+1+bl;
//									jp=ip+dp;
//=======
									jp=localj-1-bl;
									ip=jp-dp;
//>>>>>>> 1.30
									jpf=( (jp<=number)?jp:jp-number);
//<<<<<<< pfunction.cpp
//									jf=((localj<=number)?localj:localj-number);
//									if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) )

//									if (inc[ct->numseq[ip]][ct->numseq[jp]])
//									if (abs(ip-locali+jp-localj)<=maxasym) {
//		  								localrarray +=(erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf))) * (v->f(ip,jp));
//										//locali or localj is modified
//		   								if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
//											localrarray +=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
//                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
//=======
									jf=((localj<=number)?localj:localj-number);
									if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) )
									if (inc[ct->numseq[ip]][ct->numseq[jp]])
									if (abs(ip-locali+jp-localj)<=maxasym) {
										localrarray += erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) * v->f(ip,jp);
			   							if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
											localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//									}

//									jp=localj-1-bl;
//									ip=jp-dp;
//									jpf=( (jp<=number)?jp:jp-number);
//									jf=((localj<=number)?localj:localj-number);
//									if ( (ip<=number&&jp>number) || (ip<=number&&localj<=number) )
//									if (inc[ct->numseq[ip]][ct->numseq[jp]])
//									if (abs(ip-locali+jp-localj)<=maxasym) {
//										localrarray += erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) * v->f(ip,jp);
//			   							if ((mod[ip]||mod[jp])&&inc[ct->numseq[ip]][ct->numseq[jp]]&&!(fce->f(ip,jp)&SINGLE))
//											localrarray+=erg2(locali,localj,ip,jp,ct,data,fce->f(locali,ip),fce->f(jpf,jf)) *
//                  								v->f(ip+1,jp-1) * erg1(ip,jp,ip+1,jp-1,ct,data);

//									}
//								}
//=======
									}
								}
//>>>>>>> 1.30

							}
						}
//<<<<<<< pfunction.cpp
//					#else
						//If using SMP, resort to the O(N^4) recursions for internal loops: 
								
//						if (((localj-locali-1)>=(minloop+3))||(localj>(number))) {
							

//							int maxsize;
//							if (localj<=number) maxsize = min(maxinter,localj-locali-minloop-3);//interior fragment
//							else maxsize = min(localj - locali - 3,maxinter); //exterior fragment
//							for (int size=1;size<=maxsize;++size) {
								
//								for (localip=locali+1,localjp=localj-size-1;localip<=locali+size+1&&localip<=number;++localip,++localjp) {
									
									
									/*if (localip>(number)) {
											
										localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip-number),
																fce->f(localjp-number,localj-number))*
																v->f(localip-(number),localjp-(number));

										if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
											//localip or localjp is modified

											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip-number),
																	fce->f(localjp-number,localj-number))*
														v->f(localip-(number)+1,localjp-(number)-1)*
														erg1(localip-number,localjp-number,localip+1-number,localjp-1-number,ct,data);

										}
										//}
									}*/
									//else {
//										if (localjp<=number) {




//											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//																	fce->f(localjp,localj))*
//																v->f(localip,localjp);

//											if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
//												//locali or localj is modified
//												localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//																		fce->f(localjp,localj))*
//																	v->f(localip+1,localjp-1)*erg1(localip,localjp,localip+1,localjp-1,ct,data);
//
//											}
									


									
//										}
//										else {

									 
//											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number))*
//																v->f(localip,localjp);


//											if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
												//locali or localj is modified
//												localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number))*
//																	v->f(localip+1,localjp-1)*erg1(localip,localjp,localip+1,localjp-1,ct,data);

//											}

//										}
//=======
					#else
						//If using SMP, resort to the O(N^4) recursions for internal loops: 
								
						if (((localj-locali-1)>=(minloop+3))||(localj>(number))) {
							

							int maxsize;
							if (localj<=number) maxsize = min(maxinter,localj-locali-minloop-3);//interior fragment
							else maxsize = min(localj - locali - 3,maxinter); //exterior fragment
							for (int size=1;size<=maxsize;++size) {
								int localip, localjp;
								for (localip=locali+1,localjp=localj-size-1;localip<=locali+size+1&&localip<=number;++localip,++localjp) {
									
									
									/*if (localip>(number)) {
											
										localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip-number),
																fce->f(localjp-number,localj-number))*
																v->f(localip-(number),localjp-(number));

										if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
											//localip or localjp is modified

											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip-number),
																	fce->f(localjp-number,localj-number))*
														v->f(localip-(number)+1,localjp-(number)-1)*
														erg1(localip-number,localjp-number,localip+1-number,localjp-1-number,ct,data);

										}
										//}
									}*/
									//else {
										if (localjp<=number) {




											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
																	fce->f(localjp,localj))*
																v->f(localip,localjp);

											if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
												//locali or localj is modified
												localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
																		fce->f(localjp,localj))*
																	v->f(localip+1,localjp-1)*erg1(localip,localjp,localip+1,localjp-1,ct,data);

											}
									


									
										}
										else {

									 
											localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number))*
																v->f(localip,localjp);


											if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]) {
												//locali or localj is modified
												localrarray+=erg2(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),fce->f(localjp-number,localj-number))*
																	v->f(localip+1,localjp-1)*erg1(localip,localjp,localip+1,localjp-1,ct,data);

											}

										}
//>>>>>>> 1.30




									//}
									


						   
								}//for size
							}//for localip
						}//if ((localj-locali-1)>=(minloop+3))||(localj>(number))
						

					#endif //SMP not defined



				}//end of condition to make sure this wasn't the end of the sequence

				//sub1:


				//Perhaps locali,localj closes a multibranch or exterior loop, enumerate all possibilities

				if (((localj-locali-1)>=(2*minloop+4))||(localj>(number))) {


					//consider the exterior loop closed by locali,localj
					if (localj>number) {
						#ifdef SIMPLEMBLOOP
						//5' and 3' dangling ends
						if (locali!=number&&localj!=number+1) localrarray+= erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1])*
							penalty(locali,localj,ct,data)*w3[locali+1]*w5[localj-number-1]*twoscaling;
						//3' dangling end
						else if (locali!=number) localrarray+= erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*
							penalty(locali,localj,ct,data)*w3[locali+1]*w5[localj-number-1]*twoscaling;
						//5' dangling ends
						else if (localj!=number+1) localrarray+= erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*
							penalty(locali,localj,ct,data)*w3[locali+1]*w5[localj-number-1]*twoscaling;
						//no dangling ends
						else localrarray+= w3[locali+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*twoscaling;

						#else //not simplembloop
         				localrarray+= w3[locali+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*twoscaling;


						if (locali!=number) localrarray+= erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*penalty(locali,localj,ct,data)*w3[locali+2]*w5[localj-number-1]*twoscaling;
						if (localj!=(number+1)) localrarray+= erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1])*penalty(locali,localj,ct,data)*w3[locali+1]*w5[localj-number-2]*twoscaling;
						if ((locali!=number)&&(localj!=(number+1))) {
            				localrarray+= data->tstack[ct->numseq[locali]][ct->numseq[localj]][ct->numseq[locali+1]][ct->numseq[localj-1]]*pfchecknp(lfce[locali+1],lfce[localj-1])*w3[locali+2]*
								w5[localj-number-2]*penalty(locali,localj,ct,data)*twoscaling;

						}


						//consider the coaxial stacking of a helix from locali to localip onto helix localip+1 or localip+2 to localj:
						#ifndef disablecoax //a flag that can turn of coaxial stacking
						//first consider a helix stacking from the 5' sequence fragment:
						for (int localip=localj-number-minloop-1;localip>0;localip--) {
							//first consider flush stacking
							localrarray+=
								w3[locali+1]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localj-number-1,localip,ct,data)*
								ergcoaxflushbases(localip,localj-number-1,localj-number,locali,ct,data)*v->f(localip,localj-number-1)*twoscaling;

							if ((mod[localip]||mod[localj-number-1])) if (localj-number-2>0&&notgu(localip,localj-number-1,ct)&&!(fce->f(localip,localj-number-1)&SINGLE)) {
								if (inc[ct->numseq[localip+1]][ct->numseq[localj-number-2]]) {
									localrarray+=
										w3[locali+1]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localj-number-1,localip,ct,data)*
										ergcoaxflushbases(localip,localj-number-1,localj-number,locali,ct,data)*v->f(localip+1,localj-number-2)*
										erg1(localip,localj-number-1,localip+1,localj-number-2,ct,data)*twoscaling;
								}

							}


							if (localj-number-2>0) {
								//now consider an intervening nuc
								if(locali<number) {
									localrarray+=
										w3[locali+2]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localip,localj-number-2,ct,data)*
										ergcoaxinterbases2(localip,localj-number-2,localj-number,locali,ct,data)*v->f(localip,localj-number-2)*twoscaling
										*pfchecknp(lfce[localj-number-1],lfce[locali+1]);


									if ((mod[localip]||mod[localj-number-2])) if (inc[ct->numseq[localip+1]][ct->numseq[localj-number-3]]&&notgu(localip,localj-number-2,ct)
										&&!(fce->f(localip,localj-number-2)&SINGLE)) {
										localrarray+=
											w3[locali+2]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localip,localj-number-2,ct,data)*
											ergcoaxinterbases2(localip,localj-number-2,localj-number,locali,ct,data)*v->f(localip+1,localj-number-3)*
											erg1(localip,localj-number-2,localip+1,localj-number-3,ct,data)*twoscaling*pfchecknp(lfce[localj-number-1],lfce[locali+1]);


									}
								}


								//consider the other possibility for an intervening nuc
								localrarray+=
									w3[locali+1]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localip+1,localj-number-2,ct,data)*
									ergcoaxinterbases1(localip+1,localj-number-2,localj-number,locali,ct,data)*v->f(localip+1,localj-number-2)*twoscaling
									*pfchecknp(lfce[localj-number-1],lfce[localip]);


								if ((mod[localip+1]||mod[localj-number-2])) if (inc[ct->numseq[localip+2]][ct->numseq[localj-number-3]]&&notgu(localip+1,localj-number-2,ct)
									&&!(fce->f(localip+1,localj-number-2)&SINGLE)) {
									localrarray+=
										w3[locali+1]*w5[localip-1]*penalty(locali,localj,ct,data)*penalty(localip+1,localj-number-2,ct,data)*
										ergcoaxinterbases1(localip+1,localj-number-2,localj-number,locali,ct,data)*v->f(localip+2,localj-number-3)
										*erg1(localip+1,localj-number-2,localip+2,localj-number-3,ct,data)*twoscaling
										*pfchecknp(lfce[localj-number-1],lfce[localip]);
								}


							}


						}

						//now consider a helix stacking from the 3' sequence fragment:
						for (int localip=locali+minloop+1;localip<=number;localip++) {
							//first consider flush stacking

							localrarray+=
								w3[localip+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*penalty(localip,locali+1,ct,data)*
								ergcoaxflushbases(localj-number,locali,locali+1,localip,ct,data)*v->f(locali+1,localip)*twoscaling;


							if ((mod[locali+1]||mod[localip])) if (inc[ct->numseq[locali+2]][ct->numseq[localip-1]]&&notgu(locali+1,localip,ct)
								&&!(fce->f(locali+1,localip)&SINGLE)) {

								localrarray+=
									w3[localip+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*penalty(localip,locali+1,ct,data)*
									ergcoaxflushbases(localj-number,locali,locali+1,localip,ct,data)*v->f(locali+2,localip-1)
									*erg1(locali+1,localip,locali+2,localip-1,ct,data)*twoscaling;

							}

							//now consider an intervening nuc
							if (localj-number>1) {
								localrarray+=
									w3[localip+1]*w5[localj-number-2]*penalty(locali,localj,ct,data)*penalty(localip,locali+2,ct,data)*
									ergcoaxinterbases1(localj-number,locali,locali+2,localip,ct,data)*v->f(locali+2,localip)*twoscaling
									*pfchecknp(lfce[locali+1],lfce[localj-number-1]);


								if ((mod[locali+2]||mod[localip])) if (inc[ct->numseq[locali+3]][ct->numseq[localip-1]]&&notgu(locali+2,localip,ct)
									&&!(fce->f(locali+2,localip)&SINGLE)) {

									localrarray+=
										w3[localip+1]*w5[localj-number-2]*penalty(locali,localj,ct,data)*penalty(localip,locali+2,ct,data)*
										ergcoaxinterbases1(localj-number,locali,locali+2,localip,ct,data)*v->f(locali+3,localip-1)
										*erg1(locali+2,localip,locali+3,localip-1,ct,data)*twoscaling
										*pfchecknp(lfce[locali+1],lfce[localj-number-1]);

								}
							}


							//consider the other possibility for an intervening nuc
							localrarray+=
								w3[localip+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*penalty(localip-1,locali+2,ct,data)*
								ergcoaxinterbases2(localj-number,locali,locali+2,localip-1,ct,data)*v->f(locali+2,localip-1)*twoscaling
								*pfchecknp(lfce[locali+1],lfce[localip]);

							if ((mod[locali+2]||mod[localip-1])) if(inc[ct->numseq[locali+3]][ct->numseq[localip-2]]&&notgu(locali+2,localip-1,ct)
								&&!(fce->f(locali+2,localip-1)&SINGLE)) {

								localrarray+=
									w3[localip+1]*w5[localj-number-1]*penalty(locali,localj,ct,data)*penalty(localip-1,locali+2,ct,data)*
									ergcoaxinterbases2(localj-number,locali,locali+2,localip-1,ct,data)*v->f(locali+3,localip-2)
									*erg1(locali+2,localip-1,locali+3,localip-2,ct,data)*twoscaling*pfchecknp(lfce[locali+1],lfce[localip]);


							}



						}
						#endif //ifndef disablecoax
						#endif //simplembloop

					}//end of consider exterior loop closed by locali, localj






					//consider the multiloop closed by locali,localj
					if ((localj-locali)>(2*minloop+4)&&locali!=number) {
						#ifdef SIMPLEMBLOOP
						if (localj-1!=number&&ilocal!=number) {
							//5' and 3' dangling ends
							localrarray+=erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*
								erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1])*penalty(locali,localj,ct,data)*
            					wmb->f(locali+1,localj-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}
						else if (localj-1!=number) {//locali==number
							//5' dangling end
							localrarray+=
								erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1])*penalty(locali,localj,ct,data)*
            					wmb->f(locali+1,localj-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}
						if (locali!=number) {
							//3' dangling ends
							localrarray+=erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*
								penalty(locali,localj,ct,data)*
            					wmb->f(locali+1,localj-1)* data->eparam[5] * data->eparam[10]*twoscaling;
						}

						else {
							//no dangling end
							localrarray+=wmb->f(locali+1,localj-1)*
								data->eparam[5]*data->eparam[10]
            					*penalty(locali,localj,ct,data)*twoscaling;
						}
						#else //!SIMPLEMBLOOP


          				//no dangling ends on locali-localj pair:
						if (localj-1!=number) {
							localrarray+=wmb->f(locali+1,localj-1)*data->eparam[5]*data->eparam[10]
            					*penalty(locali,localj,ct,data)*twoscaling;


							//locali+1 dangles on locali-localj pair:

							if (locali+1!=number) localrarray+=erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1])*penalty(locali,localj,ct,data)*
            					wmb->f(locali+2,localj-1)* data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
						}
						if (localj-2!=number) {
							//localj-1 dangles
							if (localj!=(number+1))localrarray+=erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]) * penalty(locali,localj,ct,data) *
            					wmb->f(locali+1,localj-2) * data->eparam[5] * data->eparam[6] * data->eparam[10]*twoscaling;
							//both locali+1 and localj-1 dangle
							if ((locali+1!=number)&&(localj!=(number+1))) {
            					localrarray+=data->tstkm[ct->numseq[locali]][ct->numseq[localj]][ct->numseq[locali+1]][ct->numseq[localj-1]]*
											pfchecknp(lfce[locali+1],lfce[localj-1])*
											wmb->f(locali+2,localj-2) * data->eparam[5] * data->eparam[6] * data->eparam[6]* data->eparam[10]
											*penalty(locali,localj,ct,data)*twoscaling;
							}
						}





						//consider the coaxial stacking of a helix from locali to localj onto helix locali+1 or locali+2 to localip:
						#ifndef disablecoax //a flag to turn off coaxial stacking
						for (int localip=locali+1;(localip<localj);localip++) {
							//first consider flush stacking




							//conditions guarantee that the coaxial stacking isn't considering an exterior loop
							//if ((locali!=number)/*&&(locali+1!=number)*//*&&((localj>number)||(localip!=number)&&(localip+1!=number))&&(localj-1!=number)*/) {
							if (locali!=number&&localip!=number&&localj-1!=number) {
								if (inc[ct->numseq[locali+1]][ct->numseq[localip]]) {

									localrarray+=penalty(locali,localj,ct,data)*v->f(locali+1,localip)*
										penalty(locali+1,localip,ct,data)*data->eparam[5]
										*data->eparam[10]*data->eparam[10]*(w->f(localip+1,localj-1)+wmb->f(localip+1,localj-1))*ergcoaxflushbases(localj,locali,locali+1,localip,ct,data)
										*twoscaling;


									if((mod[locali+1]||mod[localip])) if (inc[ct->numseq[locali+2]][ct->numseq[localip-1]]&&notgu(locali+1,localip,ct)&&!(fce->f(locali+1,localip)&SINGLE)) {

										localrarray+=penalty(locali,localj,ct,data)*v->f(locali+2,localip-1)*
											penalty(locali+1,localip,ct,data)*data->eparam[5]
											*data->eparam[10]*data->eparam[10]*(w->f(localip+1,localj-1)+wmb->f(localip+1,localj-1))*ergcoaxflushbases(localj,locali,locali+1,localip,ct,data)
											*erg1(locali+1,localip,locali+2,localip-1,ct,data)*twoscaling;

									}

								}


								//check if locali+2 and localip can pair
								if (inc[ct->numseq[locali+2]][ct->numseq[localip]]) {
									//if ((localip<localj-1)&&(locali+2!=number)) {
									if (localip+2<localj-1&&locali+1!=number&&localip+1!=number) {
									//now consider an intervening nuc
										if ((localip+2<localj-1)/*&&(localj>number||localip+2!=number)*/) {
									
										
										
											localrarray+=penalty(locali,localj,ct,data)*v->f(locali+2,localip)*
											penalty(locali+2,localip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(localip+2,localj-1)+wmb->f(localip+2,localj-1))
											*ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data)*twoscaling*pfchecknp(lfce[locali+1],lfce[localip+1]);

										
										
											if((mod[locali+2]||mod[localip])) if (inc[ct->numseq[locali+3]][ct->numseq[localip-1]]&&notgu(locali+2,localip,ct)
												&&!(fce->f(locali+2,localip)&SINGLE)) {

												localrarray+=penalty(locali,localj,ct,data)*v->f(locali+3,localip-1)*
													penalty(locali+2,localip,ct,data)*data->eparam[5]
													*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
													(w->f(localip+2,localj-1)+wmb->f(localip+2,localj-1))
													*ergcoaxinterbases2(localj,locali,locali+2,localip,ct,data)
													*erg1(locali+2,localip,locali+3,localip-1,ct,data)*twoscaling*pfchecknp(lfce[locali+1],lfce[localip+1]);

											}

										

										}


										if (localip+1<localj-2&&localj-2!=number) {

										
												localrarray+=penalty(locali,localj,ct,data)*v->f(locali+2,localip)*
													penalty(locali+2,localip,ct,data)*data->eparam[5]
													*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
													(w->f(localip+1,localj-2)+wmb->f(localip+1,localj-2))
													*ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data)*twoscaling
													*pfchecknp(lfce[locali+1],lfce[localj-1]);

												if((mod[locali+2]||mod[localip])) if (inc[ct->numseq[locali+3]][ct->numseq[localip-1]]&&notgu(locali+2,localip,ct)
													&&!(fce->f(locali+2,localip)&SINGLE)) {

													localrarray+=penalty(locali,localj,ct,data)*v->f(locali+3,localip-1)*
														penalty(locali+2,localip,ct,data)*data->eparam[5]
														*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
														(w->f(localip+1,localj-2)+wmb->f(localip+1,localj-2))
														*ergcoaxinterbases1(localj,locali,locali+2,localip,ct,data)
														*erg1(locali+2,localip,locali+3,localip-1,ct,data)*twoscaling*pfchecknp(lfce[locali+1],lfce[localj-1]);

												}
										

										
										}




									}


								}//end of check whether locali+2 and loaclip can pair


							}


						}



						//consider the coaxial stacking of a helix from locali to localj onto helix localip to localj-2 or localj-1:
						for (int localip=localj-1;localip>locali;localip--) {


							//conditions guarantee that the coaxial stacking isn't considering an exterior loop
							//if ((locali!=number)&&(locali+1!=number)&&((localj>number)||(localip!=number)&&(localip-1!=number))&&(localj-1!=number)) {
							if (localj-1!=number&&localip-1!=number&&locali!=number) {
								//first consider flush stacking
								localrarray+=penalty(locali,localj,ct,data)*v->f(localip,localj-1)*
									penalty(localj-1,localip,ct,data)*data->eparam[5]
									*data->eparam[10]*data->eparam[10]*(w->f(locali+1,localip-1)+wmb->f(locali+1,localip-1))*ergcoaxflushbases(localip,localj-1,localj,locali,ct,data)
									*twoscaling;


								if((mod[localip]||mod[localj-1])) if(inc[ct->numseq[localip+1]][ct->numseq[localj-2]]&&notgu(localip,localj-1,ct)&&!(fce->f(localip,localj-1)&SINGLE)) {
									localrarray+=penalty(locali,localj,ct,data)*v->f(localip+1,localj-2)*
										penalty(localj-1,localip,ct,data)*data->eparam[5]
										*data->eparam[10]*data->eparam[10]*(w->f(locali+1,localip-1)+wmb->f(locali+1,localip-1))
										*ergcoaxflushbases(localip,localj-1,localj,locali,ct,data)
										*erg1(localip,localj-1,localip+1,localj-2,ct,data)*twoscaling;

								}





								if (localj-2!=number) {
									//now consider an intervening nuc
									//if ((localip>locali+1)&&(localj>number||localip-2!=number))
									if (localip-2>locali+1&&localip-2!=number) {
										localrarray+=penalty(locali,localj,ct,data)*v->f(localip,localj-2)*
											penalty(localj-2,localip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(locali+1,localip-2)+wmb->f(locali+1,localip-2))
											*ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data)*twoscaling*pfchecknp(lfce[localj-1],lfce[localip-1]);



										if((mod[localip]||mod[localj-2])) if(inc[ct->numseq[localip+1]][ct->numseq[localj-3]]&&notgu(localip,localj-2,ct)&&!(fce->f(localip,localj-2)&SINGLE)) {
											localrarray+=penalty(locali,localj,ct,data)*v->f(localip+1,localj-3)*
												penalty(localj-2,localip,ct,data)*data->eparam[5]
												*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
												(w->f(locali+1,localip-2)+wmb->f(locali+1,localip-2))
												*ergcoaxinterbases1(localip,localj-2,localj,locali,ct,data)
												*erg1(localip,localj-2,localip+1,localj-3,ct,data)*twoscaling*pfchecknp(lfce[localj-1],lfce[localip-1]);

										}
									}



									if ((localip-1>locali+2)&&locali+1!=number) {
										localrarray+=penalty(locali,localj,ct,data)*v->f(localip,localj-2)*
											penalty(localj-2,localip,ct,data)*data->eparam[5]
											*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
											(w->f(locali+2,localip-1)+wmb->f(locali+2,localip-1))
											*ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data)*twoscaling*pfchecknp(lfce[localj-1],lfce[locali+1]);

										if((mod[localip]||mod[localj-2])) if(inc[ct->numseq[localip+1]][ct->numseq[localj-3]]&&notgu(localip,localj-2,ct)
											&&!(fce->f(localip,localj-2)&SINGLE)) {
											localrarray+=penalty(locali,localj,ct,data)*v->f(localip+1,localj-3)*
												penalty(localj-2,localip,ct,data)*data->eparam[5]
												*data->eparam[6]*data->eparam[6]*data->eparam[10]*data->eparam[10]*
												(w->f(locali+2,localip-1)+wmb->f(locali+2,localip-1))
												*ergcoaxinterbases2(localip,localj-2,localj,locali,ct,data)
												*erg1(localip,localj-2,localip+1,localj-3,ct,data)*twoscaling*pfchecknp(lfce[localj-1],lfce[locali+1]);

										}
									}


								}



							}




						}
						#endif //ifndef disablecoax
						#endif //SIMPLEMBLOOP


						/*if (ct->intermolecular) {

            				//intermolecular, so consider wmb2,
							//don't add the multiloop penalties because this is a exterior loop

            				e[1] = min(e[1],wmb2->f(locali+1,localj-1) + penalty(locali,localj,ct,data)+infinity);


            				//locali+1 dangles on locali-localj pair:
            				if (locali!=number) e[2] = min(e[2],erg4(locali,localj,locali+1,1,ct,data,lfce[locali+1]) + penalty(locali,localj,ct,data) +
            					wmb2->f(locali+2,localj-1)+infinity);
            				//localj-1 dangles
            				if (localj!=(number+1)) e[3] = min(e[3],erg4(locali,localj,localj-1,2,ct,data,lfce[localj-1]) + penalty(locali,localj,ct,data) +
            					wmb2->f(locali+1,localj-2)+infinity);
            				//both locali+1 and localj-1 dangle
            				if ((locali!=number)&&(localj!=(number+1))) {
            					e[4] = min(e[4],
            					data->tstkm[ct->numseq[locali]][ct->numseq[localj]][ct->numseq[locali+1]][ct->numseq[localj-1]] +
									pfchecknp(lfce[locali+1],lfce[localj-1]) +
               						wmb2->f(locali+2,localj-2) + penalty(locali,localj,ct,data)+infinity);

							}




						}*/


					}//end of consider the multibranch loop closed by locali, localj



				}//end of large block for multibranch or exterior closure
				

			}//end of the large block for calculating V if all conditions are met
			//sub2:
			//Assign V:
//<<<<<<< pfunction.cpp
//			v->f(locali,localj) = localrarray;

			//cout << locali << " " << localj << " " << localrarray << "\n";
//=======
			v->f(locali,localj) = localrarray;
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
			//Apply constant (an euilibrium constant for formation of locali-localj pair), if being used:
//=======
			//cout << locali << " " << localj << " " << localrarray << "\n";

			//Apply constant (an euilibrium constant for formation of locali-localj pair), if being used:
//>>>>>>> 1.30
			if (ct->constant!=NULL) {

				//use ii and jj as the indices for accessing constant array:
//<<<<<<< pfunction.cpp

				if (locali>ct->GetSequenceLength()) ii = locali - ct->GetSequenceLength();
				else ii = locali;
				if (localj>ct->GetSequenceLength()) jj = localj - ct->GetSequenceLength();
				else jj = localj;

//=======
//				if (locali>ct->numofbases) ii = locali - ct->numofbases;
//				else ii = locali;
//				if (localj>ct->numofbases) jj = localj - ct->numofbases;
//				else jj = localj;
//>>>>>>> 1.30
				if (jj<ii) {
					p = jj;
      				jj = ii;
					ii = p;
				}

				v->f(locali,localj) = v->f(locali,localj)*ct->constant[jj][ii];

			}


			#ifndef SMP	
				//This code is for O(N^3) internal loops, and is not used for SMP calculations

				/*prefill curE[locali] and prev[locali] for the first two diognals as ll =4 and 5
					As d =10, only fill curE (ll=4, d=10)
					As d =11, only fill prevE (ll=4||5, d=11)
					As d>11, fill curE(ll=4||5, d>11)
					exchange curE and prevE after d >11  as curE[h][locali][dp]=curE=curE[h-2][locali+1][dp]
					(curE[locali][localj][ll]=curE[locali+1][localj-1][ll-2])
				*/


//<<<<<<< pfunction.cpp
//				if ((d-1)>=(minloop+3)||localj>number)
//				for (dp=d-3;dp>=((localj>number)?1:minloop+1);dp--) {
//					ll=d-dp-2;
					//calculate every localip>localip+1,localjp<localjp-1 when ll ==5 ||4
//					if (ll==4||ll==5) {
//						for (localip=locali+2;localip<=localj-2-dp;localip++) {
//							localjp=localip+dp;
//							if (inc[ct->numseq[localip]][ct->numseq[localjp]])
//							if ( (localip<=number&&localjp>number) || (localip<=number&&localj<=number) ) {
//								if(d==( (localj>number)?7:10 )||d>( (localj>number)?8:11 ))
								//fill the first diagonal of d and first two of ll for every larger d
//								{
//									curE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//              							fce->f(localj,localjp)) * v->f(localip,localjp);
//=======
				if ((d-1)>=(minloop+3)||localj>number)
				for (dp=d-3;dp>=((localj>number)?1:minloop+1);dp--) {
					ll=d-dp-2;
					//calculate every localip>localip+1,localjp<localjp-1 when ll ==5 ||4
					if (ll==4||ll==5) {
						for (int localip=locali+2;localip<=localj-2-dp;localip++) {
							int localjp=localip+dp;
							if (inc[ct->numseq[localip]][ct->numseq[localjp]])
							if ( (localip<=number&&localjp>number) || (localip<=number&&localj<=number) ) {
								if(d==( (localj>number)?7:10 )||d>( (localj>number)?8:11 ))
								//fill the first diagonal of d and first two of ll for every larger d
								{
									curE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp)) * v->f(localip,localjp);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
									//locali or localj is modified
//									if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))
//										curE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//              							fce->f(localj,localjp)) * v->f(localip+1,localjp-1) * erg1(localip,localjp,localip+1,localjp-1,ct,data);

//								}
//=======
									//locali or localj is modified
									if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))
										curE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp)) * v->f(localip+1,localjp-1) * erg1(localip,localjp,localip+1,localjp-1,ct,data);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//=======
								}
//>>>>>>> 1.30

//								else if ( d==((localj>number)?8:11) )  //fill the second diagonal of d
//								{
//				  					prevE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//              							fce->f(localj,localjp)) * v->f(localip,localjp);

//<<<<<<< pfunction.cpp
//				  					if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))
//=======
								else if ( d==((localj>number)?8:11) )  //fill the second diagonal of d
								{
				  					prevE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp)) * v->f(localip,localjp);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//			      						prevE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//              							fce->f(localj,localjp))* v->f(localip+1,localjp-1)* erg1(localip,localjp,localip+1,localjp-1,ct,data);
//								}
//=======
				  					if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))

			      						prevE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp))* v->f(localip+1,localjp-1)* erg1(localip,localjp,localip+1,localjp-1,ct,data);
								}
//>>>>>>> 1.30



							}

						}
					}

					//when size >=6 and <=30;
					//   else if (ll>=6&&ll<=(data->eparam[7]))
					else if (ll>=6&&ll<=maxinter)
					{
//<<<<<<< pfunction.cpp
						//calculate minimum curE[dp][locali] of 1 x (n-1) for next step's 2 x (n-2)
//						localip=locali+2;
//						localjp=localip+dp;
//						if ( (localip<=number&&localjp>number) || (localip<=number&&localj<=number) )
//						if (abs(localip-locali+localjp-localj)<=maxasym)
//						if (inc[ct->numseq[localip]][ct->numseq[localjp]])
//						{
//							curE[dp][locali] += erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
//              							fce->f(localj,localjp))*v->f(localip,localjp) ;
//=======
						//calculate minimum curE[dp][locali] of 1 x (n-1) for next step's 2 x (n-2)
						int localip=locali+2;
						int localjp=localip+dp;
						if ( (localip<=number&&localjp>number) || (localip<=number&&localj<=number) )
						if (abs(localip-locali+localjp-localj)<=maxasym)
						if (inc[ct->numseq[localip]][ct->numseq[localjp]])
						{
							curE[dp][locali] += erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp))*v->f(localip,localjp) ;
//>>>>>>> 1.30

							if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))
								curE[dp][locali]+=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp)) *
			     					v->f(localip+1,localjp-1) * erg1(localip,localjp,localip+1,localjp-1,ct,data);


						}

						localjp=localj-2;
						localip=localjp-dp;
						if ( (localip<=number&&localjp>number) || (localip<=number&&localj<=number) )
						if (abs(localip-locali+localjp-localj)<=maxasym)
						if (inc[ct->numseq[localip]][ct->numseq[localjp]])
						{

							curE[dp][locali] += erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp))*v->f(localip,localjp) ;

							if ((mod[localip]||mod[localjp])&&inc[ct->numseq[localip]][ct->numseq[localjp]]&&!(fce->f(localip,localjp)&SINGLE))
								curE[dp][locali] +=erg2in(locali,localj,localip,localjp,ct,data,fce->f(locali,localip),
              							fce->f(localj,localjp))*
			     					v->f(localip+1,localjp-1)*erg1(localip,localjp,localip+1,localjp-1,ct,data);


						}


					}
				}
			

			//also block propagation of interior loops that contain nucleotides that need to be double-stranded:
			if (lfce[locali]||lfce[localj]) for (dp=1;dp<=d;dp++) curE[dp][locali] = 0.0;
			#endif //SMP not defined

			if (fce->f(locali,localj)&PAIR)  {//force a pair between locali and localj
	  			w->f(locali,localj) = v->f(locali,localj)*data->eparam[10]*penalty(locali,localj,ct,data);
				wl->f(locali,localj) = v->f(locali,localj)*data->eparam[10]*penalty(locali,localj,ct,data);
				wmb->f(locali,localj) = 0.0;
				wmbl->f(locali,localj) = 0.0;
				wcoax->f(locali,localj) = 0.0;
				if (localj<=number) wca[locali][localj] = (PFPRECISION) 0;
	  			//goto sub3;
			}
			else {//not a forced pair

				////fill wmb:
				localrarray = (PFPRECISION) 0;
				if (((localj-locali)>(2*minloop+2))||localj>number) {


					#ifdef pfdebugmode
						ofstream dump;
						if (locali==10&&localj==22) {

							dump.open("dump.out");

						}
					#endif



					//also consider the coaxial stacking of two helixes in wv
					locale = (PFPRECISION) 0;
					#ifndef SIMPLEMBLOOP
					#ifndef disablecoax //a flag to diable coaxial stacking
					for (int localip=locali+minloop+1;localip<localj-minloop-1;localip++) {
						//first consider flush stacking

//<<<<<<< pfunction.cpp
//						if (localip!=number) {
//							localrarray+=v->f(locali,localip)*v->f(localip+1,localj)*penalty(locali,localip,ct,data)
//								*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data);

							

//							if ((mod[locali]||mod[localip]||mod[localip+1]||mod[localj])) {
//=======
						if (localip!=number) {
							localrarray+=v->f(locali,localip)*v->f(localip+1,localj)*penalty(locali,localip,ct,data)
								*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data);
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//								if ((mod[locali]||mod[localip])&&(mod[localip+1]||mod[localj])&&inc[ct->numseq[locali+1]][ct->numseq[localip-1]]
//									&&inc[ct->numseq[localip+2]][ct->numseq[localj-1]]&&notgu(locali,localip,ct)&&notgu(localip+1,localj,ct)
//										&&!(fce->f(localip+1,localj)&SINGLE)&&!(fce->f(locali,localip)&SINGLE)) {
//=======
							
//>>>>>>> 1.30

//<<<<<<< pfunction.cpp
//									localrarray+=v->f(locali+1,localip-1)*v->f(localip+2,localj-1)*penalty(locali,localip,ct,data)
//										*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data)
//										*erg1(locali,localip,locali+1,localip-1,ct,data)*erg1(localip+1,localj,localip+2,localj-1,ct,data);
//=======
							if ((mod[locali]||mod[localip]||mod[localip+1]||mod[localj])) {

								if ((mod[locali]||mod[localip])&&(mod[localip+1]||mod[localj])&&inc[ct->numseq[locali+1]][ct->numseq[localip-1]]
									&&inc[ct->numseq[localip+2]][ct->numseq[localj-1]]&&notgu(locali,localip,ct)&&notgu(localip+1,localj,ct)
										&&!(fce->f(localip+1,localj)&SINGLE)&&!(fce->f(locali,localip)&SINGLE)) {

									localrarray+=v->f(locali+1,localip-1)*v->f(localip+2,localj-1)*penalty(locali,localip,ct,data)
										*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data)
										*erg1(locali,localip,locali+1,localip-1,ct,data)*erg1(localip+1,localj,localip+2,localj-1,ct,data);
//>>>>>>> 1.30


								}

								if ((mod[locali]||mod[localip])&&inc[ct->numseq[locali+1]][ct->numseq[localip-1]]&&notgu(locali,localip,ct)&&!(fce->f(locali,localip)&SINGLE)) {

									localrarray+=v->f(locali+1,localip-1)*v->f(localip+1,localj)*penalty(locali,localip,ct,data)
										*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data)
										*erg1(locali,localip,locali+1,localip-1,ct,data);


								}

								if ((mod[localip+1]||mod[localj])&&inc[ct->numseq[localip+2]][ct->numseq[localj-1]]&&notgu(localip+1,localj,ct)&&!(fce->f(localip+1,localj)&SINGLE)) {


									localrarray+=v->f(locali,localip)*v->f(localip+2,localj-1)*penalty(locali,localip,ct,data)
										*penalty(localip+1,localj,ct,data)*ergcoaxflushbases(locali,localip,localip+1,localj,ct,data)
										*erg1(localip+1,localj,localip+2,localj-1,ct,data);

								}


							}



							if (localip+1!=number&&localj!=number+1) {
								if (!lfce[localip+1]&&!lfce[localj]) {
									//now consider an intervening mismatch
									locale+=v->f(locali,localip)*v->f(localip+2,localj-1)*penalty(locali,localip,ct,data)
										*penalty(localip+2,localj-1,ct,data)*ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data);

							

									if (mod[locali]||mod[localip]||mod[localip+2]||mod[localj-1]) {
										if ((mod[locali]||mod[localip])&&(mod[localip+2]||mod[localj-1])&&inc[ct->numseq[locali+1]][ct->numseq[localip-1]]
											&&inc[ct->numseq[localip+3]][ct->numseq[localj-2]]&&notgu(locali,localip,ct)&&notgu(localip+2,localj-1,ct)
												&&!(fce->f(locali,localip)&SINGLE)&&!(fce->f(localip+2,localj-1)&SINGLE)) {

											 locale+=v->f(locali+1,localip-1)*v->f(localip+3,localj-2)*penalty(locali,localip,ct,data)
												*penalty(localip+2,localj-1,ct,data)*ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data)
												*erg1(locali,localip,locali+1,localip-1,ct,data)*erg1(localip+2,localj-1,localip+3,localj-2,ct,data);


										}

										if ((mod[locali]||mod[localip])&&inc[ct->numseq[locali+1]][ct->numseq[localip-1]]&&notgu(locali,localip,ct)&&!(fce->f(locali,localip)&SINGLE)) {

											locale+=v->f(locali+1,localip-1)*v->f(localip+2,localj-1)*penalty(locali,localip,ct,data)
												*penalty(localip+2,localj-1,ct,data)*ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data)
												*erg1(locali,localip,locali+1,localip-1,ct,data);


										}

										if ((mod[localip+2]||mod[localj-1])&&inc[ct->numseq[localip+3]][ct->numseq[localj-2]]&&notgu(localip+2,localj-1,ct)&&!(fce->f(localip+2,localj-1)&SINGLE)) {


											locale+=v->f(locali,localip)*v->f(localip+3,localj-2)*penalty(locali,localip,ct,data)
												*penalty(localip+2,localj-1,ct,data)*ergcoaxinterbases2(locali,localip,localip+2,localj-1,ct,data)
												*erg1(localip+2,localj-1,localip+3,localj-2,ct,data);

										}
									}
								}

								if(!lfce[locali]&&!lfce[localip+1]&&locali!=number) {
									locale+=v->f(locali+1,localip)*v->f(localip+2,localj)*penalty(locali+1,localip,ct,data)
										*penalty(localip+2,localj,ct,data)*ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data);

							
									if (mod[locali+1]||mod[localip]||mod[localip+2]||mod[localj]) {
										if ((mod[locali+1]||mod[localip])&&(mod[localip+2]||mod[localj])&&inc[ct->numseq[locali+2]][ct->numseq[localip-1]]
											&&inc[ct->numseq[localip+3]][ct->numseq[localj-1]]&&notgu(locali+1,localip,ct)&&notgu(localip+2,localj,ct)
											&&!(fce->f(locali+1,localip)&SINGLE)&&!(fce->f(localip+2,localj)&SINGLE)	) {

											locale+=v->f(locali+2,localip-1)*v->f(localip+3,localj-1)*penalty(locali+1,localip,ct,data)
												*penalty(localip+2,localj,ct,data)*ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data)
												*erg1(locali+1,localip,locali+2,localip-1,ct,data)*erg1(localip+2,localj,localip+3,localj-1,ct,data);



										}
										if ((mod[locali+1]||mod[localip])&&inc[ct->numseq[locali+2]][ct->numseq[localip-1]]&&notgu(locali+1,localip,ct)&&!(fce->f(locali+1,localip)&SINGLE)) {

											locale+=v->f(locali+2,localip-1)*v->f(localip+2,localj)*penalty(locali+1,localip,ct,data)
												*penalty(localip+2,localj,ct,data)*ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data)
												*erg1(locali+1,localip,locali+2,localip-1,ct,data);


										}

										if ((mod[localip+2]||mod[localj])&&inc[ct->numseq[localip+3]][ct->numseq[localj-1]]&&notgu(localip+2,localj,ct)&&!(fce->f(localip+2,localj)&SINGLE)) {


											locale+=v->f(locali+1,localip)*v->f(localip+3,localj-1)*penalty(locali+1,localip,ct,data)
												*penalty(localip+2,localj,ct,data)*ergcoaxinterbases1(locali+1,localip,localip+2,localj,ct,data)
												*erg1(localip+2,localj,localip+3,localj-1,ct,data);

										}
									}
								}
							}
						}





					}
					#endif //ifndef disablecoax
					#endif //ifndef SIMPLEMBLOOP
					

					if (localj<=number) wca[locali][localj] = localrarray+locale;

					localrarray =(localrarray+locale*data->eparam[6]*data->eparam[6])*data->eparam[10]*data->eparam[10];
					wcoax->f(locali,localj) = localrarray;

					//search for an open bifurcation:
					for (int localk=locali+1;localk<localj;localk++) {
						//e = 0;
						if (localk!=number) {
							if (!lfce[locali]&&locali!=number)
								localrarray+=(wl->f(locali,localk)-wl->f(locali+1,localk)*data->eparam[6]*data->scaling+wcoax->f(locali,localk))*(wl->f(localk+1,localj)+wmbl->f(localk+1,localj));

							else localrarray+=(wl->f(locali,localk)+wcoax->f(locali,localk))*(wl->f(localk+1,localj)+wmbl->f(localk+1,localj));

						}
         			}




					if (locali!=number)
						if (!lfce[locali]) localrarray+=wmbl->f(locali+1,localj)*data->eparam[6]*data->scaling;

					wmbl->f(locali,localj) = localrarray;

					wmb->f(locali,localj) = localrarray;
					if (localj!=number+1)
						if (!lfce[localj]) wmb->f(locali,localj)+=wmb->f(locali,localj-1)*data->eparam[6]*data->scaling;


					/*if (ct->intermolecular) {
         				//intermolecular folding:




         			//search for an open bifurcation:
         			for (k=locali;k<=localj;k++) {

						if (k!=number) wmb2->f(locali,localj) = min(wmb2->f(locali,localj),w2->f(locali,k)+work2[k+1][jmt]);


         			}



					if (locali!=number)
						if (!(fce->f(locali,locali)&INTER)) wmb2->f(locali,localj) = min(wmb2->f(locali,localj) ,wmb2->f(locali+1,localj) );
						else  wmb2->f(locali,localj) = min(wmb2->f(locali,localj) ,wmb2->f(locali+1,localj) + data->init - infinity);

					if (localj!=number+1)
						if (!(fce->f(localj,localj)&INTER)) wmb2->f(locali,localj)  = min(wmb2->f(locali,localj) ,wmb2->f(locali,localj-1));
						else wmb2->f(locali,localj)  = min(wmb2->f(locali,localj) ,wmb2->f(locali,localj-1) +data->init-infinity);



					w2->f(locali,localj) = min(w2->f(locali,localj),wmb2->f(locali,localj) );

				}*/


				}

				//Compute w[locali][localj]
				if (localj>number||(localj-locali>minloop)) {

					#ifdef SIMPLEMBLOOP

  					//calculate the energy of localj stacked onto the pair of locali,localj-1
					if (localj!=N&&locali!=1) {
						//5' and 3' dangling ends
     					wl->f(locali,localj)= v->f(locali,localj)* data->eparam[10] *
     						erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1])*
							erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1])*penalty(locali,localj,ct,data);

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {

							wl->f(locali,localj)+= v->f(locali+1,localj-1) * data->eparam[10] *
     							erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1])*
								erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1])*penalty(locali,localj-1,ct,data)
								*erg1(locali,localj,locali+1,localj-1,ct,data);



						}


					}
					if (localj!=N) {//locali then ==1
						//3' dangling end
     					wl->f(locali,localj)= v->f(locali,localj)* data->eparam[10] *
     						erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1])*
							penalty(locali,localj,ct,data);

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {

							wl->f(locali,localj)+= v->f(locali+1,localj-1) * data->eparam[10] *
     							erg4(localj,locali,localj+1,1,ct,data,lfce[localj+1])*
								penalty(locali,localj-1,ct,data)
								*erg1(locali,localj,locali+1,localj-1,ct,data);



						}


					}
					else if (locali!=1) {//then localj==N
						//5' dangling end
     					wl->f(locali,localj)= v->f(locali,localj)* data->eparam[10] *
							erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1])*penalty(locali,localj,ct,data);

						if ((mod[locali]||mod[localj])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)&&!(fce->f(locali,localj)&SINGLE)) {

							wl->f(locali,localj)+= v->f(locali+1,localj-1) * data->eparam[10] *
								erg4(localj,locali,locali-1,2,ct,data,lfce[locali-1])*penalty(locali,localj-1,ct,data)
								*erg1(locali,localj,locali+1,localj-1,ct,data);



						}


					}

					else {//locali==1 and localj==N
						wl->f(locali,localj)= data->eparam[10]*v->f(locali,localj)*penalty(localj,locali,ct,data);

						if ((mod[locali]||mod[localj])) if (inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)) {

							wl->f(locali,localj)+= data->eparam[10]*v->f(locali+1,localj-1)*penalty(localj,locali,ct,data)*erg1(locali,localj,locali+1,localj-1,ct,data);

						}


					}
					#else  //!SIMPLEMBLOOP

					wl->f(locali,localj)= data->eparam[10]*v->f(locali,localj)*penalty(localj,locali,ct,data);

					if ((mod[locali]||mod[localj])) if (inc[ct->numseq[locali+1]][ct->numseq[localj-1]]&&notgu(locali,localj,ct)) {

						wl->f(locali,localj)+= data->eparam[10]*v->f(locali+1,localj-1)*penalty(localj,locali,ct,data)*erg1(locali,localj,locali+1,localj-1,ct,data);

					}



					if (locali!=number) {
      					//calculate the energy of locali stacked onto the pair of locali+1,localj

						wl->f(locali,localj)+= v->f(locali+1,localj)*data->eparam[10]*data->eparam[6]*
         					erg4(localj,locali+1,locali,2,ct,data,lfce[locali])*penalty(locali+1,localj,ct,data);

						if ((mod[locali+1]||mod[localj])) if(inc[ct->numseq[locali+2]][ct->numseq[localj-1]]&&notgu(locali+1,localj,ct)&&!(fce->f(locali+1,localj)&SINGLE)) {

							wl->f(locali,localj)+= v->f(locali+2,localj-1) * data->eparam[10] *data->eparam[6] *
         						erg4(localj,locali+1,locali,2,ct,data,lfce[locali])*penalty(locali+1,localj,ct,data)
								*erg1(locali+1,localj,locali+2,localj-1,ct,data);


						}

					}
					if (localj!=((number)+1)) {
      					//calculate the energy of localj stacked onto the pair of locali,localj-1
						if (localj!=1) {
         					wl->f(locali,localj)+= v->f(locali,localj-1)* data->eparam[10] * data->eparam[6] *
         						erg4(localj-1,locali,localj,1,ct,data,lfce[localj])*penalty(locali,localj-1,ct,data);

							if ((mod[locali]||mod[localj-1])) if(inc[ct->numseq[locali+1]][ct->numseq[localj-2]]&&notgu(locali,localj-1,ct)&&!(fce->f(locali,localj-1)&SINGLE)) {

								wl->f(locali,localj)+= v->f(locali+1,localj-2) * data->eparam[10] * data->eparam[6] *
         							erg4(localj-1,locali,localj,1,ct,data,lfce[localj])*penalty(locali,localj-1,ct,data)
									*erg1(locali,localj-1,locali+1,localj-2,ct,data);



							}


						}
					}
					if ((locali!=(number))&&(localj!=((number)+1))) {
      					//calculate locali and localj stacked onto the pair of locali+1,localj-1
						if (localj!=1&&!lfce[locali]&&!lfce[localj]) {
         					wl->f(locali,localj)+= v->f(locali+1,localj-1) *data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         						data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]][ct->numseq[localj]][ct->numseq[locali]]
								*penalty(localj-1,locali+1,ct,data);



							if ((mod[locali+1]||mod[localj-1])) if((localj-2>0)&&!(fce->f(locali+1,localj-1)&SINGLE)) {
								if(inc[ct->numseq[locali+2]][ct->numseq[localj-2]]&&notgu(locali+1,localj-1,ct)) {

									wl->f(locali,localj)+= v->f(locali+2,localj-2) * data->eparam[10] * (data->eparam[6]*data->eparam[6]) *
         								data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]][ct->numseq[localj]][ct->numseq[locali]]
										*penalty(localj-1,locali+1,ct,data)*erg1(locali+1,localj-1,locali+2,localj-2,ct,data);

								}
							}
						}
					}

					#endif  //SIMPLEMBLOOP


					if (locali!=number&&!lfce[locali]) {
         				//if (!(fce->f(locali,locali)&INTER))
					   //add a nuc to an existing loop:
         				wl->f(locali,localj)+=  wl->f(locali+1,localj)*data->eparam[6]*data->scaling;
            			//this is for when locali represents the center of an intermolecular linker:
						// else e[4] = w->f(locali+1,localj) + data->eparam[6] + infinity;
					}

					w->f(locali,localj) = wl->f(locali,localj);
					if (localj!=number+1&&!lfce[localj]) {
             			//if (!(fce->f(localj,localj)&INTER)) {
               			//add a nuc to an existing loop:
               			w->f(locali,localj)+= w->f(locali,localj-1) * data->eparam[6]*data->scaling;
					   //}
					   //else e[5] = w->f(locali,localj-1) + data->eparam[6] + infinity;

					}
				}

				/* if (ct->intermolecular) {

      				//wmb2[locali][localj%3] = infinity;
      				//keep track of w2:
					for (ii=1;ii<=5;ii++) e[ii] = 2*infinity;


					if (locali!=number) {
      					//calculate the energy of locali stacked onto the pair of locali+1,localj

         				e[1] = v->f(locali+1,localj) +
         					erg4(localj,locali+1,locali,2,ct,data,lfce[locali])+penalty(locali+1,localj,ct,data);

						if ((mod[locali+1]||mod[localj])) if(inc[ct->numseq[locali+2]][ct->numseq[localj-1]]) {

							e[1] = min(e[1],v->f(locali+2,localj-1) +
         						erg4(localj,locali+1,locali,2,ct,data,lfce[locali])+penalty(locali+1,localj,ct,data)
								+erg1(locali+1,localj,locali+2,localj-1,ct,data));

						}



         				e[4] = w2->f(locali+1,localj);

					}
      				if (localj!=((number)+1)) {
      				//calculate the energy of localj stacked onto the pair of locali,localj-1
         				if (localj!=1) {
         					e[2] = v->f(locali,localj-1)   +
         						erg4(localj-1,locali,localj,1,ct,data,lfce[localj])+penalty(locali,localj-1,ct,data);

							if ((mod[locali]||mod[localj-1])&&inc[ct->numseq[locali+1]][ct->numseq[localj-2]]) {

								e[2] = min(e[2],v->f(locali+1,localj-2) +
         							erg4(localj-1,locali,localj,1,ct,data,lfce[localj])+penalty(locali,localj-1,ct,data)
									+erg1(locali,localj-1,locali+1,localj-2,ct,data));

							}


               				e[5] = w2->f(locali,localj-1);

         				}
      				}
      				if ((locali!=(number))&&(localj!=((number)+1))) {
      					//calculate locali and localj stacked onto the pair of locali+1,localj-1
         				if (localj!=1) {
         					e[3] = v->f(locali+1,localj-1)   +
         						data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]]
											[ct->numseq[localj]][ct->numseq[locali]]
							+pfchecknp(lfce[locali+1],lfce[localj-1])
               				+penalty(localj-1,locali+1,ct,data);



							if ((mod[locali+1]||mod[localj-1])&&inc[ct->numseq[locali+2]][ct->numseq[localj-2]]) {

								e[3] = min(e[3],v->f(locali+2,localj-2) +
         							data->tstkm[ct->numseq[localj-1]][ct->numseq[locali+1]][ct->numseq[localj]][ct->numseq[locali]]
									+pfchecknp(lfce[locali+1],lfce[localj-1])
									+penalty(localj-1,locali+1,ct,data)+erg1(locali+1,localj-1,locali+2,localj-2,ct,data));

							}
         				}
      				}

					e[1] = min(e[1],(v->f(locali,localj)+penalty(localj,locali,ct,data)));

					if (mod[locali]||mod[localj]&&inc[ct->numseq[locali+1]][ct->numseq[localj-1]]) {

						e[1] = min((v->f(locali+1,localj-1)+penalty(localj,locali,ct,data)+erg1(locali,localj,locali+1,localj-1,ct,data)),e[1]);

					}


      				w2->f(locali,localj) = min(e[1],e[2]);
      				w2->f(locali,localj) = min(w2->f(locali,localj),e[3]);
      				w2->f(locali,localj) = min(w2->f(locali,localj),e[4]);
      				w2->f(locali,localj) = min(w2->f(locali,localj),e[5]);




				}*/



			}//end of else "not a forced pair"



		

		//sub3:







		


			

		}
		

		#ifdef pfdebugmode
		if (twoscaling>PFMAX||twoscaling<PFMIN) {
			//look for positions that will underflow


			ofstream *ufout = new ofstream();
			ufout->open("c:/underflow_alert.out",ios::app);
			*ufout<<"locali= "<<locali<<" localj= "<<localj<<" twoscaling= "<<twoscaling<<"\n";
			ufout->close();
			delete ufout;

		}
		#endif//pfdebugmode
	}//end for locali


	for (i=((h<=(number-1))?1:(2*number-h));i<=((h<=(number-1))?(number-h):number);i++){
		j=i+d;

		
		//cout << i << " " << j << " " << v->f(i,j) << "\n";

		//check to see if any of the 2-D arrays are about to go out of bounds
		//(not checking wca[][],curE[][],prevE[][] although they need to be rescaled too)
		if (v->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (w->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wl->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wcoax->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmb->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (wmbl->f(i,j)>PFMAX) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
			twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
		}
		else if (v->f(i,j)<PFMIN&&v->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (w->f(i,j)<PFMIN&&w->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wl->f(i,j)<PFMIN&&wl->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,curE,prevE,wca,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wcoax->f(i,j)<PFMIN&&wcoax->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmb->f(i,j)<PFMIN&&wmb->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
		else if (wmbl->f(i,j)<PFMIN&&wmbl->f(i,j)>0) {
			rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
			twoscaling = twoscaling*SCALEUP*SCALEUP;
		}
	}

		//Compute w5[i], the energy of the best folding from 1->i, and

      		//w3[i], the energy of the best folding from i-->numofbases

	if (h<=(number-1)) {
		i = 1;
		j=i+d;

//<<<<<<< pfunction.cpp
//	if (h<=(number-1)) {
//		i = 1;
//		j=i+d;

//		if (j<=minloop+1) {
//			if (lfce[j]) w5[j]= (PFPRECISION) 0;
//			else  w5[j] = w5[j-1]*data->scaling;
//		}
//=======
		if (j<=minloop+1) {
			if (lfce[j]) w5[j]= (PFPRECISION) 0;
			else  w5[j] = w5[j-1]*data->scaling;
		}
//>>>>>>> 1.30

		else {
      		if (lfce[j]) rarray = (PFPRECISION) 0;

			else rarray = w5[j-1]*data->scaling;



      		for (k=0;k<=(j-4);k++) {

				#ifdef SIMPLEMBLOOP
				//5' and 3' dangling ends
				if (k>0&&j!=N) {
					rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
					*erg4(j,k+1,k,2,ct,data,lfce[k])
					*v->f(k+1,j)*penalty(j,k+1,ct,data);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
							*erg4(j,k+1,k,2,ct,data,lfce[k])*v->f(k+2,j-1)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
					}
				}
				//5' dangling end
				if (k>0) {//j==N
					rarray+=w5[k]*erg4(j,k+1,k,2,ct,data,lfce[k])
					*v->f(k+1,j)*penalty(j,k+1,ct,data);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*erg4(j,k+1,k,2,ct,data,lfce[k])*v->f(k+2,j-1)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
					}
				}
				//3' dangling end
				if (j!=N) {//k==0
					rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
					*v->f(k+1,j)*penalty(j,k+1,ct,data);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*erg4(j,k+1,j+1,1,ct,data,lfce[j+1])
							*v->f(k+2,j-1)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
					}
				}
				//no dangling ends
				else {//k==0 and j==N
					rarray+=w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data);

					if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

						rarray+=w5[k]*v->f(k+2,j-1)
							*penalty(j-1,k+1,ct,data)*erg1(k+1,j,k+2,j-1,ct,data);
					}
				}
				#else //!SIMPLEMBLOOP

      			rarray+=w5[k]*v->f(k+1,j)*penalty(j,k+1,ct,data);

				if ((mod[k+1]||mod[j])) if(inc[ct->numseq[k+2]][ct->numseq[j-1]]&&notgu(k+1,j,ct)&&!(fce->f(k+1,j)&SINGLE)) {

					rarray+=w5[k]*v->f(k+2,j-1)*penalty(j,k+1,ct,data)
						*erg1(k+1,j,k+2,j-1,ct,data);
				}



				rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+2,j)*penalty(j,k+2,ct,data);

				if((mod[k+2]||mod[j])) if(inc[ct->numseq[k+3]][ct->numseq[j-1]]&&notgu(k+2,j,ct)
					&&!(fce->f(k+2,j)&SINGLE)) {
					rarray+=w5[k]*erg4(j,k+2,k+1,2,ct,data,lfce[k+1])*v->f(k+3,j-1)
						*penalty(j,k+2,ct,data)*erg1(k+2,j,k+3,j-1,ct,data);

				}


         		rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+1,j-1)*penalty(j-1,k+1,ct,data);

				if ((mod[k+1]||mod[j-1])) if(inc[ct->numseq[k+2]][ct->numseq[j-2]]&&notgu(k+1,j-1,ct)&&!(fce->f(k+1,j-1)&SINGLE)) {

					rarray+=w5[k]*erg4(j-1,k+1,j,1,ct,data,lfce[j])*v->f(k+2,j-2)
						*penalty(j-1,k+1,ct,data)*erg1(k+1,j-1,k+2,j-2,ct,data);
				}



				rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
								*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+2,j-1)*
								penalty(j-1,k+2,ct,data);

				if ((mod[k+2]||mod[j-1])) if(inc[ct->numseq[k+3]][ct->numseq[j-2]]&&notgu(k+2,j-1,ct)&&!(fce->f(k+2,j-1)&SINGLE)) {

					rarray+=w5[k]*data->tstack[ct->numseq[j-1]][ct->numseq[k+2]][ct->numseq[j]][ct->numseq[k+1]]
								*pfchecknp(lfce[j],lfce[k+1]) * v->f(k+3,j-2)*
								penalty(j-1,k+2,ct,data)*erg1(k+2,j-1,k+3,j-2,ct,data);

				}








				rarray+=w5[k]*wca[k+1][j];

				#endif  //SIMPLEMBLOOP
			}


			w5[j] = rarray;

			//check to see if w5 is about to go out of bounds:
			if (w5[j]>PFMAX) {
				rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
				twoscaling = twoscaling*SCALEDOWN*SCALEDOWN;
			}
			else if (w5[j]<PFMIN&&w5[j]>0) {
				rescale(i,j,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
				twoscaling=twoscaling*SCALEUP*SCALEUP;

			}
		}

	}//end if (h<=number-1)

	if (h==number-1) {

		//w3[0] = 0;
		//w3[number+1] = 0;
		for (ii=(number);ii>=(number-minloop);ii--) {    //number+1 ... number-minloop
      		if (lfce[ii]) w3[ii] = (PFPRECISION) 0;
			else w3[ii]=w3[ii+1]*data->scaling;
		}
		//w3[i]=0;
   		for (ii=((number)-minloop-1);ii>=1;ii--) {

      		if (lfce[ii]) rarray = (PFPRECISION) 0;

   			else rarray = w3[ii+1]*data->scaling;



			for (k=((number)+1);k>=(ii+4);k--) {
				#ifdef SIMPLEMBLOOP
				if (ii>1&&k!=N+1) {
					//5' and 3' dangling ends
					rarray+=v->f(ii,k-1)*erg4(k-1,ii,k,1,ct,data,lfce[k])
						*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])
						*penalty(k-1,ii,ct,data)*w3[k];

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray+=v->f(ii+1,k-2)*erg4(k-1,ii,k,1,ct,data,lfce[k])
						*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])*
						penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

					}
				}
				else if (k!=N+1) {//i==1
					//3' dangling end
					rarray+=v->f(ii,k-1)*erg4(k-1,ii,k,1,ct,data,lfce[k])
						*penalty(k-1,ii,ct,data)*w3[k];

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray+=v->f(ii+1,k-2)*erg4(k-1,ii,k,1,ct,data,lfce[k])*
						penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

					}
				}
				else if (ii>1) {//k==N+1
					//5' dangling end
					rarray+=v->f(ii,k-1)
						*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])
						*penalty(k-1,ii,ct,data)*w3[k];

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray+=v->f(ii+1,k-2)
						*erg4(k-1,ii,ii-1,2,ct,data,lfce[ii-1])*
						penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

					}
				}
				else {//ii==1&&k==N+1
					//No dangling ends
					rarray+=v->f(ii,k-1)
						*penalty(k-1,ii,ct,data)*w3[k];

					if((mod[ii]||mod[k-1]))if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
						rarray+=v->f(ii+1,k-2)*
						penalty(k-1,ii,ct,data)*w3[k]*erg1(ii,k-1,ii+1,k-2,ct,data);

					}
				}
				#else //!SIMPLEMBLOOP


      			rarray+=v->f(ii,k-1)*w3[k]*penalty(k-1,ii,ct,data);

				if((mod[ii]||mod[k-1])) if(inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {
					rarray+=v->f(ii+1,k-2)*w3[k]*penalty(k-1,ii,ct,data)*erg1(ii,k-1,ii+1,k-2,ct,data);

				}


				rarray+=v->f(ii+1,k-1)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*penalty(k-1,ii+1,ct,data) * w3[k];

				if((mod[ii+1]||mod[k-1])) if(inc[ct->numseq[ii+2]][ct->numseq[k-2]]&&notgu(ii+1,k-1,ct)&&!(fce->f(ii+1,k-1)&SINGLE)) {

					rarray+=v->f(ii+2,k-2)*erg4(k-1,ii+1,ii,2,ct,data,lfce[ii])*
						penalty(k-1,ii+1,ct,data) *w3[k]*erg1(ii+1,k-1,ii+2,k-2,ct,data);

				}


				rarray+=v->f(ii,k-2)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])*penalty(k-2,ii,ct,data)*w3[k];

				if((mod[ii]||mod[k-2]))if(inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {
					rarray+=v->f(ii+1,k-3)*erg4(k-2,ii,k-1,1,ct,data,lfce[k-1])*
						penalty(k-2,ii,ct,data)*w3[k]*erg1(ii,k-2,ii+1,k-3,ct,data);

				}

				if (!lfce[ii]&&!lfce[k-1]) {
					rarray+=v->f(ii+1,k-2)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
						[ct->numseq[k-1]][ct->numseq[ii]]
						*w3[k]*
						penalty(k-2,ii+1,ct,data);



					if((mod[ii+1]||mod[k-2]))if(inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
						rarray+=v->f(ii+2,k-3)*data->tstack[ct->numseq[k-2]][ct->numseq[ii+1]]
							[ct->numseq[k-1]][ct->numseq[ii]]
							*pfchecknp(lfce[k-1],lfce[ii])*w3[k]*
							penalty(k-2,ii+1,ct,data)*erg1(ii+1,k-2,ii+2,k-3,ct,data);


					}
				}

				//also consider coaxial stacking:
				#ifndef disablecoax //a flag to disable coaxial stacking
				for (ip=k+minloop+1;ip<=number+1;ip++) {


					//first consider flush stacking:
					rarray+=v->f(ii,k-1)*v->f(k,ip-1)*w3[ip]*
						penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
						ergcoaxflushbases(ii,k-1,k,ip-1,ct,data);

					if(mod[ii]||mod[k-1]||mod[k]||mod[ip-1]) {

						if ((mod[ii]||mod[k-1])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]
							&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii,k-1,ct)&&notgu(k,ip-1,ct)
							&&!(fce->f(ii,k-1)&SINGLE)&&!(fce->f(k,ip-1)&SINGLE)) {

							rarray+=v->f(ii+1,k-2)*v->f(k+1,ip-2)*w3[ip]*
								penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
								*erg1(ii,k-1,ii+1,k-2,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);
						}
						if((mod[ii]||mod[k-1])&&inc[ct->numseq[ii+1]][ct->numseq[k-2]]&&notgu(ii,k-1,ct)&&!(fce->f(ii,k-1)&SINGLE)) {

							rarray+=v->f(ii+1,k-2)*v->f(k,ip-1)*w3[ip]*
								penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
								*erg1(ii,k-1,ii+1,k-2,ct,data);

						}

						if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {

							rarray+=v->f(ii,k-1)*v->f(k+1,ip-2)*w3[ip]*
								penalty(ii,k-1,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxflushbases(ii,k-1,k,ip-1,ct,data)
								*erg1(k,ip-1,k+1,ip-2,ct,data);
						}

					}


					//now consider an intervening mismatch:
					if (ii>0&&!lfce[ii]&&!lfce[k-1]) {
						rarray+=v->f(ii+1,k-2)*v->f(k,ip-1)*w3[ip]*
							penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
							ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data);

						if(mod[ii+1]||mod[k-2]||mod[k]||mod[ip-1]){

							if((mod[ii+1]||mod[k-2])&&(mod[k]||mod[ip-1])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]
								&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(ii+1,k-2,ct)&&notgu(k,ip-1,ct)
								&&!(fce->f(k,ip-1)&SINGLE)&&!(fce->f(ii+1,k-2)&SINGLE)){
								rarray+=v->f(ii+2,k-3)*v->f(k+1,ip-2)*w3[ip]*
									penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
									ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
									*erg1(ii+1,k-2,ii+2,k-3,ct,data)*erg1(k,ip-1,k+1,ip-2,ct,data);

							}

							if((mod[ii+1]||mod[k-2])&&inc[ct->numseq[ii+2]][ct->numseq[k-3]]&&notgu(ii+1,k-2,ct)&&!(fce->f(ii+1,k-2)&SINGLE)) {
								rarray+=v->f(ii+2,k-3)*v->f(k,ip-1)*w3[ip]*
								penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
								ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
								*erg1(ii+1,k-2,ii+2,k-3,ct,data);

							}
							if((mod[k]||mod[ip-1])&&inc[ct->numseq[k+1]][ct->numseq[ip-2]]&&notgu(k,ip-1,ct)&&!(fce->f(k,ip-1)&SINGLE)) {
								rarray+=v->f(ii+1,k-2)*v->f(k+1,ip-2)*w3[ip]*
									penalty(ii+1,k-2,ct,data)*penalty(k,ip-1,ct,data)*
									ergcoaxinterbases1(ii+1,k-2,k,ip-1,ct,data)
									*erg1(k,ip-1,k+1,ip-2,ct,data);


							}


						}

					}
					if (!lfce[k-1]&&!lfce[ip-1]) {

						rarray+=v->f(ii,k-2)*v->f(k,ip-2)*w3[ip]*
							penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
							ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data);

						if (mod[ii]||mod[k-2]||mod[k]||mod[ip-2]) {

							if ((mod[ii]||mod[k-2])&&(mod[k]||mod[ip-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]
								&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(ii,k-2,ct)&&notgu(k,ip-2,ct)
								&&!(fce->f(ii,k-2)&SINGLE)&&!(fce->f(k,ip-2)&SINGLE)) {

								rarray+=v->f(ii+1,k-3)*v->f(k+1,ip-3)*w3[ip]*
									penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
									ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
									*erg1(ii,k-2,ii+1,k-3,ct,data)*erg1(k,ip-2,k+1,ip-3,ct,data);
							}

							if ((mod[ii]||mod[k-2])&&inc[ct->numseq[ii+1]][ct->numseq[k-3]]&&notgu(ii,k-2,ct)&&!(fce->f(ii,k-2)&SINGLE)) {

								rarray+=v->f(ii+1,k-3)*v->f(k,ip-2)*w3[ip]*
									penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
									ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
									*erg1(ii,k-2,ii+1,k-3,ct,data);
							}

							if ((mod[k]||mod[ip-2])&&inc[ct->numseq[k+1]][ct->numseq[ip-3]]&&notgu(k,ip-2,ct)&&!(fce->f(k,ip-2)&SINGLE)) {

								rarray+=v->f(ii,k-2)*v->f(k+1,ip-3)*w3[ip]*
									penalty(ii,k-2,ct,data)*penalty(k,ip-2,ct,data)*
									ergcoaxinterbases2(ii,k-2,k,ip-2,ct,data)
									*erg1(k,ip-2,k+1,ip-3,ct,data);
							}

						}
					}

				}
				#endif //ifndef disablecoax
				#endif //SIMPLEMBLOOP


			}

			w3[ii] = rarray;
			//check to see if w5 is about to go out of bounds:
			if (w3[ii]>PFMAX) {
				rescale(1,number,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEDOWN);
				twoscaling=twoscaling*SCALEDOWN*SCALEDOWN;
			}
			else if (w3[ii]<PFMIN&&w3[ii]>0) {
				rescale(1,number,ct,data,v,w,wl,wcoax,wmb,wmbl,w5,w3,wca,curE,prevE,SCALEUP);
				twoscaling = twoscaling*SCALEUP*SCALEUP;
			}
		}
	}

	#ifndef SMP
		//This is O(N^3) internal loop code that is not used with SMP
		if (d>(j>number?8:11))
		{
			tempE=curE;
			curE=prevE;
			prevE=tempE;
		}
		if(d> (j>number?7:10))
		
		for (dp=1;dp<=d-1;dp++) {
			for (i=((h<=(number-2))?1:(2*number-h-1));i<=((h<=(number-2))?(number-h-1):number);i++){
				if (i<number) curE[dp][i]=curE[dp][i+1];
			}
		}
	#endif




}

for(ii=0;ii<=number;ii++) {
	delete[] wca[ii];

	#ifndef SMP
	delete[] curE[ii];
	delete[] prevE[ii];
	#endif

}
delete[] wca;
#ifndef SMP
delete[] curE;
delete[] prevE;
#endif



#ifdef timer
timeout << time(NULL)<<"\n";
timeout << time(NULL) - seconds;
timeout.close();
#endif

//////////////////////////
//output V, W, WMB, and W2V:
#if defined (pfdebugmode)
	ofstream foo;
	foo.open("arrays.out");
	foo << "i" << "\t"<<"j"<<"\t"<<"v->f(i,j)"<<"\t"<<"w->f(i,j)"<<"\t"<<"wmb->f(i,j)\twmbl->f(i,j)\twcoax->f(i,j)"<<"\t"<<"wl->f(i,j)"<<"\t"<<"v->f(j,i+number)"<<"\t"<<"w->f(j,i+number)"<<"\t"<<"wmb->f(j,i+number)"<<"\t"<<"wl->f(j,i+number)"<<"\t"<<"wmbl->f(j,i+numer)\twcoax->f(j,i+number)"<<"\n";
	for (j=1;j<=number;j++) {
		for (i=1;i<=j;i++) {

			foo << i << "\t"<<j<<"\t"<<v->f(i,j)<<"\t"<<w->f(i,j)<<"\t"<<wmb->f(i,j)<<"\t"<<wmbl->f(i,j)<<"\t"<<wcoax->f(i,j)<<"\t"<<wl->f(i,j)<<"\t"<<v->f(j,i+number)<<"\t"<<w->f(j,i+number)<<"\t"<<wmb->f(j,i+number)<<"\t"<<wl->f(j,i+number)<<"\t"<<wmbl->f(j,i+number)<<"\t"<<wcoax->f(j,i+number)<<"\n";

		}
	}

	foo <<"\n\n\n";
	foo << "i" << "\t" << "w5[i]" << "\t" << "w3[i]" << "\n";
	for (i=0;i<=number;i++) foo << i << "\t" << w5[i] << "\t" << w3[i] << "\n";

	foo.close();

#endif

}




//This function cacluates a partition function for the sequence in CT
//If quickQ == true, return the partition function value in Q
//	otherwise, save the partial partition functions in the datafile named save
//If updates on progress are unwanted, set update=NULL
void pfunction(structure* ct,pfdatatable* data, TProgressDialog* update, char* save, bool quickQ, PFPRECISION *Q)
{



int i,j;
bool *lfce,*mod;//[maxbases+1][maxbases+1];
PFPRECISION *w5,*w3;
int number;



#ifdef equiout
	ofstream kout;
	kout.open("k.out");
	kout << "sequence length = "<<ct->GetSequenceLength()<<"\n";
	for (i=1;i<=ct->GetSequenceLength();i++) {
		kout << tobase(ct->numseq[i]);
		if (i%20==0) kout << "\n";
	}
	kout << "\n";
	kout.close();
#endif





//array *v,*w;
//inc is an array that saves time by showing which bases can pair before
//	erg is called

//number is the number of bases being folded
//v[i][j] is the best energy for subsequence i to j when i and j are paired
//	for i<j<n, it is the interior framgment between nucleotides i and j
//	for i<n<j, it is the exterior structure fragmnet that is 5' to j-n and 3' to i
//w[i][j] is the best energy for subsequence i to j

number = (ct->GetSequenceLength());//place the number of bases in an integer

//scaling is a per nucleotide scale factor for which W and V are divided
//This is necessary to keep the doubles from overflowing:

//scaling = 0.2;//this factor assumes about 1 kcal/mol/base



//allocate space for the v and w arrays:
pfunctionclass w(number);
pfunctionclass v(number);
pfunctionclass wmb(number);
pfunctionclass wl(number);
pfunctionclass wmbl(number);
pfunctionclass wcoax(number);
forceclass fce(number);

if (ct->intermolecular) {
	//take advantage of templating to prevent intramolecular base pairs

	ct->allocatetem();
	for (i=1;i<ct->inter[0];i++) {
		for (j=i+1;j<=ct->inter[2];j++) {

			ct->tem[j][i]=false;

		}
	}
	for (i=ct->inter[2]+1;i<ct->GetSequenceLength();i++) {
		for (j=i+1;j<=ct->GetSequenceLength();j++) {

			ct->tem[j][i]=false;

		}
	}


}

//This code converts the SHAPE array of data to equilibrium constants.  This is
//needed for the partition function.  NOTE, however, there is no going back so a structure
//that has been used for partition functions cannot then be used to predict a structure
//by free energy minimization. This is a compromise for efficiency, but clearly something
//that could cause a problem for programmers.


if (ct->shaped) {
	for (i=1;i<=2*ct->GetSequenceLength();i++) ct->SHAPE[i]=boltzman(ct->SHAPE[i], data->temp);


}

//add a second array for intermolecular folding:

/*if (ct->intermolecular) {
	w2 = new pfunctionclass(number);
	wmb2 = new pfunctionclass(number);




}

else {

	wmb2=NULL;
	w2=NULL;

}*/




lfce = new bool [2*number+1];
mod = new bool [2*number+1];

for (i=0;i<=2*number;i++) {
	lfce[i] = false;
	mod[i] = false;
}


//Register modified nucleotides in the mod array for fast access
for (i=0;i<ct->GetNumberofModified();i++) {

	if (ct->GetModified(i)!=1&&ct->GetModified(i)!=ct->GetSequenceLength()) {
		mod[ct->GetModified(i)]=true;
		mod[ct->GetModified(i)+ct->GetSequenceLength()]=true;
	}
}




w5 = new PFPRECISION [number+1];
w3 = new PFPRECISION [number+2];
//wca = new PFPRECISION *[number+1];



//The next section handles the case where base pairs are not
				//not allowed to form between nucs more distant
				//than ct->GetPairingDistanceLimit()
if (ct->DistanceLimited()) {

	if (!ct->templated) ct->allocatetem();

	for (j=minloop+2;j<=ct->GetSequenceLength();j++) {
		for (i=1;i<j;i++) {
			if (j-i>=ct->GetPairingDistanceLimit()) ct->tem[j][i]=false;
		}
	}



}



calculatepfunction(ct,data,update,save,quickQ,Q,&w,&v,&wmb,&wl,&wmbl,&wcoax,&fce,w5,w3,mod,lfce);




if (save!=0) {
	writepfsave(save,ct,w5,w3,&v,&w,&wmb,&wl,&wmbl,&wcoax,&fce,mod,lfce,data);
}

if (quickQ) *Q = w5[ct->GetSequenceLength()];


delete[] lfce;
delete[] mod;




delete[] w5;
delete[] w3;


/*if (ct->intermolecular) {
	delete w2;
	delete wmb2;
}*/



return;
}


//writepfsave writes a save file with partition function data.
void writepfsave(char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
			 forceclass *fce, bool *mod, bool *lfce, pfdatatable *data) {

	int i,j,k,l,m,n,o,p;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};

	ofstream sav(filename,ios::binary);

	//write the save file information so that the sequence can be re-folded,
		//include thermodynamic data to prevent traceback errors

	//start with a version number of the savefile
	int vers=pfsaveversion;
	write(&sav,&vers); //save a version of the save file

	//start with structure information
	int SequenceLength = ct->GetSequenceLength();
	write(&sav,&(SequenceLength));
	write(&sav,&(ct->intermolecular));
	write(&sav,&data->scaling);

	int constraint;

	//Write out the number of forced pairs
	constraint = ct->GetNumberofPairs();
	write(&sav,&(constraint));
	//Write out the list of forced pairs
	for (i=0;i<ct->GetNumberofPairs();i++) {
		constraint = ct->GetPair5(i);
		write(&sav,&(constraint));
		constraint = ct->GetPair5(i);
		write(&sav,&(constraint));
	}
	for (i=0;i<=ct->GetSequenceLength();i++) {

		write(&sav,&(ct->hnumber[i]));
		sav.write(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->numseq[i]));

	//Write out the number of nucleotides forced double-stranded
	constraint = ct->GetNumberofDoubles();
	write(&sav,&(constraint));

	//Write out the nucleotides forced double stranded
	for (i=0;i<ct->GetNumberofDoubles();i++) {
		
		constraint = ct->GetDouble(i);
		write(&sav,&(constraint));

	}


	if (ct->intermolecular) {
		for (i=0;i<3;i++) write(&sav,&(ct->inter[i]));

	}

	//Write the number of nucleotides forced single stranded
	constraint = ct->GetNumberofSingles();
	write(&sav,&(constraint));

	//Write the list of nucleotides forced single stranded
	for (i=0;i<ct->GetNumberofSingles();i++) {
		
		constraint = ct->GetSingle(i);
		write(&sav,&(constraint));

	}

	//Write the number of nucleotides that are modified 
	constraint = ct->GetNumberofModified();
	write(&sav,&(constraint));

	//Write the list of modified nucleotides
	for (i=0;i<ct->GetNumberofModified();i++) {
		
		constraint = ct->GetModified(i);
		write(&sav,&(constraint));

	}

	//Write the number of Us in GU pairs
	constraint = ct->GetNumberofGU();
	write(&sav,&(constraint));

	//Write the list of Us in GU pairs
	for (i=0;i<ct->GetNumberofGU();i++) {
		
		constraint = ct->GetGUpair(i);
		write(&sav,&(constraint));

	}
	
	string label=ct->GetSequenceLabel();
	write(&sav,&(label));

	write(&sav,&(ct->templated));
	if (ct->templated) {
		for (i=0;i<=ct->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) write(&sav,&(ct->tem[i][j]));

		}

	}

	//write the SHAPE data (for pseudo-free energy constraints)
	write(&sav,&(ct->shaped));
	if (ct->shaped) {
		for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->SHAPE[i]));
		for (i=0;i<=2*ct->GetSequenceLength();i++) write(&sav,&(ct->SHAPEss[i]));

	}


	//now write the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		write(&sav,&(w3[i]));
		write(&sav,&(w5[i]));
//<<<<<<< pfunction.cpp
//		for (j=0;j<=ct->GetSequenceLength();j++) {
//			write(&sav,&(v->dg[i][j]));
//			write(&sav,&(w->dg[i][j]));
//			write(&sav,&(wmb->dg[i][j]));
//			write(&sav,&(wmbl->dg[i][j]));
//			write(&sav,&(wl->dg[i][j]));
//			write(&sav,&(wcoax->dg[i][j]));
//=======
		for (j=0;j<=ct->GetSequenceLength();j++) {
			write(&sav,&(v->dg[i][j+i]));
			write(&sav,&(w->dg[i][j+i]));
			write(&sav,&(wmb->dg[i][j+i]));
			write(&sav,&(wmbl->dg[i][j+i]));
			write(&sav,&(wl->dg[i][j+i]));
			write(&sav,&(wcoax->dg[i][j+i]));
//>>>>>>> 1.30
			writesinglechar(&sav,&(fce->dg[i][j]));
			//if (ct->intermolecular) {
			//	write(&sav,&(w2->dg[i][j]));
			//	write(&sav,&(wmb2->dg[i][j]));

			//}


		}


	}


	// I should replace this with functionality to convert .pfs files to .txt files -- rhiju.
	// this is inelegant -- quick way to get at base pair probabilities...
	if ( false ) { //turn this off before CVS checkin!
	  string bpp_file( "bpp.txt" );
	  ofstream bpp_out( bpp_file.c_str() );
	  for (i=1;i<=ct->GetSequenceLength();i++) {
	    for (j=1;j<=ct->GetSequenceLength();j++) {
	      bpp_out << ' ' << calculateprobability(i,j,v,w5,ct,data,lfce,mod,data->scaling,fce);
	    }
	    bpp_out << endl;
	  }
	  bpp_out.close();
	}


	write(&sav,&(w3[ct->GetSequenceLength()+1]));
	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		write(&sav,&(lfce[i]));
		write(&sav,&(mod[i]));

	}





	//now write the thermodynamic data:
	write(&sav,&(data->temp));
	for (i=0;i<5;i++) write(&sav,&(data->poppen[i]));
	write(&sav,&(data->maxpen));
	for (i=0;i<11;i++) write(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		write(&sav,&(data->inter[i]));
		write(&sav,&(data->bulge[i]));
		write(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					write(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<6;l++) {
					write(&sav,&(data->stack[i][j][k][l]));
					write(&sav,&(data->tstkh[i][j][k][l]));
					write(&sav,&(data->tstki[i][j][k][l]));
					write(&sav,&(data->coax[i][j][k][l]));
					write(&sav,&(data->tstackcoax[i][j][k][l]));
					write(&sav,&(data->coaxstack[i][j][k][l]));
					write(&sav,&(data->tstack[i][j][k][l]));
					write(&sav,&(data->tstkm[i][j][k][l]));
					write(&sav,&(data->tstki23[i][j][k][l]));
					write(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							write(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) write(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										write(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}


						}
					}
				}
			}
		}
	}
	write(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		write(&sav,&(data->itloop[i]));
		write(&sav,&(data->tloop[i]));

	}
	write(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		write(&sav,&(data->itriloop[i]));
		write(&sav,&(data->triloop[i]));

	}
	write(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		write(&sav,&(data->ihexaloop[i]));
		write(&sav,&(data->hexaloop[i]));

	}
	write(&sav,&(data->auend));
	write(&sav,&(data->gubonus));
	write(&sav,&(data->cint));
	write(&sav,&(data->cslope));
	write(&sav,&(data->c3));
	write(&sav,&(data->efn2a));
	write(&sav,&(data->efn2b));
	write(&sav,&(data->efn2c));
	write(&sav,&(data->init));
	write(&sav,&(data->mlasym));
	write(&sav,&(data->strain));
	write(&sav,&(data->prelog));
	write(&sav,&(data->singlecbulge));
	write(&sav,&(data->maxintloopsize));



	sav.close();


}



////////////////////////////////////////////////////////////////////////
//pfunctionclass encapsulates the large 2-d arrays of w and v, used by the
//	partition function

      //the constructor allocates the space needed by the arrays
pfunctionclass::pfunctionclass(int size) {
	//zero indicates whether the array should be set to zero as opposed
		//to being set to infinity, it is false by default



	infinite = (PFPRECISION) 0;

    Size = size;
    /*register*/ int i,j;
    dg = new PFPRECISION *[size+1];

	for (i=0;i<=(size);i++)  {
   		dg[i] = new PFPRECISION [size+1];
   	}
    for (i=0;i<=size;i++) {
         for (j=0;j<size+1;j++) {

             dg[i][j] = (PFPRECISION) 0;
         }
    }
	

	//Move pointers for faster access:
	for (i=1;i<=size;++i) {
		dg[i]-=i;

	}

}

//the destructor deallocates the space used
pfunctionclass::~pfunctionclass() {


	int i;

    for (i=0;i<=Size;i++) {
		//restore the pointer position before deleting
		dg[i]+=i;
        delete[] dg[i];
    }
     delete[] dg;
}

      //f is an integer function that references the correct element of the array
//inline PFPRECISION &pfunctionclass::f(int i, int j) {



//   if (i>j) {
 //       return infinite;
//    }
 //  else if (i>Size) return f(i-Size,j-Size);
 //  else return dg[i][j-i];

//}

//When considering mismatch at the end of a helix, consult this function to check
//	whether the nucs are required to pair
PFPRECISION pfchecknp(bool lfce1,bool lfce2) {

	if (lfce1||lfce2) return 0;
	else return 1;
}



pfdatatable::pfdatatable() {
}

pfdatatable::pfdatatable(datatable *data, PFPRECISION Scaling, PFPRECISION  Temp) {
	//the partition function datatable needs to be initialized from the datatable
	int i,j,k,l,m,n,o,p;

	scaling = Scaling;

	//store the temperature in the pfdatatable
	temp = Temp;

	for (i=1;i<5;i++) poppen[i]=boltzman(data->poppen[i],temp);
	maxpen = boltzman(data->maxpen,temp);
	for (i=1;i<11;i++) eparam[i] = boltzman(data->eparam[i],temp);
	maxintloopsize = data->eparam[7];
	for (i=1;i<31;i++) {
		inter[i] = boltzman(data->inter[i],temp)*pow(scaling,i+2);
		bulge[i] = boltzman(data->bulge[i],temp)*pow(scaling,i+2);
		hairpin[i] = boltzman(data->hairpin[i],temp)*pow(scaling,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=1;l<3;l++) {
					#ifdef SIMPLEMBLOOP
					//In the case of simple multibranch loops, dangles no longer
						//occupy a position in the sequence and do not need a scaling
						//factor
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],temp);
					#else //!SIMPLEMBLOOP
					dangle[i][j][k][l] = boltzman(data->dangle[i][j][k][l],temp)*scaling;

					#endif
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=boltzman(data->stack[i][j][k][l],temp)*pow(scaling,2);
					tstkh[i][j][k][l]=boltzman(data->tstkh[i][j][k][l],temp);
					tstki[i][j][k][l]=boltzman(data->tstki[i][j][k][l],temp);
					coax[i][j][k][l]=boltzman(data->coax[i][j][k][l],temp);
					tstackcoax[i][j][k][l]=boltzman(data->tstackcoax[i][j][k][l],temp)*pow(scaling,2);
					coaxstack[i][j][k][l] = boltzman(data->coaxstack[i][j][k][l],temp);
					tstack[i][j][k][l]=boltzman(data->tstack[i][j][k][l],temp)*pow(scaling,2);
					tstkm[i][j][k][l]=boltzman(data->tstkm[i][j][k][l],temp)*pow(scaling,2);
					tstki23[i][j][k][l]=boltzman(data->tstki23[i][j][k][l],temp);
					tstki1n[i][j][k][l]=boltzman(data->tstki1n[i][j][k][l],temp);
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=boltzman(data->iloop11[i][j][k][l][m][n],temp)*pow(scaling,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=boltzman(data->iloop21[i][j][k][l][m][n][o],temp)*pow(scaling,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=boltzman(data->iloop22[i][j][k][l][m][n][o][p],temp)*pow(scaling,6);
								}
							}


						}
					}
				}
			}
		}
	}
	numoftloops = data->numoftloops;
	for (i=1;i<=data->numoftloops;i++) {
		itloop[i]=data->tloop[i][0];
		tloop[i] = boltzman(data->tloop[i][1],temp)*pow(scaling,6);

	}
	numoftriloops=data->numoftriloops;
	for (i=1;i<=data->numoftriloops;i++) {
		itriloop[i] = data->triloop[i][0];
		triloop[i] = boltzman(data->triloop[i][1],temp)*pow(scaling,5);
	}
	numofhexaloops=data->numofhexaloops;
	for (i=1;i<=data->numofhexaloops;i++) {
		ihexaloop[i]=data->hexaloop[i][0];
		hexaloop[i] = boltzman(data->hexaloop[i][1],temp)*pow(scaling,8);

	}
	auend = boltzman(data->auend,temp);
	gubonus = boltzman(data->gubonus,temp);
	cint = boltzman(data->cint,temp);
	cslope = boltzman(data->cslope,temp);
	c3 = boltzman(data->c3,temp);
	efn2a = boltzman(data->efn2a,temp);
	efn2b = boltzman(data->efn2b,temp);
	efn2c = boltzman(data->efn2c,temp);
	init = boltzman(data->init,temp);
	mlasym = boltzman(data->mlasym,temp);
	strain = boltzman(data->strain,temp);
	prelog = data->prelog/conversionfactor;
	singlecbulge = boltzman(data->singlecbulge,temp);





}





//calculate the energy of stacked base pairs

//oriented by:
//5' i ip 3'
//   |  |
//   j jp

PFPRECISION erg1(int i,int j,int ip,int jp,structure *ct, pfdatatable *data)
{

		PFPRECISION energy;

		 if ((i==(ct->GetSequenceLength()))||(j==((ct->GetSequenceLength())+1))) {
      		//this is not allowed because n and n+1 are not cavalently attached
			energy = (PFPRECISION) 0;
		}
		else {
      		energy = data->stack[(ct->numseq[i])][(ct->numseq[j])]
				[(ct->numseq[ip])][(ct->numseq[jp])]*data->eparam[1];

				if (ct->shaped) {
					energy=energy*ct->SHAPE[i];
					energy=energy*ct->SHAPE[j];
					energy=energy*ct->SHAPE[ip];
					energy=energy*ct->SHAPE[jp];
				}

				if ( ct->experimentalPairBonusExists ) {

				    energy = energy * ct->EX[i][j] * ct->EX[ip][jp];

				}
		}

		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k base pair stack ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}


//calculate the energy of a bulge/internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,
	char a, char b)
{

	int size,size1,size2, lopsid, count;
	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return 0;
      }



      size1 = ip-i-1;
		size2 = j - jp - 1;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return 0;//the loop contains a nuc that
      		//should be double stranded



	}


      //a typical internal or bulge loop:
		//size1 = ip-i-1;
		//size2 = j - jp - 1;
		if (size1==0||size2==0) {//bulge loop


			size = size1+size2;
			if (size==1) {
				count = 1;
				energy = data->stack[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[ip]][ct->numseq[jp]]
						*data->bulge[size]*data->eparam[2]/(data->scaling*data->scaling);
				if (size1==1)  {

					//count the number of alternative bulges that exist:

					//k = i;
					//while (ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k--;
					//}
					//k=ip;
					//while(ct->numseq[k]==ct->numseq[i+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[i+1]==2&&(ct->numseq[i+2]==2||ct->numseq[i]==2)) energy= energy*data->singlecbulge;

				}

				else {
					//size2 == 1

					//count the number of alternative bulges that exist:

					//k = jp;
					//while (ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k--;
					//}
					//k=j;
					//while(ct->numseq[k]==ct->numseq[jp+1]) {
					//	count++;
					//	k++;
					//}
					//give bonus to C adjacent to single C bulge
					if (ct->numseq[j-1]==2&&(ct->numseq[j-2]==2||ct->numseq[j]==2)) energy=energy*data->singlecbulge;

				}
				//do not apply a correction for the number of equivalent states because
					//the bulge can move to adjacent sites in the partition function calc
				//energy-= (int) (rt*conversionfactor* log ((double) count));
			}
			else if (size>30) {

				loginc = ((data->prelog)*log(PFPRECISION ((size)/30.0)));
				energy = data->bulge[30]*exp(-loginc/(RKC*data->temp))*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data)*pow(data->scaling,size-30);

			}
			else {
         		energy = data->bulge[size]*data->eparam[2];
				energy = energy*penalty(i,j,ct,data)*penalty(jp,ip,ct,data);
			}
		}
		else {//internal loop
			size = size1 + size2;
			lopsid = abs(size1-size2);

			if (size>30) {

				loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));
				if (size1==1||size2==1) {
            		energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp)) *data->eparam[3]*
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);

				}

				else {
					energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->tstki[ct->numseq[jp]][ct->numseq[ip]]
						[ct->numseq[jp+1]][ct->numseq[ip-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp))*data->eparam[3] *
						max(data->maxpen,
						pow(data->poppen[min(2,min(size1,size2))],lopsid))
						*pow(data->scaling,size-30);
				}
			}
			else if ((size1==2)&&(size2==2))//2x2 internal loop
			    energy = data->iloop22[ct->numseq[i]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[jp]]
					[ct->numseq[i+1]][ct->numseq[i+2]]
					[ct->numseq[j-1]][ct->numseq[j-2]];


			else if ((size1==1)&&(size2==2)) {//2x1 internal loop
				energy = data->iloop21[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]]
					[ct->numseq[j-1]][ct->numseq[jp+1]][ct->numseq[ip]][ct->numseq[jp]];


			}
			else if ((size1==2)&&(size2==1)) {//1x2 internal loop
				energy = data->iloop21[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]]
					[ct->numseq[ip-1]][ct->numseq[i+1]][ct->numseq[j]][ct->numseq[i]];

			}

			else if (size==2) //a single mismatch

				energy = data->iloop11[ct->numseq[i]][ct->numseq[i+1]][ct->numseq[ip]]
					[ct->numseq[j]][ct->numseq[j-1]][ct->numseq[jp]];
			else if (size1==1||size2==1) { //this loop is lopsided
         	//this is a 1xn loop:
				energy = data->tstki1n[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
						data->tstki1n[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
						data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}


			else if ((size1==2&&size2==3)||(size1==3&&size2==2)) {
			//this is a 2x3 loop
				energy = data->tstki23[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki23[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));


			}
			else {



         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
					data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					data->inter[size] * data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
			}
		}
		#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		kout << "k internal/bulge ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<energy<< "\n";
		kout.close();
		#endif

		return energy;
}

//calculate the energy of the intorior part of a internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2in(int i,int j,int ip,int jp,structure *ct, pfdatatable *data,char a,char b)
{

	int size,size1,size2, lopsid;
	PFPRECISION energy;

	if ((a>0)||(b>0)) {
      	if ((a&DUBLE)||(b&DUBLE)) return 0;//the loop contains a nuc that
      		//should be double stranded



	}

   	if (((i<=(ct->GetSequenceLength()))&&(ip>(ct->GetSequenceLength())))||((
      	jp<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength())))) {
         //A loop cannot contain the ends of the sequence

         return 0;
      }
      size1 = ip-i-1;
	  size2 = j - jp - 1;
      size = size1 + size2;
	  lopsid = abs(size1-size2);
	  if (size1 != 0 && size2 !=0)
	 	{
        energy = 	data->tstki[ct->numseq[jp]][ct->numseq[ip]][ct->numseq[jp+1]][ct->numseq[ip-1]] *
					 data->eparam[3] *
					max(data->maxpen,pow(data->poppen[min(2,min(size1,size2))],lopsid));
		}

		return energy;
}


//calculate the energy of the exterior part of internal loop
//where i is paired to j; ip is paired to jp; ip > i; j > jp
PFPRECISION erg2ex(int i,int j,int size, structure *ct, pfdatatable *data)
{

	PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/
				if (size>30) {

				loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));
				energy = data->tstki[ct->numseq[i]][ct->numseq[j]]
						[ct->numseq[i+1]][ct->numseq[j-1]]*
						data->inter[30]* exp(-loginc/(RKC*data->temp))*pow(data->scaling,size-30);
			}
						else
         		energy = data->tstki[ct->numseq[i]][ct->numseq[j]][ct->numseq[i+1]][ct->numseq[j-1]] *
						 data->inter[size]	;
		return energy;
}




//calculate the energy of a hairpin loop:
#ifndef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#else
PFPRECISION erg3indirect(int i,int j,structure *ct, pfdatatable *data,char dbl) {
#endif
int size,count,key,k;
PFPRECISION energy,loginc;
	/* size,size1,size2 = size of a loop
		energy = energy calculated
		loginc = the value of a log used in large hairpin loops
	*/

	if ((i<=(ct->GetSequenceLength()))&&(j>(ct->GetSequenceLength()))) {
      	//A hairpin cannot contain the ends of the sequence

         return 0;
    }

	if (dbl&DUBLE) return 0;//the loop contains a base that should be
      										//double stranded

    else if (dbl&INTER) {//intermolecular interaction
      	//intermolecular "hairpin" free energy is that of intermolecular
         //	initiation times the stacked mismatch

         energy =  data->tstack[ct->numseq[i]][ct->numseq[j]]
         	[ct->numseq[i+1]][ct->numseq[j-1]]*penalty(i,j,ct,data)*pow(data->scaling,j-i-3);


		//Add other states that can exist, i.e. 5' and 3' dangles:
		if (ct->numseq[i+1]!=5&&ct->numseq[j-1]!=5) {
			//Add a 3' dangling end

			energy+=erg4(i,j,i+1,1,ct,data,false)*pow(data->scaling,j-i-2)*penalty(i,j,ct,data);


			//Add a 5' dangling end

			energy+=erg4(i,j,j-1,2,ct,data,false)*pow(data->scaling,j-i-2)*penalty(i,j,ct,data);


			//Add the case where nothing stacks

			energy+=pow(data->scaling,j-i-1)*penalty(i,j,ct,data);

		}
		else if (ct->numseq[i+1]!=5||ct->numseq[j-1]!=5) {

			//Add the case where nothing stacks

			energy+=pow(data->scaling,j-i-1)*penalty(i,j,ct,data);

		}


         return data->init*energy;
    }




		size = j-i-1;



		if (size>30) {

			loginc = ((data->prelog)*log((PFPRECISION ((size))/30.0)));



			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[30]*exp(-loginc/(RKC*data->temp))*data->eparam[4]*pow(data->scaling,size-30);
		}
		else if (size<3) {
      		energy = data->hairpin[size]*data->eparam[4];
				if (ct->numseq[i]==4||ct->numseq[j]==4) energy = energy*exp(-.6/(RKC*data->temp));
		}
		else if (size==4) {

			key = (ct->numseq[j])*3125 + (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftloops;count++) {
					if (key==data->itloop[count]) return data->tloop[count];
			}
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}
		else if (size==3) {

			key = (ct->numseq[j])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numoftriloops;count++) {
				if (key==data->itriloop[count]) return data->triloop[count];
			}

			energy =	data->hairpin[size] * data->eparam[4]
         	*penalty(i,j,ct,data);
		}
		else if (size==6) {
			key = (ct->numseq[j])*78125 + (ct->numseq[i+6])*15625 + (ct->numseq[i+5])*3125
				+ (ct->numseq[i+4])*625 +
				(ct->numseq[i+3])*125 + (ct->numseq[i+2])*25+(ct->numseq[i+1])*5+(ct->numseq[i]);
			for (count=1;count<=data->numofhexaloops;count++) {
				if (key==data->ihexaloop[count]) return data->hexaloop[count];
			}

			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}

		else {
			energy = data->tstkh[ct->numseq[i]][ct->numseq[j]]
				[ct->numseq[i+1]][ct->numseq[j-1]]
				* data->hairpin[size] * data->eparam[4];
		}




		//check for GU closeure preceded by GG
      if (ct->numseq[i]==3&&ct->numseq[j]==4) {
      	if ((i>2&&i<ct->GetSequenceLength())||(i>ct->GetSequenceLength()+2))
       		if (ct->numseq[i-1]==3&&ct->numseq[i-2]==3) {

         		energy = energy * data->gubonus;



         	}
      }

      //check for an oligo-c loop

      for (k=1;(k<=size);k++) {
       	if (ct->numseq[i+k] != 2) return energy;//this is not an oligo-C loop
      }
      //this is a poly c loop so penalize
      if (size==3) return (energy *data->c3);
      else return (energy * data->cint * pow(data->cslope,size));


}

#ifdef equiout
PFPRECISION erg3(int i,int j,structure *ct, pfdatatable *data,char dbl) {
	PFPRECISION energy;

	energy = erg3indirect(i,j,ct, data,dbl,temp);
	ofstream kout;
	kout.open("k.out",ofstream::app);
	kout << "k hairpin ("<<i<<","<<j<<") ="<<energy<< "\n";
	kout.close();
	return energy;
}
#endif

//calculate the energy of a dangling end:
PFPRECISION erg4(int i,int j,int ip,int jp,structure *ct, pfdatatable *data, bool lfce)
{

//dangling base
		// jp = 1 => 3' dangle
		// jp = 2 => 5' dangle



      if (lfce) return 0;//stacked nuc should be double stranded

	  //commented out 11/8/99
      //if (ip==5) return 0;//dangling nuc is an intermolecular linker
	    #ifdef equiout
			ofstream kout;
			kout.open("k.out",ofstream::app);
			kout << "k dangle ("<<i<<","<<j<<","<<ip<<","<<jp<<") ="<<data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp]<< "\n";
			kout.close();
		#endif

		return data->dangle[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][jp];

}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty(int i,int j,structure* ct, pfdatatable *data) {

	#ifdef equiout
		ofstream kout;
		kout.open("k.out",ofstream::app);
		if (ct->numseq[i]==4||ct->numseq[j]==4) kout << "k penalty ("<<i<<","<<j<<") ="<<data->auend<< "\n";
		else kout << "k penalty ("<<i<<","<<j<<") ="<<1.0<< "\n";
		kout.close();
	#endif

	if (ct->numseq[i]==4||ct->numseq[j]==4)
   	return data->auend;
	else return 1;//no end penalty




}

//this function calculates whether a terminal pair i,j requires the end penalty
PFPRECISION penalty2(int i,int j, pfdatatable *data) {


   if (i==4||j==4)
   	return data->auend;
   else return 1;//no end penalty


}



//this function calculates the free energy for coaxial stacking of pair i-j onto ip-jp
//	require that the backbone continues directly between j and ip without a nucleotide in
//	a canonical pair
//k indicates the intervening nuc in intervening mismatch that is not sandwiched between
//	j and ip
//#ifdef equiout
//PFPRECISION ergcoaxindirect(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//#else
PFPRECISION ergcoaxflushbases(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif

		//flush stacking

		return data->coax[ct->numseq[i]][ct->numseq[j]][ct->numseq[ip]][ct->numseq[jp]];





}

PFPRECISION ergcoaxinterbases1(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
		//coaxial stacking with an intervening mismatch
		//(k==i-1)

			return data->tstackcoax[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[i-1]] *
				data->coaxstack[ct->numseq[j+1]][ct->numseq[i-1]]
				[ct->numseq[ip]][ct->numseq[jp]];


}

PFPRECISION ergcoaxinterbases2(int i, int j, int ip, int jp, structure *ct, pfdatatable *data) {
//#endif
	//coaxial stacking with an intervening mismatch
	/*(k==jp+1) {*/

			return data->tstackcoax[ct->numseq[jp]][ct->numseq[ip]]
				[ct->numseq[jp+1]][ct->numseq[ip-1]] *
				data->coaxstack[ct->numseq[j]][ct->numseq[i]]
				[ct->numseq[j+1]][ct->numseq[jp+1]];

}

//#ifdef equiout
//PFPRECISION ergcoax(int i, int j, int ip, int jp, int k, structure *ct, pfdatatable *data) {
//PFPRECISION energy;

//	energy = ergcoaxindirect(i,j,ip,jp,k,ct,data);

//	ofstream kout;

//	if (((i<j&&j<ip&&ip<jp)||(i>j&&j>ip&&ip>jp))&&k!=i&&k!=j&&k!=ip&&k!=jp) {
//		kout.open("k.out",ofstream::app);
//		kout << "k coax ("<<i<<","<<j<<","<<ip<<","<<jp<<","<<k<<") ="<<energy<< "\n";
//		kout.close();
//	}
//	return energy;
//}
//#endif


void readpfsave(const char *filename, structure *ct,
			 PFPRECISION *w5, PFPRECISION *w3,
			 pfunctionclass *v, pfunctionclass *w, pfunctionclass *wmb, pfunctionclass *wl, pfunctionclass *wmbl, pfunctionclass *wcoax,
			 forceclass *fce, PFPRECISION *scaling, bool *mod, bool *lfce, pfdatatable *data) {
	 int i,j,k,l,m,n,o,p;
	ifstream sav(filename,ios::binary);
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};//a mask array indicating the identity of canonical pairs
	int vers;

	//read the save file

	//read the file version first
	read(&sav,&vers);

	//start with structure information
	int SequenceLength;
	read(&sav,&(SequenceLength));
	//ct->numofbases = SequenceLength;
	read(&sav,&(ct->intermolecular));
	read(&sav,scaling);

	data->scaling=*scaling;

	int constraint,constraint2,numberofconstraints;

	//Read information about the forced pairs
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		read(&sav,&(constraint));
		read(&sav,&(constraint2));

		ct->AddPair(constraint,constraint2);
	}
	for (i=0;i<=ct->GetSequenceLength();i++) {

		read(&sav,&(ct->hnumber[i]));
		sav.read(&(ct->nucs[i]),1);

	}

	for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->numseq[i]));


	//Read information about nucleotides forced to be double stranded.
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		
		read(&sav,&(constraint));

		ct->AddDouble(constraint);

	}
	if (ct->intermolecular) {
		for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));

	}

	//Read information about nucleotides not allowed to pair
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		
		read(&sav,&(constraint));

		ct->AddSingle(constraint);

	}

	//Read information about nucleotides that are accessible to chemical modification:
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		
		read(&sav,&(constraint));

		ct->AddModified(constraint);

	}

	//Read information about Us in GU pairs:
	read(&sav,&(numberofconstraints));
	for (i=0;i<numberofconstraints;i++) {
		
		read(&sav,&(constraint));

		ct->AddGUPair(constraint);

	}

	string label;
	read(&sav,&(label));
	ct->SetSequenceLabel(label);

	read(&sav,&(ct->templated));
	if (ct->templated) {
		ct->allocatetem();
		for (i=0;i<=ct->GetSequenceLength();i++) {
			for (j=0;j<=i;j++) read(&sav,&(ct->tem[i][j]));

		}

	}

	read(&sav,&(ct->shaped));
	if (ct->shaped) {
		ct->SHAPE = new double [2*ct->GetSequenceLength()+1];
		for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->SHAPE[i]));
		ct->SHAPEss = new double [2 * ct->GetSequenceLength() + 1];
		for (i=0;i<=2*ct->GetSequenceLength();i++) read(&sav,&(ct->SHAPEss[i]));

	}


	//now read the array class data for v, w, and wmb:
	for (i=0;i<=ct->GetSequenceLength();i++) {
		read(&sav,&(w3[i]));
		read(&sav,&(w5[i]));
//<<<<<<< pfunction.cpp
//		for (j=0;j<=ct->GetSequenceLength();j++) {
//			read(&sav,&(v->dg[i][j]));
//			read(&sav,&(w->dg[i][j]));
//			read(&sav,&(wmb->dg[i][j]));
//			read(&sav,&(wmbl->dg[i][j]));
//			read(&sav,&(wl->dg[i][j]));
//			read(&sav,&(wcoax->dg[i][j]));
//=======
		for (j=0;j<=ct->GetSequenceLength();j++) {
			read(&sav,&(v->dg[i][j+i]));
			read(&sav,&(w->dg[i][j+i]));
			read(&sav,&(wmb->dg[i][j+i]));
			read(&sav,&(wmbl->dg[i][j+i]));
			read(&sav,&(wl->dg[i][j+i]));
			read(&sav,&(wcoax->dg[i][j+i]));
//>>>>>>> 1.30
			readsinglechar(&sav,&(fce->dg[i][j]));
			//if (ct->intermolecular) {
			//	read(&sav,&(w2->dg[i][j]));
			//	read(&sav,&(wmb2->dg[i][j]));

			//}


		}


	}

	read(&sav,&(w3[ct->GetSequenceLength()+1]));
	for (i=0;i<=2*ct->GetSequenceLength();i++) {
		read(&sav,&(lfce[i]));
		read(&sav,&(mod[i]));

	}





	//now read the thermodynamic data:
	read(&sav,&(data->temp));
	for (i=0;i<5;i++) read(&sav,&(data->poppen[i]));
	read(&sav,&(data->maxpen));
	for (i=0;i<11;i++) read(&sav,&(data->eparam[i]));
	for (i=0;i<31;i++) {
		read(&sav,&(data->inter[i]));
		read(&sav,&(data->bulge[i]));
		read(&sav,&(data->hairpin[i]));

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					read(&sav,&(data->dangle[i][j][k][l]));
				}
				for (l=0;l<6;l++) {
					read(&sav,&(data->stack[i][j][k][l]));
					read(&sav,&(data->tstkh[i][j][k][l]));
					read(&sav,&(data->tstki[i][j][k][l]));
					read(&sav,&(data->coax[i][j][k][l]));
					read(&sav,&(data->tstackcoax[i][j][k][l]));
					read(&sav,&(data->coaxstack[i][j][k][l]));
					read(&sav,&(data->tstack[i][j][k][l]));
					read(&sav,&(data->tstkm[i][j][k][l]));
					read(&sav,&(data->tstki23[i][j][k][l]));
					read(&sav,&(data->tstki1n[i][j][k][l]));
					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							read(&sav,&(data->iloop11[i][j][k][l][m][n]));
							for (o=0;o<6;o++) {
								if (inc[i][j]&&inc[n][o]) read(&sav,&(data->iloop21[i][j][k][l][m][n][o]));
								for (p=0;p<6;p++) {
									if (inc[i][k]&&inc[j][l])
										read(&sav,&(data->iloop22[i][j][k][l][m][n][o][p]));
								}
							}


						}
					}
				}
			}
		}
	}
	read(&sav,&(data->numoftloops));
	for (i=0;i<=data->numoftloops;i++) {
		read(&sav,&(data->itloop[i]));
		read(&sav,&(data->tloop[i]));

	}
	read(&sav,&(data->numoftriloops));
	for (i=0;i<=data->numoftriloops;i++) {
		read(&sav,&(data->itriloop[i]));
		read(&sav,&(data->triloop[i]));

	}
	read(&sav,&(data->numofhexaloops));
	for (i=0;i<=data->numofhexaloops;i++) {
		read(&sav,&(data->ihexaloop[i]));
		read(&sav,&(data->hexaloop[i]));

	}
	read(&sav,&(data->auend));
	read(&sav,&(data->gubonus));
	read(&sav,&(data->cint));
	read(&sav,&(data->cslope));
	read(&sav,&(data->c3));
	read(&sav,&(data->efn2a));
	read(&sav,&(data->efn2b));
	read(&sav,&(data->efn2c));
	read(&sav,&(data->init));
	read(&sav,&(data->mlasym));
	read(&sav,&(data->strain));
	read(&sav,&(data->prelog));
	read(&sav,&(data->singlecbulge));
	read(&sav,&(data->maxintloopsize));



	sav.close();
}

//return the pairing probability of the i=j pair, where i<j.
PFPRECISION calculateprobability(int i, int j, pfunctionclass *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce) {
	PFPRECISION interior, exterior;
	int before,after;
	int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},
	{0,1,0,1,0,0},{0,0,0,0,0,0}};
	bool adjacentgu;

	if (!mod[i]&&!mod[j]) {
		if (ct->constant==NULL) return (v->f(i,j)*v->f(j,i+ct->GetSequenceLength()))/(w5[ct->GetSequenceLength()]*scaling*scaling);
		else {
			//constant is being used.
			if (ct->constant[j][i]<EPSILON) return 0.0;
			return (v->f(i,j)*v->f(j,i+ct->GetSequenceLength()))/(w5[ct->GetSequenceLength()]*scaling*scaling*ct->constant[j][i]);
		}
	}
	else {
		if (!(fce->f(i,j)&SINGLE)) {
			before =0;
			if ((i>1&&j<(2*ct->GetSequenceLength())&&j!=ct->GetSequenceLength())) {
				if ((j>ct->GetSequenceLength()&&((i-j+ct->GetSequenceLength())>minloop+2))||j<ct->GetSequenceLength()) {
					before = inc[ct->numseq[i-1]][ct->numseq[j+1]];
				}
			}


			//after = 0 if a stacked pair cannot form 3' to i
			if ((((j-i)>minloop+2)&&(j<=ct->GetSequenceLength())||(j>ct->GetSequenceLength()+1))&&(i!=ct->GetSequenceLength())) {
				after = inc[ct->numseq[i+1]][ct->numseq[j-1]];

			}
			else after = 0;

			adjacentgu = false;
			//check for preceding or following GU or whether the pair itself is gu
			if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) adjacentgu = true;
			else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) adjacentgu = true;
			else if (ct->numseq[i]==3&&ct->numseq[j]==4) adjacentgu = true;
			else if (ct->numseq[i]==4&&ct->numseq[j]==3) adjacentgu = true;
			else if (i-1>0&&j+1<=ct->GetSequenceLength()) {
				if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) adjacentgu = true;
				else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) adjacentgu = true;

			}

			//if there are no stackable pairs to i.j then don't allow a pair i,j
			if ((before!=0)||(after!=0)) {

				if (i+1<j-1&&!adjacentgu) interior = erg1(i,j,i+1,j-1,ct,data)*v->f(i+1,j-1);
				else interior = (PFPRECISION) 0;
				if (j+1<=ct->GetSequenceLength()&&!adjacentgu) exterior = erg1(j,i+ct->GetSequenceLength(),j+1,i+ct->GetSequenceLength()-1,ct,data)*v->f(j+1,i+ct->GetSequenceLength()-1);
				else exterior = (PFPRECISION) 0;
				return ((v->f(i,j)+interior)*(v->f(j,i+ct->GetSequenceLength())+exterior)-interior*exterior)/(w5[ct->GetSequenceLength()]*scaling*scaling);
			}
			else return 0;


		}
		else return 0;

	}

}

//function to rescale all arrays when partition function calculation is headed out of bounds


void rescale(int i, int j, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
			 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, PFPRECISION **wca, 
			 PFPRECISION **curE, PFPRECISION **prevE, PFPRECISION rescalefactor) {



	int h,d,dp,ii,jj,hh,nucs,index,lowi,highi,number;
	double multiplier;
	number=ct->GetSequenceLength();
	d=j-i;
	h=(j<=number)?d:(d+number-1);

#ifdef pfdebugmode
	ofstream ufout;
	ufout.open("c:/underflow_alert.out",ios::app);
	ufout<<"rescale factor = "<<rescalefactor<<"\n";
	ufout.close();
#endif
	//rescale v,w,wl,wcoax,wmb,wmbl,wca
	for (hh=0;hh<=h;hh++){
		//determine the boundary of ii
		if (h < number) {
			if (hh == h) {
				highi=i;
				lowi=1;
			}
			else {
				highi=number-hh;
				lowi=1;
			}
		}
		else {
			if (hh==h) {
				highi=i;
				lowi= 2*number-hh;
			}
			else if(hh >= number) {
				highi=number;
				lowi=2*number-hh;
			}
			else if(hh<number) {
				highi=number-hh;
				lowi=1;
			}
		}
		for (ii=lowi;ii<=highi;ii++){
			jj=(hh<=(number-1))?hh+ii:(hh+ii-number+1);
			nucs = jj-ii+1; //this work even jj> number
			multiplier = pow(rescalefactor,(PFPRECISION) nucs);


#ifdef pfdebugmode
			//look for positions that will underflow

			if (multiplier<0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)<log10(PFMIN)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}



			}

			if (multiplier>0) {
				//values will get smaller...
				if (v->f(ii,jj)>0) {
					if (log10(v->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" v= "<<v->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (w->f(ii,jj)>0) {
					if (log10(w->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" w= "<<w->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wl->f(ii,jj)>0) {
					if (log10(wl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wl= "<<wl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wcoax->f(ii,jj)>0) {
					if (log10(wcoax->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wcoax= "<<wcoax->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmb->f(ii,jj)>0) {
					if (log10(wmb->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmb= "<<wmb->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}
				if (wmbl->f(ii,jj)>0) {
					if (log10(wmbl->f(ii,jj))+log10(multiplier)>log10(PFMAX)) {
					ufout.open("c:/underflow_alert.out",ios::app);
					ufout<<"of i= "<<ii<<" j= "<<jj<<" wmbl= "<<wmbl->f(ii,jj)<<"multiplier = "<<multiplier<<"\n";
					ufout.close();
					}
				}



			}

#endif //debug
			//rescale v,w,wl,wcoax,wmb,wmbl,wca
			v->f(ii,jj)=v->f(ii,jj)*multiplier;
			w->f(ii,jj)=w->f(ii,jj)*multiplier;
			wl->f(ii,jj)=wl->f(ii,jj)*multiplier;
			wcoax->f(ii,jj)=wcoax->f(ii,jj)*multiplier;
			wmb->f(ii,jj)=wmb->f(ii,jj)*multiplier;
			wmbl->f(ii,jj)=wmbl->f(ii,jj)*multiplier;
			if (jj<=number) wca[ii][jj]=wca[ii][jj]*multiplier;
			if (ii==1 && jj<=number) {
				//rescale w5
				w5[jj]=w5[jj]*pow(rescalefactor,(PFPRECISION) jj);

				if (jj==number) {
					//rescale w3
					for (index=1;index<=number;index++) w3[index]=w3[index]*pow(rescalefactor,(PFPRECISION) (number-index+1));

				}
			}
		}
	}

	if (curE!=NULL) {
			//curE and preE are not used in SMP mode:
		//rescale curE and prevE
		for (ii=((h<=(number-2))?1:(2*number-h-1));ii<=((h<=(number-2))?(number-h):number);ii++)
		for (dp=1;dp<=d-1;dp++){
			if (ii<number) {
				curE[dp][ii] = curE[dp][ii]*pow(rescalefactor,(PFPRECISION)(dp+1));
				prevE[dp][ii+1] = prevE[dp][ii+1]*pow(rescalefactor,(PFPRECISION)(dp+1));
			}
		}
	}

	//rescale datatable
	data->rescaledatatable(rescalefactor);
}


//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w3
//void rescaleatw3(int ii, structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {



//}

//function to rescale all arrays when partition function calculation is headed out
	//of bounds when claculating w5
//void rescaleatw5(int jj,structure *ct, pfdatatable *data, pfunctionclass *v, pfunctionclass *w, pfunctionclass *wl, pfunctionclass *wcoax,
//				 pfunctionclass *wmb,pfunctionclass *wmbl, PFPRECISION *w5, PFPRECISION *w3, double rescalefactor) {

	//rescale the previously filled arrays
//	rescale(1, jj, ct, data, v, w, wl, wcoax, wmb, wmbl, w5, w3, rescalefactor);

//	w5[jj]=w5[jj]*pow(rescalefactor,(double) jj);

//}

//rescale the entries in datatable
void pfdatatable::rescaledatatable(PFPRECISION rescalefactor) {

	scaling=scaling*rescalefactor;
	int i,j,k,l,m,n,o,p;


	for (i=0;i<31;i++) {
		inter[i] = inter[i]*pow(rescalefactor,i+2);
		bulge[i] = bulge[i]*pow(rescalefactor,i+2);
		hairpin[i] = hairpin[i]*pow(rescalefactor,i+2);

	}
	for (i=0;i<6;i++) {
		for (j=0;j<6;j++) {
			for (k=0;k<6;k++) {
				for (l=0;l<3;l++) {
					dangle[i][j][k][l] = dangle[i][j][k][l]*rescalefactor;
				}
				for (l=0;l<6;l++) {
					stack[i][j][k][l]=stack[i][j][k][l]*pow(rescalefactor,2);

					tstackcoax[i][j][k][l]=tstackcoax[i][j][k][l]*pow(rescalefactor,2);

					tstack[i][j][k][l]=tstack[i][j][k][l]*pow(rescalefactor,2);
					tstkm[i][j][k][l]=tstkm[i][j][k][l]*pow(rescalefactor,2);

					for (m=0;m<6;m++) {
						for (n=0;n<6;n++) {
							iloop11[i][j][k][l][m][n]=iloop11[i][j][k][l][m][n]*pow(rescalefactor,4);
							for (o=0;o<6;o++) {
								iloop21[i][j][k][l][m][n][o]=iloop21[i][j][k][l][m][n][o]*pow(rescalefactor,5);
								for (p=0;p<6;p++) {
									iloop22[i][j][k][l][m][n][o][p]=iloop22[i][j][k][l][m][n][o][p]*pow(rescalefactor,6);
								}
							}


						}
					}
				}
			}
		}
	}

	for (i=0;i<=numoftloops;i++) {

		tloop[i] = tloop[i]*pow(rescalefactor,6);

	}

	for (i=0;i<=numoftriloops;i++) {

		triloop[i] = triloop[i]*pow(rescalefactor,5);
	}

	for (i=0;i<=numofhexaloops;i++) {

		hexaloop[i] = hexaloop[i]*pow(rescalefactor,8);

	}



}

//thresh-structure builds a structure containing all base pairs above the probability thresh (expressed as a fraction from 0 to 1).
//Note that thresh must be 0.5 or larger for the resulting structure to be a valid secondary structure.
void thresh_structure(structure *ct, char *pfsfile, double thresh) {


	int i,j;
	int vers;
	PFPRECISION *w5, *w3, scaling;
	pfunctionclass *v, *w,*wmb,*wmbl,*wl,*wcoax;
	forceclass *fce;
	bool *mod,*lfce;
	pfdatatable *data;

	//allocate the ct file by reading the save file:
	ifstream sav(pfsfile,ios::binary);


	read(&sav,&(vers));//read the version of the save file
		//right now there is no infrastructure to indicate the wrong version is being read.
		//This should be changed in the future...

	int SequenceLength;
	read(&sav,&(SequenceLength));
	//ct->numofbases = SequenceLength;

	sav.close();

	//allocate everything

	//array = new PFPRECISION *[ct.GetSequenceLength()+1];
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	array[i] = new PFPRECISION [i+1];
	//}
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	for (j=1;j<=i;j++) {
	//		array[i][j]=0;
	//	}
	//}

	ct->allocate(SequenceLength);

	w = new pfunctionclass(ct->GetSequenceLength());
	v = new pfunctionclass(ct->GetSequenceLength());
	wmb = new pfunctionclass(ct->GetSequenceLength());
	wmbl = new pfunctionclass(ct->GetSequenceLength());
	wcoax = new pfunctionclass(ct->GetSequenceLength());
	wl = new pfunctionclass(ct->GetSequenceLength());
	fce = new forceclass(ct->GetSequenceLength());

	w5 = new PFPRECISION [ct->GetSequenceLength()+1];
	w3 = new PFPRECISION [ct->GetSequenceLength()+2];

	lfce = new bool [2*ct->GetSequenceLength()+1];
    mod = new bool [2*ct->GetSequenceLength()+1];

	data = new pfdatatable();

	//load all the data from the pfsavefile:
	readpfsave(pfsfile, ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&scaling,mod,lfce,data);

	//reset the base pairing info:
	//for (i=1;i<=ct->GetSequenceLength();i++) ct->basepr[1][i]=0;

	ct->AddStructure();
	//fill array with the values for the plot:
	for (i=1;i<ct->GetSequenceLength();i++) {
		for (j=i+1;j<=ct->GetSequenceLength();j++) {



			if(calculateprobability(i,j,v,w5,ct,data,lfce,mod,scaling,fce)>thresh) {

				ct->SetPair(i,j);


			}



		}
	}

	//now build the structure

	

	delete w;
	delete v;
	delete wmb;
	delete fce;
	delete[] w5;
	delete[] w3;
	//for (i=1;i<=ct.GetSequenceLength();i++) {
	//	delete[] array[i];
	//}
	//delete[] array;
	delete[] lfce;
	delete[] mod;
	delete data;
	delete wmbl;
	delete wl;
	delete wcoax;


}
