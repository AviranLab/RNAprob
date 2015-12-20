#ifndef PROBSCAN_CLASS
#define PROBSCAN_CLASS

#include "RNA.h"
#include <vector>

const static int max_internal_loop=30;
const static int inc[6][6]={{0,0,0,0,0,0},{0,0,0,0,1,0},{0,0,0,1,0,0},{0,0,1,0,1,0},{0,1,0,1,0,0},{0,0,0,0,0,0}};//array representing legal base pairs

//types which contain nucleotide indices which specify a structure
//a hairpin closed by i and j
typedef struct hp {PFPRECISION probability;
                   int i;
                   int j;} hairpin_t;
//an internal loop closed by i and j on the exterior and k and l in the interior, where i<k<l<j
typedef struct il {PFPRECISION probability;
                   int i;
                   int j;
                   int k;
                   int l;} internal_loop_t;
//a branch of a multibranch loop
typedef std::pair<int,int> branch;
//a multibranch loop, represented by a vector of branches.
//branches[0] must close the multibranch loop
typedef struct mb{PFPRECISION probability;
                  std::vector<branch> branches;} multibranch_loop_t;

//functions to make loops
//return a loop with probability and closing nuc indices provided in arguments
inline hairpin_t hairpin(PFPRECISION p,int i, int j);
inline internal_loop_t internal_loop(PFPRECISION p,int i, int j, int k, int l);
//to build a multibranch loop, call multibranch_loop with the closing base pair
//then call add_branch with the branches, from 5' to 3'
inline multibranch_loop_t multibranch_loop(int i, int j);
inline void add_branch(multibranch_loop_t& mb,int k, int l);

//display loops and their probabilities to stdout
void show_hairpins(vector<hairpin_t>);
void show_internal_loops(vector<internal_loop_t>);
void show_mbl(multibranch_loop_t mbl);

class element;//a pair or an unpaired nucleotide,used in multibranch loop probability calculation

class ProbScan : public RNA
{
 public:
  ProbScan(const char[],bool,bool);//constructor takes partition function save file 
//OR sequence file in which case it will calculate partition function

//return probability of a hairpin closed by (i,j)
  PFPRECISION probability_of_individual_hairpin(int i,int j);
//search over all possible hairpins, return a vector of hairpin_t
//for all hairpins with probability>threshold
  std::vector<hairpin_t> probability_of_all_hairpins(
                           int min,int max,PFPRECISION threshold);
//return probability of a iloop closed by (i,j)
  PFPRECISION probability_of_internal_loop(int,int,int,int);
//search over all possible iloops, return a vector of hairpin_t
//for all iloops with probability>threshold
  std::vector<internal_loop_t> probability_of_all_internal_loops(PFPRECISION);

//return probability of a multibranch loop defined by a multibranch_loop_t
  PFPRECISION probability_of_multibranch_loop(const multibranch_loop_t&);
 private:
//calculate equilibrium constant for a multibranch loop defined by 
//a multibranch_loop_t for use in probability calculation
  PFPRECISION equilibrium_constant_for_multibranch_loop(const multibranch_loop_t&);
//helper functions for Kmb calculation
  std::vector<element> construct_element_array(const multibranch_loop_t&);
  int count_alternative_bulge_loops(const int,int);
};
//print element array for debugging multibranch calculation
void show_element_array(vector<element>);

#endif
