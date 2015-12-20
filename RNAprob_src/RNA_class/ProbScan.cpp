#include "../src/pfunction.h"
#include "ProbScan.h"
#include "../src/structure.h"
#include <algorithm>
#include <functional>
#include <numeric>
#include <iomanip>
#include <assert.h>
using std::vector;

ProbScan::ProbScan(const char filename[], bool from_sequence_file, bool isRNA):RNA(filename,3,true)
{
  if (from_sequence_file)
    PartitionFunction();//calculate the partition function if it hasn't been
                        //done yet
}

PFPRECISION ProbScan::probability_of_individual_hairpin(int i,int j) 
{
  return(v->f(j,i+GetSequenceLength()) //V'(i,j)
         * erg3(i,j,GetStructure(),pfdata,0)) //K for hairpin
         / (w5[GetSequenceLength()]*pfdata->scaling*pfdata->scaling); //Q
//divide by scaling^2 so the closing nucs aren't double counted
}

vector<hairpin_t> ProbScan::probability_of_all_hairpins(int min_size, int max_size,PFPRECISION threshold) 
{
  vector<hairpin_t> hairpins;
  structure* st = GetStructure();
//#pragma omp parallel for
//parallelization causes datarace because of push_back method
  for(int i=1;i<GetSequenceLength()-min_size-1;i++){//search over all 0<i<j<n
    for(int j=i+min_size+1;j<std::max(i+max_size,GetSequenceLength());j++){
      if(inc[st->numseq[i]][st->numseq[j]]){//if i and j can pair
        //get probability
        PFPRECISION probability = probability_of_individual_hairpin(i,j);
        if (probability>threshold){ //add to the list if p>threshold
          hairpins.push_back(hairpin(probability,i,j));
        }
      }
    }
  }
  //hairpins now contains every hairpin where p>threshold
  return hairpins;
}

int ProbScan::count_alternative_bulge_loops(const int i, const int ip)
{//count the number of isoenergetic alternative bulge loops
 //for explanation see Mathews et al. 2004 PNAS 101:19 pp. 7287-92
  int count = 1;
  int k;
  int n = GetSequenceLength();
  int * s = GetStructure()->numseq;//the RNA sequence
  if (i!=n){
    k = i;//first walk towards the 5' end
    while (k>=1 && s[k]==s[i+1]) {
      count++;
      k--;
    }
    k = ip;//then walk towards the 3' end
    while (k<=n && (s[k]==s[i+1])) {
      count++;
      k++;
    }
  }
  return count;
}

bool is_single_bulge(int i,int j,int k,int l)
{
  return ((k-i==1 && j-l==2) || (k-i==2 && j-l==1));
}

PFPRECISION ProbScan::probability_of_internal_loop(int i,int j, int k, int l) 
{
  PFPRECISION probability = (v->f(k,l) //V(k,l)
         * v->f(j,i+GetSequenceLength()) //V'(i,j)
         * erg2(i,j,k,l,GetStructure(),pfdata,0,0))//K for iloop
         / (w5[GetSequenceLength()]*pfdata->scaling*pfdata->scaling);//Q
//we divide by scaling^2 so the closing nucs aren't double counted
 //we perform a correction for slipping by isoenergetic single nucleotide bulges 
 if (is_single_bulge(i,j,k,l)){
    if (k-i==1) probability = probability * (double) count_alternative_bulge_loops(l,j);
    if (j-l==1) probability = probability * (double) count_alternative_bulge_loops(i,k);
  }
  return probability;
}

vector<internal_loop_t> ProbScan::probability_of_all_internal_loops(PFPRECISION threshold) 
{
  vector<internal_loop_t> iloops;//holds internal loops that we find
  int n = GetSequenceLength();
  structure* st = GetStructure();

//#pragma omp parallel for //can't parallelize because of push_back
//search over all i,j,k,l with < max_internal_loop unpaired nucs
  for(int i=1;i<n-3;i++)
    for(int k=i+1;k<std::min(i+max_internal_loop,n-2);k++)
      for(int l=k+minloop+1;l<n-1;l++)
        for(int j=l+1;j<std::min(l+(max_internal_loop-(k-i+1)),n);j++){
          //if i can pair to j and k can pair to l
          if(inc[st->numseq[i]][st->numseq[j]]&&inc[st->numseq[k]][st->numseq[l]]){
            assert(i<k && l<j && k<l);//input validation
            assert(k-i>1 || j-l>1);
            //get probability of the internal loop
            PFPRECISION probability=probability_of_internal_loop(i,j,k,l);
            if (probability>threshold) {//add to list if prob>threshold
              iloops.push_back(internal_loop(probability,i,j,k,l));
            }
          }
        }
  return iloops;//vector now holds all possible iloops with p>threshold
}

//element class represents an element of a multibranch loop
//either a helix closed by i and j or an unpaired nucleotide at i
class element{
 public:
  int i; 
  int j; 
  bool is_a_pair;
  element(branch h) : i(h.first),j(h.second),is_a_pair(true) {}
  element(int nuc) : i(nuc),j(0),is_a_pair(false) {}
};

vector<element> ProbScan::construct_element_array(const multibranch_loop_t& mb)
{
//construct the element array containing hairpins and nucleotides
//first make the closing hairpin,
//swapping its indices so its the same orientation as the others
  vector<element> element_array;
  element closing_hairpin = element(std::make_pair(mb.branches[0].second,mb.branches[0].first));
  element_array.push_back(closing_hairpin);
  bool first=true;
  bool last=false;
//mb is a vector of pairs
//for each pair in the multibranch loop after the first one, add any unpaired 
//nucleotides between the last pair and this pair, then add the new pair
  for(vector<branch>::const_iterator it=mb.branches.begin()+1;it!=mb.branches.end();++it){
    element last_pair = element(*(it-1));//get last pair
    element next_pair = element(*it);//and next pair
    assert(last_pair.is_a_pair && next_pair.is_a_pair);//input validation
    for(int x = first?last_pair.i+1:last_pair.j+1;x<next_pair.i;x++){
      element_array.push_back(element(x));//insert unpaired nucs between
                                          //last_pair and next_pair
    }
    first=false;
    element_array.push_back(next_pair); //add next_pair, which will be the new
                                        //last_pair
  }
//add the last few unpaired nucs
  for(int x = element_array.back().j+1;x<element_array.front().i;x++){
    element_array.push_back(element(x));
  }

  for(int i=0;i<4;i++){//duplicate the first 4 elements of the array on the end
    element_array.push_back(element_array[i]);
  }
  return element_array;
}

PFPRECISION prev_val(int index,int offset, vector<vector<PFPRECISION> >& arr)
{
  const PFPRECISION initial_value = 1.0;
  if(index>=offset) return arr[index][offset];
  else return initial_value;
}

PFPRECISION ProbScan::equilibrium_constant_for_multibranch_loop(const multibranch_loop_t& mb)
{
  vector<element> elements = construct_element_array(mb);
  if(elements.size()<=8) return 0.0;//Rahul did this.. need to think more about justification
//N+4 by 4 array for accumulating the partition function
  vector<vector<PFPRECISION> > arr(elements.size(),vector<PFPRECISION>(4,1.0));
  int* s = GetStructure()->numseq;//the nucleotide sequence
  int offset = 0;
  const int n = elements.size()-4;
  int stems=0,unpaired_nucs=0;
//first calculate the partition function around the circle starting at 4
//different start points "offsets". we'll handle interactions across the //starting points later
  for(int offset=0;offset<4;offset++){
    for(int x=offset;x<n+offset;x++){
      PFPRECISION result = 0.0;
      if(elements[x].is_a_pair){
        if(offset==0) stems += 1;
  //case where new helix has no coaxial stack or danging end
        result += prev_val(x-1,offset,arr);

  //case where the previous nuc is making a 5' dangle on new helix
        if(!elements[x-1].is_a_pair && x>offset){
          result += prev_val(x-2,offset,arr)*erg4(elements[x].j,elements[x].i,elements[x-1].i,2,
                                  GetStructure(),pfdata,0);
        }
#ifndef disablecoax
  //case where new helix is coaxially stacking flush with helix at element x-1
        if(elements[x-1].is_a_pair && x>offset){
          result += prev_val(x-2,offset,arr) * 
                     ergcoaxflushbases(elements[x-1].i,elements[x-1].j,
                                       elements[x].i,elements[x].j,
                                       GetStructure(),pfdata);
        }
  //case where new helix is coaxially stacking on element x-2 with an intervening mismatch
  //coax 1: , || , ||
        if(x>2+offset && elements[x-2].is_a_pair && 
            !elements[x-1].is_a_pair && !elements[x-3].is_a_pair){
          result += prev_val(x-4,offset,arr) *
                     ergcoaxinterbases1(elements[x-2].i,elements[x-2].j,
                                        elements[x].i,elements[x].j,
                                        GetStructure(),pfdata);
        }
#endif //disablecoax 
      }
      else if(!elements[x].is_a_pair){
        if(offset==0) unpaired_nucs += 1;
  //case where nuc is not dangling or participating in a mismatch
  //we have to scale for the unpaired nuc
        result += prev_val(x-1,offset,arr)*pfdata->scaling;
  //case where nuc is dangling 3' on helix at element x-1
        if(x>offset && elements[x-1].is_a_pair){
          result += prev_val(x-2,offset,arr)*erg4(elements[x-1].j,elements[x-1].i,elements[x].i,1,
                                  GetStructure(),pfdata,0);
        }
  //case where element x and x-2 are nucs forming a mismatch on helix at x-1
  // , || ,
        if(x>1+offset && elements[x-1].is_a_pair && !elements[x-2].is_a_pair){
          result += prev_val(x-3,offset,arr) * 
                    pfdata->tstkm[s[elements[x-1].j]][s[elements[x-1].i]]
                               [s[elements[x].i]][s[elements[x-2].i]];

        }
  //case where new nuc is participating in a mismatch coaxial stack between helices at x-1 and x-3
  //coax 2: || , || ,
#ifndef disablecoax
        if(x>2+offset && !elements[x-2].is_a_pair && 
            elements[x-1].is_a_pair && elements[x-3].is_a_pair){
          result += prev_val(x-4,offset,arr) * 
                     ergcoaxinterbases2(elements[x-3].i,elements[x-3].j,
                                        elements[x-1].i,elements[x-1].j,
                                        GetStructure(),pfdata);
        } 
#endif //disablecoax
      }
      arr[x][offset] = result;
    }
  }
//now let's paste the ends of the sequence together
//there are only a few cases because element 0 is always 
//a helix in current implementation
  PFPRECISION initiation = 1 
                  * pow(pfdata->eparam[10],stems) //per stem penalty "c"
                  * pow(pfdata->eparam[6],unpaired_nucs) //per nuc penalty "b"
                  * pfdata->eparam[5];//closing multibranch loop penality "a"
  PFPRECISION pfunc = arr[n-1][0];
#ifndef disablecoax
  if(elements[0].is_a_pair && elements[n-1].is_a_pair){
    pfunc += arr[n-2][1]*ergcoaxflushbases(elements[n-1].i,elements[n-1].j,
                                       elements[0].i,elements[0].j,
                                       GetStructure(),pfdata);
  }
#endif//disablecoax
  if(elements[0].is_a_pair && !elements[n-1].is_a_pair){//5' dangle
    pfunc += arr[n-2][1]*erg4(elements[0].j,elements[0].i,elements[n-1].i,2,
                                  GetStructure(),pfdata,0);
  }
//terminal mismatch
  if(elements[0].is_a_pair && !elements[1].is_a_pair && !elements[n-1].is_a_pair){
    pfunc += arr[n-2][2]*pfdata->tstkm[s[elements[n].j]][s[elements[n].i]]
                               [s[elements[n+1].i]][s[elements[n-1].i]];
   }
//3 mismatch coax possibilities
// case like this: 5' || , ...  || ,  3'
#ifndef disablecoax
  if(elements[0].is_a_pair && elements[n-2].is_a_pair &&
     !elements[1].is_a_pair && !elements[n-1].is_a_pair){
    pfunc += arr[n-3][2] * ergcoaxinterbases2(elements[n-2].i,elements[n-2].j,
                                        elements[n].i,elements[n].j,
                                        GetStructure(),pfdata);
  }
// case like this: 5' || ... , || ,  3'
  if(!elements[n-3].is_a_pair && !elements[n-1].is_a_pair &&
     elements[n-2].is_a_pair && elements[n].is_a_pair){
    pfunc += arr[n-4][1] * ergcoaxinterbases1(elements[n-2].i,elements[n-2].j,
                                      elements[n].i,elements[n].j,
                                      GetStructure(),pfdata);
    }
// case like this: 5' || , ||   ...   ,  3'
  if(!elements[n-1].is_a_pair && !elements[n+1].is_a_pair &&
     elements[n].is_a_pair && elements[n+2].is_a_pair){
    pfunc += arr[n-2][3] * ergcoaxinterbases1(elements[n].i,elements[n].j,
                                      elements[n+2].i,elements[n+2].j,
                                      GetStructure(),pfdata);
    }
#endif//disablecoax
  return pfunc*initiation;
}

PFPRECISION ProbScan::probability_of_multibranch_loop(const multibranch_loop_t& mb)
{
  assert(mb.size()>=3);
  //holds the values from v array
  vector<PFPRECISION> vs;
  //V(j,i+numberofbases) for closing pair
  vs.push_back(v->f(mb.branches[0].second,mb.branches[0].first+GetSequenceLength())
                   *penalty(mb.branches[0].second,mb.branches[0].first,GetStructure(),pfdata));
  //V(i,j) for each branch, with AU/GU end penalty
  for(vector<branch>::const_iterator it=mb.branches.begin()+1;it!=mb.branches.end();++it){
    vs.push_back(v->f(it->first,it->second)
                   *penalty(it->first,it->second,GetStructure(),pfdata));
  }
  //calculate equilibrium constant
  PFPRECISION Kmb = equilibrium_constant_for_multibranch_loop(mb);
  //take product of values from V array
  PFPRECISION product_of_vs = std::accumulate(vs.begin(),vs.end(),1.0
                                     ,std::multiplies<PFPRECISION>());
  //return probability
  return (Kmb * product_of_vs) / w5[GetSequenceLength()];
}


//functions for dealing with structures
inline hairpin_t hairpin(PFPRECISION p,int i, int j)
{
  hairpin_t h;
  h.probability = p;
  h.i = i;
  h.j = j;
  return h;
}

inline internal_loop_t internal_loop(PFPRECISION p,int i, int j, int k, int l)
{
  internal_loop_t il;
  il.i=i;
  il.j=j;
  il.k=k;
  il.l=l;
  il.probability=p;
  return il;
}

//to build a multibranch loop, call this with the closing base pair
//then call add_branch with the branches, from 5' to 3'
inline multibranch_loop_t multibranch_loop(int i, int j)
{
  multibranch_loop_t mb;
  mb.branches.push_back(std::make_pair(i,j));
  return mb;
}

inline void add_branch(multibranch_loop_t& mb,int k, int l)
{
  mb.branches.push_back(std::make_pair(k,l));
} 

void show_hairpins(vector<hairpin_t> hairpins)//print hairpin output
{
  cout <<"--hairpins--"<<endl;
  cout << "prob i j" <<endl;
  for(vector<hairpin_t>::const_reverse_iterator it=hairpins.rbegin();it!=hairpins.rend();++it)
    cout << std::fixed<<std::setprecision(3)<<it->probability << " " << it->i << " " << it->j <<endl; 
  cout<< "--hairpins end--"<<endl <<endl;
}

void show_internal_loops(vector<internal_loop_t> internals)//print iloop output
{
  cout << "--internal loops--"<<endl;
  cout << "prob i j k l"<<endl;
  for(vector<internal_loop_t>::const_reverse_iterator it=internals.rbegin();it!=internals.rend();++it)
    cout << std::fixed<<std::setprecision(3)<< it->probability << " " << it->i << " " << it->j <<" " << it->k << " " << it->l <<endl;
  cout<< "--internal loops end--"<<endl <<endl;
}

void show_element_array(vector<element> e)//just for debugging
{
  int count = 0;
  for(vector<element>::iterator it=e.begin();it!=e.end();++it){
    cout<<count++<<" ";
    cout << (it->is_a_pair?"Pair: ":"Nuc ") << it->i<<" ";
    if(it->is_a_pair) cout <<it->j;
    cout <<std::endl;
  }
}

void show_mbl(multibranch_loop_t mb){//print multibranch loop output
  cout <<mb.probability;
  for(vector<branch>::iterator it = mb.branches.begin();it!=mb.branches.end();++it){
    cout<<'\t'<<it->first<<'-'<<it->second;
  }
  cout<<'\n';
}
