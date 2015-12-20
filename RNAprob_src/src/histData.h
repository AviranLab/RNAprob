#ifndef HISTDATA_H
#define HISTDATA_H

#include <vector>
using namespace std;

class histData 
{
private:
	int size;
	double start;
	double binSize;
	vector<double> prob;
public:
	histData();
	histData(double startPos, double bin_size);
	//~histData();
	int getSize();
	double getStart();
	double getEnd();
	double getBinSize();
	void setStart(double startPos);
	void setBinSize(double bin_size);
	double getProbability(int index);
	double getProbability(double reactivity);
	void add(double probability);
	void print();
	
};



#endif