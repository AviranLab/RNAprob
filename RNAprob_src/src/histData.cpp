/**
 * structure for a historgram data set
 */
#include <iostream>
#include "histData.h"

using namespace std;

histData::histData()
{
	size = 0;
	start = 0.0;
	binSize = 0.1;
	//allocate space for vectors
	//prob = new vector<double>();
}

histData::histData(double startPos, double bin_size)
{
	size = 0;
	start = startPos;
	binSize = bin_size;
	//allocate space for vectors
	//prob = new vector<double>();
}

/*
	get the probability for a given bin
*/
double histData::getProbability(int index)
{
	if (size == 0)
		return 0;
	//assign outlier to boundry
	if (index < 0)
		return prob.at(0);
	if (index >= size)
		return prob.at(size-1);
	return prob.at(index);
}

/*
	get the probability for a given reactivity value
*/
double histData::getProbability(double reactivity)
{
	if (size == 0)
		return 0;
	double delta = 0.000001;
	int index = (reactivity - start + delta) / binSize;
	//assign outlier to boundry
	if (index < 0)
		return prob.at(0);
	if (index >= size)
		return prob.at(size-1);
	return prob.at(index);
}

/*
	add a new data to the histogram
*/
void histData::add(double probability)
{
	prob.push_back(probability);
	size++;
}

/*
	get the starting reactivity (lower boundary) of the histogram
*/
double histData::getStart()
{
	return start;
}

void histData::setStart(double startPos)
{
	start = startPos;
}

/*
	get the last reactivity (upper boundary) of the histogram 
*/
double histData::getEnd()
{
	return start + binSize * size;
}

/*
	get the bin size of the histogram
*/
double histData::getBinSize()
{
	return binSize;
}


void histData::setBinSize(double bin_size)
{
	binSize = bin_size;
}

/*
	get the size of a histogram
*/
int histData::getSize()
{
	return size;
}

/*
	print the histogram details (For test purpose)
*/
void histData::print()
{
	if (getSize() <=0 )
		return;
	cout<<"Start printing histogram..."<<endl;
	cout<<"Reactivity starts at:"<<getStart()<<"  ends at:"<<getEnd()<<endl;
	cout<<"#bins:"<<getSize()<<"  bin size:"<<getBinSize()<<endl;
	cout<<"Probability for each bin:\n";
	for (int i = 0; i < getSize(); i++)
		cout<<getProbability(i)<<endl;
	cout<<"End of histogram\n";
}

/*
	delete the structure
*/
/*histData::~histData()
{
	//delete spaces allocated to vectors
	prob->clear();
	delete[] prob;
}*/