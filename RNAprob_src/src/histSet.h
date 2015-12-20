#ifndef HISTSET_H
#define HISTSET_H

#include <vector>
#include "histData.h"
#include "defines.h"

using namespace std;

class histSet 
{
private:
	//vector<histData *>* shapeHist;
	vector<vector<histData*>*>* shapeHist;
	vector<vector<histData*>*>* dmsHist;
	
	void init();
	void getDataTitle(string title, int& dataSource, int& strucType, int& baseType);
	void setHistParam(int dataSource, int strucType, int baseType, double startPos, double binSize);
	void add(int dataSource, int strucType, int baseType, double probability);
	
public:
	histSet();
	~histSet();	
	void readHistFile(const char* filename);
	histData* getHistData(int dataSource, int strucType, int baseType);
	void print();
};

int extractDataSource(string ds);
int extractStructureType(char* stype);
int extractBaseType(char* base);
bool isEmptyLine(string str);




#endif