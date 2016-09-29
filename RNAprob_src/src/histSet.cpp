#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <algorithm>
#include <string.h>
#include "histSet.h"
#include "defines.h"

using namespace std;

histSet::histSet()
{
	init();
}

/*
	initialize vectors for storing the histogram data
*/
void histSet::init()
{	
	//allocate space for storing shape histogram
	shapeHist = new vector<vector<histData*>*>();
	for (int i = 0; i < NUM_BASE; i++){
		vector<histData*>* temp = new vector<histData*>;
		shapeHist->push_back(temp);
		for (int j = 0; j < NUM_STYPE; j++){
			histData* ahist = new histData();
			temp->push_back(ahist);
		}
	}	
	
	//allocate space for storing dms histogram
	dmsHist = new vector<vector<histData*>*>();
	for (int i = 0; i < NUM_BASE; i++){
		vector<histData*>* temp = new vector<histData*>;
		dmsHist->push_back(temp);
		for (int j = 0; j < NUM_STYPE; j++){
			histData* ahist = new histData();
			temp->push_back(ahist);
		}
	}	
}

/*
	destructor, to be improved
*/
histSet::~histSet()
{
	if (shapeHist != NULL){
		//shapeHist->clear();
		delete[] shapeHist;
	}
	if (dmsHist != NULL){
		//dmsHist->clear();
		delete[] dmsHist;
	}
}

/*
	Test if a string is empty (cotains space only)
*/
bool isEmptyLine(string str)
{
	string newstr(str);
	std::string::iterator end_pos = std::remove(newstr.begin(), newstr.end(), ' ');
	newstr.erase(end_pos, newstr.end());
	if (newstr.length() == 0)
		return true;
	return false;
}

/*
	read histogram data file (train_param) 
*/
void histSet::readHistFile(const char* filename)
{
	int baseType = -1;
	int strucType = -1;
	int dataSource = -1;
	double startPos, binSize, probability;
	string str;
	
	ifstream in(filename);
	while (!in.eof()){
		in>>str;
		if ( isEmptyLine(str) )
			continue;
		if (str.find(">") != string::npos){
			//a new histogram data begins, extract info. from header line
			getDataTitle(str, dataSource, strucType, baseType);
			//get starting reactivity and bin size of the histogram
			in>>startPos;
			in>>binSize;
			setHistParam(dataSource, strucType, baseType, startPos, binSize);
		}
		else {
			probability = atof(str.c_str());
			add(dataSource, strucType, baseType, probability);
		}
	}
	
	in.close();
}

/*
	given a header line, get structure type and base type info.
*/
void histSet::getDataTitle(string title, int& dataSource, int& strucType, int& baseType)
{
    baseType = -1;
    strucType = -1;
	
	//get data source
	dataSource = extractDataSource(title);
	
    char* token;
    char* cstr = new char[title.length() + 1];
    strcpy(cstr, title.c_str());
    token = strtok(cstr, " |");

    //get structure type
    token = strtok(NULL, " |");
	strucType = extractStructureType(token);
	
    //get base type
    token = strtok(NULL, " |");
	if (token != NULL)
		baseType = extractBaseType(token);

    delete[] cstr;
}

/*
	set the starting value and bin size of a histogram
*/
void histSet::setHistParam(int dataSource, int strucType, int baseType, double startPos, double binSize)
{
	histData* ahist = getHistData(dataSource, strucType, baseType);
	ahist->setStart(startPos);
	ahist->setBinSize(binSize);
}

/*
	get a particular histogram data 
*/
histData* histSet::getHistData(int dataSource, int strucType, int baseType)
{
	if (dataSource == DSOURCE_SHAPE)
		return shapeHist->at(baseType)->at(strucType);
	else if(dataSource == DSOURCE_DMS)
		return dmsHist->at(baseType)->at(strucType);
	else{
		cout<<"Error: unrecognized data source"<<endl;
		return NULL;
	}
}

/*
	add a (reactivity, frequency) data to the corresponding histogram
*/
void histSet::add(int dataSource, int strucType, int baseType, double probability)
{
	histData* ahist = getHistData(dataSource, strucType, baseType);
	if (ahist != NULL)
		ahist->add(probability);
}


/*
	test function. print the histogram data 
*/
void histSet::print()
{
	cout<<"printing SHAPE histogram data"<<endl;
	for (int i=0; i < NUM_BASE; i++){
		for (int j=0; j < NUM_STYPE; j++){
			histData* ahist = getHistData(DSOURCE_SHAPE, j, i);
			
			string stype="";
			if (i==0) stype = "A";
			else if (i==1) stype = 'C';
			else if (i==2) stype = 'G';
			else if (i==3) stype = 'U';
			else if (i==4) stype = 'X'; 
			string btype = "";
			if (j==0) btype = "All";
			else if(j==1) btype= "unpaired";
			else if(j==2) btype = "paired";
			else if (j==3) btype="helix";
			else if(j==4) btype="stacked";
				
			cout<<"\n>SHAPE|"<<btype<<"|"<<stype<<endl;
			ahist->print();
				
		}
	}
	
	cout<<"\nprinting DMS histogram data"<<endl;
	for (int i=0; i < NUM_BASE; i++){
		for (int j=0; j < NUM_STYPE; j++){
			histData* ahist = getHistData(DSOURCE_DMS, j, i);
			string stype="";
			if (i==0) stype = "A";
			else if (i==1) stype = 'C';
			else if (i==2) stype = 'G';
			else if (i==3) stype = 'U';
			else if (i==4) stype = 'X'; 
			string btype = "";
			if (j==0) btype = "All";
			else if(j==1) btype= "unpaired";
			else if(j==2) btype = "paired";
			else if (j==3) btype="helix";
			else if(j==4) btype="stacked";
				
			cout<<"\n>DMS|"<<btype<<"|"<<stype<<endl;
			ahist->print();
		}
	}
}

/*
	convert structure type into corresponding int value
*/
int extractStructureType(char* stype)
{
	if (strcmp(stype, "all") == 0)
		return STYPE_ALL;
	else if (strcmp(stype, "unpaired") == 0)
		return STYPE_UNPAIRED;
	else if (strcmp(stype, "paired") == 0)
		return STYPE_PAIRED;
	else if (strcmp(stype, "helix_end") == 0)
		return STYPE_HELIXEND;
	else if (strcmp(stype, "stacked") == 0)
		return STYPE_STACKED;
	else {
		cout<<"Error: unidentified structure type!"<<endl;
		return -1;
	}
}

/*
	extract data source from header line, a header line is sth
	like >SHAPE|unpaired.
*/
int extractDataSource(string ds)
{
	if (ds.find("SHAPE") != string::npos)
		return DSOURCE_SHAPE;
	else if (ds.find("DMS") != string::npos)
		return DSOURCE_DMS;
	else {
		cout<<"Error: unidentified data source!"<<endl;
		return -1;
	}
}

/*
	convert a base to the corresponding int value
*/
int extractBaseType(char* base)
{
    if (strcmp(base, "A") == 0)
		return BASE_A;
    else if (strcmp(base, "C") == 0)
		return BASE_C;
    else if (strcmp(base, "G") == 0)
		return BASE_G;
    else if (strcmp(base, "U") == 0)
		return BASE_U;
    else if (strcmp(base, "X") == 0)
		return BASE_X;
	else {
		cout<<"Error: unidentified base type!"<<base<<endl;
		return -1;
	}
}
















