//
// AbstractEEMsearch.h
//
#ifndef AbstractEEMsearch_H
#define AbstractEEMsearch_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
using namespace std;

#include "../utility/MyMat.h"
#include "../utility/StopWatch.h"
using namespace utility;

#include "AbstractGeneSetAnalysis.h"
#include "EEMsearch.h"
#include "EEM.h"
#include "ExpressionModule.h"
#include "ExpressionModuleSet.h"

namespace eem{


///
/// inheritance relationships are:
/// CoherenceBasedEEMsearch->AbstractEEMsearch->AbstractGenesetAnalysis
///

class AbstractEEMsearch : public AbstractGeneSetAnalysis, public EEMsearch{
public:
	AbstractEEMsearch();
	AbstractEEMsearch(AbstractEEMsearch &A);
	virtual ~AbstractEEMsearch();

	//--
	//
	//--
	MyMat  originalExp;
	int    itrForPvalue2Calculation; // iteration number for calculation precise p-values (Pvalue2)
	double Pvalue1Cutoff; // cutoff for approximate p-values (Pvalue1), used to find candidates, which are subjected to Pvalue2 calculation
	StopWatch stopWatch;

	map<string, EEM*> eems;	// store EEM objects for each geneset. keys are seed  gene set ID
	vector<string> seeds; //  seed (input) gene set IDs
	vector<string> candidates; // IDs of gene set subjected to Pvalue2 calculation
	string errLog; 
	string timeLog;
	map<int, vector<int> > nullDistrbutionForRecycle;


	virtual void setGeneSets(map<string, vector<string> > geneSets);
	virtual void setPvalue1Cutoff(double d);
	virtual void suppressPvalue1Cutoff();
	virtual void setItrForPvalue2Calculation(int i);
	virtual void recycleNullDistribution();

protected:
	// these are used for EEM search and to be parallelized
	virtual void findModuleGenes();
	virtual void calculatePvalue1();
	virtual void findCandidates(); // find candidates, which are subjected to Pvalue2 calculation, based on Pvalue1Cutoff
	virtual void calculatePvalue2();

public:
	virtual void printResults(string outfile) throw ();
	virtual void printLog(string outfile) throw ();

protected:
	virtual void setEEM();

public:
	//
	// generate ExpressionModuleSet, ExpressionModule
	//
	virtual ExpressionModuleSet* getExpressionModuleSet(); // get result of EEM search
	virtual string getLog();

protected:
	virtual void writeTimeLog();

public:
	virtual void perform();
	virtual map<string, double> getPvalues();

	// add method
	string name(){ return "AbstractEEMsearch";}

};
} // namespace eem
#endif // AbstractEEMsearch_H
