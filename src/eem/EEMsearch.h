//
// EEMsearch.h
//
#ifndef EEMsearch_H
#define EEMsearch_H

#include <exception>
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "../utility/MyException.h"
using namespace utility;

#include "ExpressionModuleSet.h"

namespace eem{

class EEMsearch{
public:
	EEMsearch();
	virtual ~EEMsearch();

	virtual void setPvalue1Cutoff(double d){ ;}//=0;
	virtual void suppressPvalue1Cutoff(){ ;}//=0;

	////////////////virtual void perform(){ ;}//=0;

	virtual ExpressionModuleSet* getExpressionModuleSet(){ ExpressionModuleSet* expModSet=NULL; return expModSet;}//=0;
	virtual void printResults(string outfile) throw (){ ;}
	virtual void printLog(string outfile) throw (){ ;}
	virtual map<string, double> getPvalues(){ map<string,double> pval;  return pval;}//=0;

	//////////////virtual void setGeneSets(map<string, vector<string> > geneSets){ ;}//=0;
	//////////////virtual void setMaxGeneSetSize(int i){ ;}//=0;
	//////////////virtual void setMinGeneSetSize(int i){ ;}//=0;

	virtual void setItrForPvalue2Calculation(int i){ ;}//=0;
	virtual void recycleNullDistribution(){ ;}//=0;
	virtual string getLog(){ string str; return str;}//=0;
};
} // namespace eem
#endif // EEMsearch_H
