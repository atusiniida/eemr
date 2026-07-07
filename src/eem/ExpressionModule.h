//
// ExpressionModule.h
//
#ifndef ExpressionModule_H
#define ExpressionModule_H

#include <map>
#include <vector>
#include <string>
#include <cfloat>
#include <cmath>
using namespace std;

#include "../utility/MyFunc.h"
#include "../utility/MyMat.h"
using namespace utility;

#include "EEM.h"

//
// a ExpressionModule class object contains EEM search result for each gene set
//
namespace eem{
class ExpressionModuleSet;
class ExpressionModule{
private:
	ExpressionModule(){};
public:
	ExpressionModule(ExpressionModule& e);
	ExpressionModule(const string& id, EEM *eem, MyMat &Exp);
	virtual ~ExpressionModule();

private:

	string id;
	EEM  *eem;
	map<string, double> activityProfile; // activity module of expression modules
	vector<string>      seedGenes;
	vector<string>      moduleGenes;

	// for multiple hypothesis testing (used when EEM search is performed with multiple parameters)
	double minimunPvalue; //-log10 scale
	double PvalueCorrectedForMultipleTest;
	int    numberOfMultipleTesting;

	ExpressionModuleSet* expressionModuleCluster;// for ExpressionModuleClustering


public:
	string getId();
	double getPvalue();
	EEM*   getEEM();

	string toString();

	///
	/// obtain better ExpressionModules based on p-values while correcting multiple hypothesis testing
	///
	static ExpressionModule* getBetterExpressionModules(ExpressionModule& e1, ExpressionModule& e2);
};

ostream &operator <<(ostream &stream, ExpressionModule &em);

} // namespace eem

#endif // ExpressionModule_H
