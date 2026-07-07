//
// AbstractGeneSetAnalysis.h
//
#ifndef AbstractGeneSetAnalysis_H
#define AbstractGeneSetAnalysis_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
using namespace std;

#include "../utility/MyFunc.h"
using namespace utility;

#include "GeneSetAnalysis.h"

namespace eem{

class AbstractGeneSetAnalysis : public GeneSetAnalysis{
public:
	AbstractGeneSetAnalysis();
	virtual ~AbstractGeneSetAnalysis();

	vector<string> allGenes;

	//for filtering of input gene sets based on their size
	int maxSeedGeneSize;// = 2000;
	int minSeedGeneSize;// = 10 ;
	map<string, vector<string> > seedGeneSets;

	//--
	//
	//--
	virtual void setGeneSets(map<string, vector<string> > geneSets);
	virtual void setMaxGeneSetSize(int i);
	virtual void setMinGeneSetSize(int i);

	virtual void perform();

};
} // namespace eem
#endif // AbstractGeneSetAnalysis_H
