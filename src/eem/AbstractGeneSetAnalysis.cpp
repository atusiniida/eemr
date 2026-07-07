//
// AbstractGeneSetAnalysis.cpp
//
#include "AbstractGeneSetAnalysis.h"
using namespace eem;

AbstractGeneSetAnalysis::AbstractGeneSetAnalysis()
{
	maxSeedGeneSize = 2000;
	minSeedGeneSize = 10 ;
}
AbstractGeneSetAnalysis::~AbstractGeneSetAnalysis()
{
}
//--
//
//--
void AbstractGeneSetAnalysis::setGeneSets( map<string, vector<string> > geneSets )
{
	map<string, vector<string> >::iterator it;

	for(it=geneSets.begin(); it != geneSets.end(); it++){
		vector<string> tmp = MyFunc::isect( it->second, allGenes );

		if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
			std::cerr << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!" << endl;
			continue;
		}
		//
		// map => multimap ?
		//
		seedGeneSets[it->first] = tmp;// Java: seedGeneSets.put(it->first, tmp);
	}
}
void AbstractGeneSetAnalysis::setMaxGeneSetSize(int i)
{
	maxSeedGeneSize = i;
}
void AbstractGeneSetAnalysis::setMinGeneSetSize(int i)
{
	minSeedGeneSize = i;
}
void AbstractGeneSetAnalysis::perform()
{
	cerr << "Unsupported Operation  Exception : AbstractGeneSetAnalysis::perform()" << endl;
}
