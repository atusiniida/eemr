//
// GeneSetAnalysis.h
//
#ifndef GeneSetAnalysis_H
#define GeneSetAnalysis_H

#include <map>
#include <string>
#include <vector>
using namespace std;

namespace eem{

class GeneSetAnalysis{
public:
	GeneSetAnalysis();
	virtual ~GeneSetAnalysis();

	virtual void perform(){ ;}//=0;
	virtual void setGeneSets(map<string, vector<string> > geneSets){ ;}//=0;
	virtual void setMaxGeneSetSize(int i){ ;}//=0;
	virtual void setMinGeneSetSize(int i){ ;}//=0;
};
} // namespace eem
#endif // GeneSetAnalysis_H
