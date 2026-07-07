//
// CoherenceBasedEEM.h
//
#ifndef CoherenceBasedEEM_H
#define CoherenceBasedEEM_H

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cmath>
using namespace std;

#include "../utility/MyException.h"
#include "../utility/MyFunc.h"
using namespace utility;

#include "EEM.h"

namespace eem{
//
//public class CoherenceBasedEEM implements EEM, Serializable{
//
class CoherenceBasedEEMsearch;

///
///an  CoherenceBasedEEM  class object performs an EEM seach for each geneset
///
///
    
class CoherenceBasedEEM : public EEM {
//protected:
public:
	CoherenceBasedEEM();
	CoherenceBasedEEM(CoherenceBasedEEMsearch *parent, vector<string>& gene);
	virtual ~CoherenceBasedEEM();

protected:

	double absoluteRadius;
	double relativeRadius;
	vector<string> seedGenes;
	vector<string> moduleGenes;

	double Pvalue1; // approximate p-value (minus log scale)
	double Pvalue2; // precise p-value (minus log scale)
	double Pvalue;  // Pvalue2 if available, or Pvalue1

	int itrForPvalue2Calculation;
	int coreGeneSize;
	vector<double> center;

	//--
	// transient(Java)
	//--
	vector<int> nullDist;
	vector<int> densityDist;

	CoherenceBasedEEMsearch *parent;

	string centerGene;

public:
	virtual string toString();
protected:
	virtual void findCenter();
	virtual void refineCenter() throw ();
	virtual void refineCenter0()throw ();

public:
	virtual void findModuleGenes() throw ();

protected:
	virtual vector<string> getGlobalCoherentGenes();
	virtual void generateNullDist(int n);

public:
	virtual void calculatePvalue1();
	virtual void calculatePvalue2();

protected:
	virtual map<int, double> getDensityDist();

public:
	//--
	// inline
	//--
	virtual vector<string> getModuleGenes() {return moduleGenes;}
	virtual void  cutParent(){parent = NULL;}

	virtual double  getAbsoluteRadius(){return absoluteRadius;};
	virtual double  getRelativeRadius(){return relativeRadius;};
	virtual vector<double>  getCenter(){return center;};
	virtual string getCenterGene(){return centerGene;};
	virtual void   setAbsoluteRadius(double r){absoluteRadius = r;}

	virtual double  getPvalue(){return Pvalue;}
	virtual double  getPvalue1(){return Pvalue1;}
	virtual double  getPvalue2(){return Pvalue2;}
	virtual void  setPvalue2(double p){Pvalue2 = p;}
	virtual void  setPvalue1(double p){Pvalue1 = p;}
	virtual void  setPvalue(double p){Pvalue = p;}
	virtual vector<string>  getSeedGenes(){return seedGenes;}

public:
	// add method
	virtual string name(){ return "CoherenceBasedEEM";}
};
} // namecpace eem
#endif // CoherenceBasedEEM_H
