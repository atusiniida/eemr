///
///EEM.h
///

#ifndef EEM_H
#define EEM_H

#include <map>
#include <vector>
#include <string>
using namespace std;

namespace eem{

///
///the EEM class is  an abstract  class of CoherenceBasedEEM
///
    
class EEM{
public:
	EEM();
	virtual ~EEM();

	virtual vector<string> getSeedGenes(){ vector<string> vstr; return vstr;}//=0;
	virtual vector<string> getModuleGenes(){ vector<string> vstr; return vstr;}//=0;
	virtual double getPvalue1(){ return 0.0;}//=0;
	virtual double getPvalue2(){ return 0.0;}//=0;
	virtual double getPvalue(){ return 0.0;}//=0;
	virtual void setPvalue2(double p){ ;}//=0;
	virtual void setPvalue1(double p){ ;}//=0;
	virtual void setPvalue(double p){ ;}//=0;
	virtual void findModuleGenes() throw (){ ;}//=0;
	virtual void calculatePvalue1(){ ;}//=0;
	virtual void calculatePvalue2(){ ;}//=0;
	virtual void cutParent(){ ;}//=0;

	//
	// add interface method
	//
public:
	virtual string toString(){ string str;   return str;}
protected:
	virtual void findCenter(){ ;}
	virtual void refineCenter() throw () { ;}
	virtual void refineCenter0()throw () { ;}

	virtual vector<string> getGlobalCoherentGenes(){ vector<string> vstr;   return vstr;}
	virtual void generateNullDist(int n){ ;}

	virtual map<int, double> getDensityDist(){ map<int, double> m;   return m;}
public:
	virtual double  getAbsoluteRadius(){ double val=0.0;  return val;}
	virtual double  getRelativeRadius(){ double val=0.0;  return val;}
	virtual vector<double>  getCenter(){ vector<double> vec; return vec;}
	virtual string getCenterGene()     { string str;  return str;}
	virtual void   setAbsoluteRadius(double r){ ;}

	//
	// add method
	//
public:
	virtual string name(){ return "EEM";}
};
} // namespace eem
#endif // EEM_H
