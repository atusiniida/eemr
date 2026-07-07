//
// CoherenceBasedEEMsearch.h
//
#ifndef CoherenceBasedEEMsearch_H
#define CoherenceBasedEEMsearch_H

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
using namespace std;

#include "../utility/Dist.h"
#include "../utility/MyException.h"
#include "../utility/MyFunc.h"
#include "../utility/MyMat.h"
using namespace utility;

#include "AbstractEEMsearch.h"

namespace eem{
class DistConverter;
class CoherenceBssedEEM;
    
///
/// an CoherenceBasedEEMsearch class objects performs EEM analysis
///
/// inheritance relationships are:
/// CoherenceBasedEEMsearch->AbstractEEMsearch->AbstractGenesetAnalysis
///
class CoherenceBasedEEMsearch : public AbstractEEMsearch{
public:
	CoherenceBasedEEMsearch();
	virtual ~CoherenceBasedEEMsearch();


	MyMat Exp; // expression matrix
	Dist  Cor; // correlation half matrix

	// these are parameters for EEM search
	double absoluteRadius;
	double relativeRadius; // the parameter to be parallelized
	int coreGeneSize;

	DistConverter *distConverter;  // use for parameter conversion (relativeRadius <-> absoluteRadius)

	//--
	// inner class
	//--
	class DistFunc {
	public:
		DistFunc(){}
		DistFunc(Dist* cor){ Cor = cor;}
		virtual ~DistFunc(){}

		Dist* Cor;//------ 2013.11.21 pointer

		virtual double get(vector<double> a, vector<double> b){ return 1;}//=0;
		virtual double get(string a, string b){ return 1;}//=0;
	};
	class Correlation: public DistFunc {
		public:
			Correlation(){}
			Correlation(Dist* cor){ Cor = cor;}
			virtual ~Correlation(){}

			virtual double get(vector<double> a, vector<double> b){
				return 1.0 - MyFunc::pearsonCorrelationForNormarizedList(a,b);
			}
			virtual double get(string a, string b){
				return 1.0 - Cor->get(a,b);
			}
	};

	//-----------------
	DistFunc *distfunc;
	//-----------------

	CoherenceBasedEEMsearch(MyMat& E);
	
protected:
	// prepare for EEM search
	// here a correlation matrix is calculated (to be parallelized) and parameter conversion performed
	virtual void setEEM();// generate CoherenceBasedEEM

public:
	virtual void setCoreGeneSize(int i);

	virtual void setCor(Dist& D);
	// a correlation matrix is calculated (to be parallelized)
	virtual void calculateCor();

	// parameter setting and conversion (absoluteRadius <-> relativeRadius)
	virtual void setAbsoluteRadius(double r);
	virtual void setRelativeRadius(double r);

	virtual string toString();

};
} // namespace eem
#endif // CoherenceBasedEEMsearch_H
