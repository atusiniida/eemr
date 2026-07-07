//
// Dist.h
//
#ifndef Dist_H
#define Dist_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
using namespace std;

#include "MyFunc.h"
#include "MyMat.h"
#include "MyException.h"

namespace utility{

class Dist{
public:
	Dist();
	/*type must be 'e'(euclideanDist), 'c' (pearsonCorrelation),
	 * or 'C' (pearsonCorrelationForNormarizedList)
	 */
	Dist(const MyMat& m, char type);

	virtual ~Dist();

private:
	static long long serialVersionUID;// = -1890276595674423636L;

	vector<vector<double> >    M;
	map<string, int>  name2index;
	int n;
	vector<string> name;

	double diagonalElement;

	bool b_argCon;//use Dist(MyMat&, char) ?

public:
	const vector<string>& getNames();
	void setDiagonalElement(double d);
	double get(const string& s, const string& t);

	static Dist readFromBinary(const string& infile) throw ();// FileNotFoundException, IOException, ClassNotFoundException

	//--
	// operator 2013.11.16
	//--
	Dist& operator=(const Dist& dist);

	//--
	// add method 2013.11.16
	//--
	bool use_argumentConstruct(){ return b_argCon;}
	void call_argumentConstruct(const MyMat& m, char type);

	void clear();//------------2013.11.20
};
} // namespace utility
#endif // Dist_H
