//
// MyFunc.h
//
#ifndef MyFunc_H
#define MyFunc_H

#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
using namespace std;

namespace utility {

class MyFunc{
public:
	MyFunc();
	virtual ~MyFunc();

	static double mean(const vector<double> &v);
	static double sd(const vector<double> &v);
	static double var(const vector<double> &v);
	static double abs(const vector<double> &v);
	static double innerProduct(const vector<double> &v, const vector<double> &u);
	static double pearsonCorrelation(const vector<double> &v, const vector<double> &u);
	static double euclideanDist(const vector<double> &v, const vector<double> &u);

	static double pearsonCorrelationForNormarizedList(const vector<double> &v, const vector<double> &u);

	static double max(const vector<double> &v);
	static double min(const vector<double> &v);
	static double median(const vector<double> &v);
	static double mode(const vector<double> &v);

	static double percentile(const vector<double> &v, double d);
	static double upperQuartile(const vector<double> &v);
	static double lowerQuartile(const vector<double> &v);
	static double sum(const vector<double> &v);
	static vector<double> normalize( const vector<double> &v);
	static double tStatistic( const vector<double> &v, const vector<double> &u );

	static double ANOVAStatistics(const vector<vector<double> > &X);

	static string join(const string &d, const vector<string> &s);
	static vector<string> uniq(const vector<string> &O);
	static vector<string> isect(const vector<string> &A, const vector<string> &B);
	static vector<string> sample(const vector<string> &S, int n);

	static double calculatePvalueForSetOverlap(int t, int a, int b, int u);
	static vector<string> sortKeysByAscendingOrderOfValues(const map<string, double>& org_map);
	static vector<string> sortKeysByDescendingOrderOfValues(const map<string, double>& org_map);
	static vector<vector<string> > getAllCombination(const vector<string>& v, int n);

	static vector<string> split(const string &line, char delim);
	static map<string, vector<string> > readGeneSetFromGmtFile(const string &gmtfile);

	// add function
	static bool contains(const vector<string>& vec_large, const vector<string>& vec_small);
	static double parseDouble(const string& str);
	static string toString(const double& val);
	static string toString(const int& val);
	static void erase(vector<string>& vstr, const string& s);
};

} // namespace utility

#endif // MyFunc_H
