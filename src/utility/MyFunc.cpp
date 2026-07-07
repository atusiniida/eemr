//
// MyFunc.cpp
//
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
using namespace std;

#include "MyFunc.h"
#include "HypergeometricDistributionImpl.h"
#include "MyException.h"
using namespace utility;

MyFunc::MyFunc()
{
}
MyFunc::~MyFunc()
{
}

double MyFunc::mean(const vector<double> &v){
	double m = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		m += v[i];
	}
	return m/v.size();
}
double MyFunc::sd(const vector<double> &v){
	return sqrt(var(v));
}
double MyFunc::var(const vector<double> &v){
	double var = 0;
	double m = mean(v);
	for(int i = 0, n = v.size(); i < n; i++){
		double tmp = v[i] - m;
		var += tmp * tmp;
	}
	var /= (double)v.size();
	return var;
}
double MyFunc::abs(const vector<double> &v){
	double tmp = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		tmp += v[i] * v[i];
	}
	return sqrt(tmp);
}
double MyFunc::innerProduct(const vector<double> &v, const vector<double> &u){
	double tmp = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		tmp += v[i] * u[i];
	}
	return tmp;
}
double MyFunc::pearsonCorrelation(const vector<double> &v, const vector<double> &u){
	double c = 0;
	double vm = mean(v);
	double um = mean(u);
	for(int i = 0, n = v.size(); i < n; i++){
		c +=  (v[i]-vm) * (u[i]-um);
	}
	c /= ((double)v.size())*sd(v)*sd(u);
	return c;
}
double MyFunc::euclideanDist(const vector<double> &v, const vector<double> &u){
	double d = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		double tmp = v[i]-u[i];
		d += tmp * tmp;
	}
	return sqrt(d);
}

double MyFunc::pearsonCorrelationForNormarizedList(const vector<double> &v, const vector<double> &u){
	double c = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		c +=  v[i]*u[i];
	}
	c /= (double)v.size();
	return c;
}

double MyFunc::max(const vector<double> &v){
	double  m = v[0];
	for(int i = 0, n = v.size(); i < n; i++){
		if(v[i] > m){
			m = v[i];
		}
	}
	return m;
}
double MyFunc::min(const vector<double> &v){
	double  m = v[0];
	for(int i = 0, n = v.size(); i < n; i++){
		if(v[i] <  m){
			m = v[i];
		}
	}
	return m;
}
double MyFunc::median(const vector<double> &v){
	vector<double> u(v);
	sort(u.begin(), u.end());
	int n = u.size();
	if(n%2 == 1){
		return u[(n-1)/2];
	}else{
		return (u[n/2] + u[(n/2)-1])/2;
	}
}
double MyFunc::mode(const vector<double> &v){
	double M = max(v);
	double m = min(v);
	double binNum = 1000;
	double delta = (M-m)/binNum;
	vector<double> u(v);
	sort(u.begin(), u.end());

	int count = 0;
	int countMax = 0;
	double  d = m+delta;
	double  dMax = d;
	for(int i = 0, n = u.size(); i < n; i++){
		double e = u[i];
		if(e <= d){
			count++;
		}else{
			if(count > countMax){
				dMax = d;
				countMax = count;
			}
			while(e > d){
				d += delta;
			}
			count = 1;
		}
	}
	dMax = dMax-(0.5*delta);
	return dMax;
}


double MyFunc::percentile(const vector<double> &v, double d){
	vector<double> u(v);
	sort(u.begin(), u.end());
	int n = u.size();
	if(d > 1.0){
		d /= 100;
	}
	int i = (int)floor(n*d)-1;
	if(i < 0){
		i = 0;
	}
	return u[i];
}
double MyFunc::upperQuartile(const vector<double> &v){
	return percentile(v,0.75);
}
double MyFunc::lowerQuartile(const vector<double> &v){
	return percentile(v,0.25);
}
double MyFunc::sum(const vector<double> &v){
	double s = 0;
	for(int i = 0, n = v.size(); i < n; i++){
		s += v[i];
	}
	return s;
}
vector<double> MyFunc::normalize( const vector<double> &v){
	vector<double> tmp(v);
	double m = mean(v);
	double s = sd(v);
	for(int i = 0, n = v.size(); i < n; i++){
		tmp[i] = (tmp[i] - m)/s;
	}
	return tmp;
}
double MyFunc::tStatistic( const vector<double> &v, const vector<double> &u ){
	double V  = sqrt(var(v)/v.size() + var(u)/u.size());
	return (mean(v) - mean(u))/V;
}


double MyFunc::ANOVAStatistics(const vector<vector<double> > &X){
	double SStotal = 0, SSbetween = 0;
	vector<double> all;
	for(int i = 0, n = X.size(); i < n; i++) {
		all.insert(all.end(), X[i].begin(), X[i].end());
	}
	double meanAll = mean(all);
	for(int i = 0, n = X.size(); i < n; i++) {
		double tmp = mean(X[i])-meanAll;
		SSbetween += X[i].size()* tmp * tmp;
		for(int j = 0, m = X[i].size(); j < m; j++){
			tmp = X[i][j]-meanAll;
			SStotal += tmp * tmp;
		}
	}
	double SSwithin = SStotal -SSbetween;
	double MSbetween = SSbetween/((double)X.size()-1);
	double MSwithin = SSwithin/((double)all.size() - X.size());
	return  MSbetween/MSwithin;
}

string MyFunc::join(const string &d, const vector<string> &s){
	string S(s[0]);
	for(int i = 1, n = s.size(); i < n; i++){
		S += d;
		S += s[i];
	}
	return S;
}

vector<string> MyFunc::uniq(const vector<string> &O){
	vector<string> tmp(O);
	sort(tmp.begin(), tmp.end());
	tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
	return tmp;
}

vector<string> MyFunc::isect(const vector<string> &A, const vector<string> &B){
	vector<string> v_isect;
	vector<string> Asort(A), Bsort(B);
	sort(Asort.begin(), Asort.end());
	sort(Bsort.begin(), Bsort.end());
	set_intersection(Asort.begin(), Asort.end(), Bsort.begin(), Bsort.end(),
			 back_inserter(v_isect));
	return v_isect;
}

vector<string> MyFunc::sample(const vector<string> &S, int n){
	vector<string> tmp(S);
	random_shuffle(tmp.begin(), tmp.end());
	tmp.erase(tmp.begin()+n, tmp.end());
	return tmp;
}

/*  Set <?> A , B, Background;
 *   calculatePvalueForSetOverlap( Background.size(), A.size(), B.size(), (isect(A,B)).size);
 */
double MyFunc::calculatePvalueForSetOverlap(int t, int a, int b, int u){
	if(!(t >  0 && a > 0 && b > 0 && u > 0  && t >= a && t >= b && a >= u && b >= u)){
		cerr << "calculatePvalueForSetOverlap: input error!\n";
		throw MyException("calculatePvalueForSetOverlap: input error!");
	}
	int s;
	int l;
	if(a>b){
		l = a;
		s = b;
	}else{
		l = b;
		s = a;
	}
	HypergeometricDistributionImpl H(t, l, s);
	return  H.upperCumulativeProbability(u);

}

vector<string> MyFunc::sortKeysByAscendingOrderOfValues(const map<string,double> &org_map){
	typedef map<string, double> org_map_t;
	typedef map<double, vector<string> > new_map_t;
	new_map_t new_map;
	for (org_map_t::const_iterator it = org_map.begin(); it != org_map.end(); it++) {
		double val = it->second;
		string key(it->first);
		new_map[val].push_back(key);
	}
	vector<string> tmp;
	for (new_map_t::iterator it2 = new_map.begin(); it2 != new_map.end(); it2++) {
		vector<string> &s = it2->second;
		tmp.insert(tmp.end(), s.begin(), s.end());
	}
	return tmp;
}

vector<string> MyFunc::sortKeysByDescendingOrderOfValues(const map<string,double> &org_map){
	vector<string> tmp(sortKeysByAscendingOrderOfValues(org_map));
	reverse(tmp.begin(), tmp.end());
	return tmp;
}

vector<vector<string> > MyFunc::getAllCombination(const vector<string> &v, int n) {
	vector<vector<string> > V;
	if (v.size() <= n) return V;
	set<vector<string> > seen;
	int i;
	for (i = 0; i < v.size(); i++) {
		vector<string> tmp;
		tmp.push_back(v[i]);
		if (seen.find(tmp) != seen.end()) {
			V.push_back(tmp);
			seen.insert(tmp);
		}
	}

	int j, k;
	for (i = 1; i < n; i++) {
		vector<vector<string> > TMP;
		for (j = 0; j < V.size(); j++) {
			for (k = 0; k < v.size(); k++) {
				vector<string> tmp(V[j]);
				tmp.push_back(v[k]);
				sort(tmp.begin(), tmp.end());
				if (seen.find(tmp) != seen.end()) continue;
				if (find(V[j].begin(), V[j].end(), v[k]) == V[j].end()) {
					TMP.push_back(tmp);
					seen.insert(tmp);
				}
			}
		}
		V = TMP;
	}
	return V;
}

vector<string> MyFunc::split(const string &line, char delim) {
	vector<string> ret;
	istringstream iss(line);
	string tmp;
	while (getline(iss, tmp, delim)) {
		ret.push_back(tmp);
	}
	return ret;
}

map<string, vector<string> > MyFunc::readGeneSetFromGmtFile(const string &gmtfile) {
	ifstream inputStream(gmtfile.c_str());
	vector<string> str;
	map<string, vector<string> > GeneSetMap;
	string line;
	while (getline(inputStream, line)) {
		if (line[0] == '#') continue;
		str = split(line, '\t');
		if (str.size() < 3) continue;
		string id = str[0];
		str = MyFunc::uniq(vector<string>(str.begin()+2, str.end()));
		GeneSetMap[id] = vector<string>(str);
	}
	inputStream.close();
	if (GeneSetMap.empty()) {
		cerr << "readGeneSetFromGmtFile: file format is wrong!" << endl;
		throw MyException("readGeneSetFromGmtFile: file format is wrong!");
	}
	return GeneSetMap;
}


//--
// add function
//--
bool MyFunc::contains(const vector<string>& vec_l, const vector<string>& vec_s)
{
	vector<string>::const_iterator it = vec_s.begin();

	for(; it != vec_s.end(); it++){
		if( find(vec_l.begin(), vec_l.end(), *it) == vec_l.end() ){
			return false;
		}
	}
	return true;
}
double MyFunc::parseDouble(const string& str)
{
	double val;

	stringstream ss;
	ss << str;
	ss >> val;

	return val;
}
string MyFunc::toString(const double& val)
{
	string str;

	stringstream ss;
	ss << setprecision(17) << val;
	ss >> str;

	return str;
}
string MyFunc::toString(const int& val)
{
	string str;

	stringstream ss;
	ss << val;
	ss >> str;

	return str;
}
void MyFunc::erase(vector<string>& vstr, const string& s)
{
	vector<string>::iterator it;

	for(it=vstr.begin(); it != vstr.end(); it++){
		if(*it == s) vstr.erase(it);
	}
}
