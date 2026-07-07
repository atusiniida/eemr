//
// MyMat.h
//
#ifndef MyMat_H
#define MyMat_H

#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "MyFunc.h"
#include "MyException.h"

namespace utility{
class MyMat{
public:
	MyMat();
	MyMat(const string &infile);// throws MyException on I/O or format errors
	MyMat(const vector<string> &rownames, const vector<string> &colnames, const vector<vector<double> > &data);// build directly from in-memory data (used by the R/Rcpp binding)
	MyMat(const MyMat& m);
	virtual ~MyMat();

private:
	static long long serialVersionUID;// = -6849257508515306759L;

protected:
	vector<vector<double> > M;
	map<string, int> colname2index;
	map<string, int> rowname2index;
	int ncol;
	int nrow;
	vector<string> colname;
	vector<string> rowname;

public:
	double get(int i, int j) const;
	const vector<string>& getRowNames() const;
	const vector<string>& getColNames() const;

	int rowSize() const;
	int colSize() const;

	vector<double> getRow(int i) const;
	vector<double> getRow(const string &i);

	void normalizeRows();

	vector<double> getRowMeans(const vector<string> &rows);

	//--
	// operator 2013.11.16
	//--
	MyMat& operator=(const MyMat& mat);
};
} // namespace utility

#endif // MyMat_H
