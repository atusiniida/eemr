//
// DistConverter.h
//
#ifndef DistConverter_H
#define DistConverter_H

#include <map>
#include <string>
#include <iostream>
#include <cmath>
using namespace std;

#include "../utility/MyFunc.h"
using namespace utility;

namespace eem{
class CoherenceBasedEEMsearch;
class CoherenceBasedEEM;
class DistConverter{
public:
	DistConverter();
	DistConverter(CoherenceBasedEEMsearch* parent);
	virtual ~DistConverter();

private:
	int maxItrForBisection;
	int maxItrForBisection2;
	double deltaForBisection;
	double upperAbsoluteDist;
	double lowerAbsoluteDist;
	CoherenceBasedEEMsearch *parent;
	CoherenceBasedEEM       *e;
	map<string, double> absolute2relative;

public:
	double convertAbsolute2relativeDist(double absoluteDist);// generate CoherenceBasedEEM
	double convertRelative2absoluteDist(double relativeDist);

};
} // namespace eem
#endif // DistConverter_H
