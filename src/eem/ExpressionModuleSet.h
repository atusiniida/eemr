//
// ExpressionModuleSet.h
//
#ifndef ExpressionModuleSet_H
#define ExpressionModuleSet_H

#include <exception>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
using namespace std;

#include "../utility/MyException.h"
using namespace utility;
//
// class ExpressionModuleSet implements Serializable
//
namespace eem{
class ExpressionModule;

///
/// an ExpressionModuleSet class object is a container of ExpressionModule class objects
///
class ExpressionModuleSet{
public:
	ExpressionModuleSet();
	virtual ~ExpressionModuleSet();

private:

	map<string, ExpressionModule*> expressionModules;
	vector<string> ids;

public:
	// add expressionModule
	// if an expressionModule with same ID exists, set the better one.
	void add(ExpressionModule* expMod);

	ExpressionModule* get(string id);
	ExpressionModule* get(int i);

	vector<string> getIds();

	void writeToFile(string file) throw ();
};
} // namespace eem
#endif // ExpressionModuleSet_H
