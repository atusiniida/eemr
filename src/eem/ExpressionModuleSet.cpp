//
// ExpressionModuleSet.cpp
//
#include "ExpressionModule.h"
#include "ExpressionModuleSet.h"
using namespace eem;


ExpressionModuleSet::ExpressionModuleSet()
{
}
ExpressionModuleSet::~ExpressionModuleSet()
{
	// TODO: safe
	//
	map<string, ExpressionModule*>::iterator it;

	for(it=expressionModules.begin(); it != expressionModules.end(); it++){
		delete it->second;
	}
}


void ExpressionModuleSet::add(ExpressionModule* expMod)
{
	string id = expMod->getId();

	if( expressionModules.find(id) == expressionModules.end()){
		expressionModules[id] = expMod;
		if(ids.empty()){
			ids.push_back(id);
		}else{
			for(int i = 0, n = ids.size(); i < n;i++){
				if(expMod->getPvalue() > get(i)->getPvalue()){
					ids.insert(ids.begin()+i, id);
					return;
				}
			}
			ids.push_back(id);
			return;
		}
	}else{
		ExpressionModule *e = ExpressionModule::getBetterExpressionModules(*expressionModules[id], *expMod);
		expressionModules[id] = e;

		MyFunc::erase(ids, id);

		for(int i = 0, n = ids.size(); i < n;i++){
			if(e->getPvalue() > get(i)->getPvalue()){
				ids.insert(ids.begin()+i, id);
				return;
			}
		}
		ids.push_back(id);
		return;
	}
}

ExpressionModule* ExpressionModuleSet::get(string id)
{
	return expressionModules[id];
}

ExpressionModule* ExpressionModuleSet::get(int i)
{
	return expressionModules[ids[i]];
}

vector<string> ExpressionModuleSet::getIds()
{
	return ids;
}


