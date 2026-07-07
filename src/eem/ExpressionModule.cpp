//
// ExpressionModule.cpp
//
#include "ExpressionModuleSet.h"
#include "ExpressionModule.h"
using namespace eem;


ExpressionModule::~ExpressionModule()
{
}

ExpressionModule::ExpressionModule(ExpressionModule& e)
{
	id = e.id;
	eem = e.eem;
	activityProfile = e.activityProfile;
	moduleGenes = e.moduleGenes;
	seedGenes = e.seedGenes;
	expressionModuleCluster = e.expressionModuleCluster;
	minimunPvalue = e.minimunPvalue;
	PvalueCorrectedForMultipleTest = e.PvalueCorrectedForMultipleTest;
	numberOfMultipleTesting = e.numberOfMultipleTesting;
}

ExpressionModule::ExpressionModule(const string& id, EEM *eem, MyMat& Exp)
{
	this->id = id;
	this->eem = eem;
	minimunPvalue = eem->getPvalue();
	PvalueCorrectedForMultipleTest = eem->getPvalue();
	numberOfMultipleTesting = 1;
	seedGenes = eem->getSeedGenes();
	moduleGenes = eem->getModuleGenes();
}


string ExpressionModule::getId()
{
	return id;
}

double ExpressionModule::getPvalue()
{
	return PvalueCorrectedForMultipleTest;
}

EEM* ExpressionModule::getEEM()
{
	return eem;
}

string ExpressionModule::toString()
{
	string str = id + "\t" + MyFunc::toString(PvalueCorrectedForMultipleTest) +
		"\t" + MyFunc::toString(numberOfMultipleTesting) +
		"\t" + getEEM()->toString();

	return str;
}

//--
// generate ExpressionModule
//--
ExpressionModule* ExpressionModule::getBetterExpressionModules(ExpressionModule& e1, ExpressionModule& e2)
{
	ExpressionModule* e3;
	if(e1.minimunPvalue > e2.minimunPvalue){
		e3 = new ExpressionModule(e1);//----------- generate ExpressionModule
		e3->minimunPvalue = e1.minimunPvalue;
	}else{
		e3 = new ExpressionModule(e2);//----------- generate ExpressionModule
		e3->minimunPvalue = e2.minimunPvalue;
	}
	e3->numberOfMultipleTesting = e1.numberOfMultipleTesting + e2.numberOfMultipleTesting;
	if(e3->minimunPvalue != DBL_MAX){
		double p = e3->minimunPvalue;
		if(p < 5){
			p = pow(10, -p);
			p = 1- pow(1-p, e3->numberOfMultipleTesting);
			p  = - log10(p);
		}else{
			p = p - log10( (double)e3->numberOfMultipleTesting );
		}
		e3->PvalueCorrectedForMultipleTest = p;
	}
	return e3;
		
}


namespace eem {
ostream &operator <<(ostream &stream, ExpressionModule &em)
{
	stream << em.toString();
	return stream;
}
} // namespace eem
