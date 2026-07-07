//
// AbstractEEMsearch.cpp
//
#include "AbstractEEMsearch.h"
using namespace eem;

AbstractEEMsearch::AbstractEEMsearch()
{
	maxSeedGeneSize = 2000;
	minSeedGeneSize = 10;
	itrForPvalue2Calculation = 500;
	Pvalue1Cutoff = 1.0;

	//stopWatch = new StopWatch();
	//seedGeneSets = new HashMap<String, List<String>>();
	//eems = new HashMap<String, EEM>();
	//seeds = new ArrayList<String>();
	//candidates  = new ArrayList<String>();

	errLog = "";
	timeLog = "";

	nullDistrbutionForRecycle.clear();
}
AbstractEEMsearch::~AbstractEEMsearch()
{
	//TODO : change safe delete
	map<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){
		EEM *eem = it->second;
		delete eem;
	}
}

AbstractEEMsearch::AbstractEEMsearch(AbstractEEMsearch &A)
{
	originalExp = A.originalExp;
	allGenes = A.allGenes;
	maxSeedGeneSize = A.maxSeedGeneSize;
	minSeedGeneSize = A.minSeedGeneSize;
	itrForPvalue2Calculation = A.itrForPvalue2Calculation;
	Pvalue1Cutoff = A.Pvalue1Cutoff;

	seedGeneSets = A.seedGeneSets;

	errLog = "";
	timeLog = "";
}



void AbstractEEMsearch::setGeneSets( map<string, vector<string> > geneSets )
{
	map<string, vector<string> >::iterator it;
	for(it=geneSets.begin(); it != geneSets.end(); it++){

		vector<string> tmp = MyFunc::isect(it->second, allGenes);

		if(tmp.size() < minSeedGeneSize || tmp.size() > maxSeedGeneSize){
			stringstream ss;
			ss << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!";
			cerr << ss.str() << endl;

			ss.str("");
			ss.clear(stringstream::goodbit);
			ss << it->first << ": seed geneset size (" << tmp.size() << ") is out of range!\n";
			errLog += ss.str();

			continue;
		}
		seedGeneSets[it->first] = tmp;
	}
}

void AbstractEEMsearch::setPvalue1Cutoff(double d)
{
	Pvalue1Cutoff = d;
}

void AbstractEEMsearch::suppressPvalue1Cutoff()
{
	Pvalue1Cutoff = 0.0;
}

void AbstractEEMsearch::setItrForPvalue2Calculation(int i) 
{
	itrForPvalue2Calculation= i;
}

void AbstractEEMsearch::recycleNullDistribution()
{
	nullDistrbutionForRecycle.clear();// = new HashMap<Integer, List<Integer>>();
}

//to be parallelized
void AbstractEEMsearch::findModuleGenes()
{
	cerr << "Findng module gene..." << endl;

	int i = 0;
	int n = eems.size();
	vector<string> new_seeds;// = new ArrayList<String>();

	for(; i < seeds.size(); i++){
		string s = seeds[i];
		try{
			eems[s]->findModuleGenes();// EEM*

			stringstream ss;
			ss << s << "(" << i+1 << "/" << n << "): succeed! " << eems[s]->getModuleGenes().size() << "/" << eems[s]->getSeedGenes().size();
			cerr << ss.str() << endl;

			new_seeds.push_back(s);
		}
		catch (MyException &err) {
			stringstream ss;
			ss << s << "(" << i+1 << "/" << n << "): failed!";

			cerr << err.what_message(ss.str()) << endl;

			errLog += s +  ": unable to find module genes!\n"; 
			continue;
		}
	}
	seeds = new_seeds;

	if(seeds.empty()){
		//throw new MyException("findModuleGenes: Any module has passed!");
		throw MyException();// TODO
	}

}

//to be parallelized
void AbstractEEMsearch::calculatePvalue1()
{
	if(Pvalue1Cutoff == 0.0){
		return;
	}else{
		cerr << "Calculating approximate P values to filter out non-significant genesets..." << endl;
		int i = 0;
		int n = seeds.size();
		vector<string> new_seeds;// = new ArrayList<String>();

		for(; i < seeds.size(); i++){
			string s = seeds[i];
			try{
				eems[s]->calculatePvalue1();

				stringstream ss;
				ss << s << "(" << i+1 << "/" << n << "): " << "succeed! " << eems[s]->getPvalue1();
				cerr << ss.str() << endl;

				new_seeds.push_back(s);
			}
			catch (MyException &err) {
				stringstream ss;
				ss << s << "(" << i << "/" << n << "): failed!";

				cerr << err.what_message(ss.str()) << endl;

				errLog +=  s +  ": unable to calculate a P value  based on hypergeometric distribution!\n";
				continue;
			}
		}

		seeds = new_seeds;

		if(seeds.empty()){
			//throw new MyException(" calculatePvalue1: Any module has passed!");
			throw MyException();// TODO
		}
	}
}

void AbstractEEMsearch::findCandidates()
{
	if(Pvalue1Cutoff == 0.0){
		candidates = seeds;// candidates = new ArrayList<String>(seeds);
	}else{
		for(size_t i=0; i < seeds.size(); i++){
			string s = seeds[i];
			if(eems[s]->getPvalue1() >= Pvalue1Cutoff){
				candidates.push_back(s);
			}
		}
		if(candidates.empty()){
			//throw MyException("findCandidates: Any module has passed!");
			cerr << "findCandidates: Any module has passed!" << endl;
		}
	}
}

//to be parallelized
void AbstractEEMsearch::calculatePvalue2()
{
	cerr << "Calculating accurate P values..." << endl;

	int i = 0;
	int n = candidates.size();
	vector<string> new_candidates;// = new ArrayList<String>();

	for(; i < n; i++){
		string s;
		try{
			s = candidates[i];
			eems[s]->calculatePvalue2();

			stringstream ss;
			ss << s << "(" << i+1 << "/"  << n << "): " << "succeed! " << eems[s]->getPvalue2();
			cerr << ss.str() << endl;

			new_candidates.push_back( s );
		}
		catch (MyException &err) {
			// err.printStackTrace();
			stringstream ss;
			ss << s << "(" << i+1 << "/" << n << "): failed!";
			cerr << err.what_message( ss.str() ) << endl;

			errLog +=  s +  ": unable to calculate a P value based on arandomization test!\n"; 
			continue;
		}
	}

	candidates = new_candidates;
	if(candidates.empty()){
		//throw MyException("calculatePvalue2: Any module has passed!");
		cerr << "calculatePvalue2: Any module has passed!" << endl;
	}

}

void AbstractEEMsearch::printResults(string outfile) throw ()
{
	map<string, double> p;

	map<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){
		p[it->first] = it->second->getPvalue();
	}
	// tmp <= p
	vector<string> tmp = MyFunc::sortKeysByDescendingOrderOfValues( p );

	// tmp => outfile
	ofstream ofs;
	ofs.open(outfile.c_str());
	size_t i;
	for( i=0; i < tmp.size(); i++ ){
		string id = tmp[i];

		stringstream ss;
		ss << id << "\t" << eems[id]->getModuleGenes().size() << "/" << eems[id]->getSeedGenes().size()
			<< "\t" << eems[id]->getPvalue() << "\t" << MyFunc::join("\t", eems[id]->getModuleGenes());

		ofs << ss.str() << endl;
	}
	ofs.close();
}

void AbstractEEMsearch::printLog(string outfile) throw ()
{
	ofstream ofs;

	ofs.open(outfile.c_str(), ios::out|ios::app);
	ofs << getLog() << endl;
	ofs.close();
}

void AbstractEEMsearch::setEEM()
{
	//
	// throw new UnsupportedOperationException();
	//
	cerr << "Unsupported Operation  Exception : void AbstractEEMsearch::setEEM()" << endl; 
}
//--
// generate ExpressionModuleSet, ExpressionModule
//--
ExpressionModuleSet* AbstractEEMsearch::getExpressionModuleSet()
{
	ExpressionModuleSet* tmp = new ExpressionModuleSet();// generate ExpressionModuleSet

	map<string, EEM*>::iterator it;

	for(it=eems.begin(); it != eems.end(); it++){
		if(it->second->getPvalue() >= 0.0){
			it->second->cutParent();
			tmp->add(new ExpressionModule(it->first, it->second, originalExp));// generate ExpressionModule
		}
	}
	return tmp;
}

string AbstractEEMsearch::getLog()
{
	return timeLog + "\n" + this->name() + "\n" + errLog;
}

void AbstractEEMsearch::writeTimeLog()
{
	timeLog = "started at " +  stopWatch.getStartDate().toString() + "\n" +
			"finished at " + stopWatch.getStopDate().toString() + "\n" +
			"It took " + stopWatch.toString() + "\n";
}

void AbstractEEMsearch::perform()
{
	try{
		setEEM();
		findModuleGenes();
		calculatePvalue1();
		findCandidates();
		calculatePvalue2();

	}catch(MyException& err){
		cerr << err.what() << endl;

		stopWatch.stop();
		writeTimeLog();
	}
}

map<string, double> AbstractEEMsearch::getPvalues()
{
	map<string, double> P;

	map<string, EEM*>::iterator it;
	for(it=eems.begin(); it != eems.end(); it++){

		if(it->second->getPvalue() >= 0.0){
			P[it->first]= it->second->getPvalue();
		}
	}

	return P;
}
