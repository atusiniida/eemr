#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>

#include "eem/CoherenceBasedEEMsearch.h"
#include "eem/CoherenceBasedEEM.h"
#include "eem/ExpressionModuleSet.h"
#include "eem/ExpressionModule.h"
#include "eem/EEM.h"
#include "utility/MyMat.h"

using namespace Rcpp;
using namespace eem;
using namespace utility;
using namespace std;

static string collapse(const vector<string> &v, const string &sep)
{
	string out;
	for (size_t i = 0; i < v.size(); i++) {
		if (i > 0) out += sep;
		out += v[i];
	}
	return out;
}

// [[Rcpp::export]]
DataFrame eem_search_cpp(NumericMatrix expr,
                          List geneSets,
                          double relativeRadius,
                          double pvalue1Cutoff,
                          int itr,
                          int minGeneSetSize,
                          int maxGeneSetSize,
                          int coreGeneSize)
{
	RObject dimnamesAttr = expr.attr("dimnames");
	if (dimnamesAttr.isNULL()) {
		stop("expr must have rownames (gene ids)");
	}
	List dimnames(dimnamesAttr);
	RObject rnAttr = dimnames[0];
	if (rnAttr.isNULL()) {
		stop("expr must have rownames (gene ids)");
	}
	CharacterVector rn(rnAttr);
	int nr = expr.nrow();
	int nc = expr.ncol();

	vector<string> rownames_vec(nr);
	for (int i = 0; i < nr; i++) rownames_vec[i] = as<string>(rn[i]);

	vector<string> colnames_vec(nc);
	RObject cnAttr = dimnames[1];
	if (!cnAttr.isNULL()) {
		CharacterVector cn(cnAttr);
		for (int j = 0; j < nc; j++) colnames_vec[j] = as<string>(cn[j]);
	} else {
		for (int j = 0; j < nc; j++) colnames_vec[j] = "V" + std::to_string(j + 1);
	}

	vector<vector<double> > data(nr, vector<double>(nc));
	for (int i = 0; i < nr; i++) {
		for (int j = 0; j < nc; j++) {
			data[i][j] = expr(i, j);
		}
	}

	MyMat mymat(rownames_vec, colnames_vec, data);

	CoherenceBasedEEMsearch CE(mymat);

	map<string, vector<string> > gs;
	CharacterVector gsNames = geneSets.names();
	if (gsNames.isNULL()) {
		stop("geneSets must be a named list (gene set id -> character vector of gene ids)");
	}
	for (int i = 0; i < geneSets.size(); i++) {
		string name = as<string>(gsNames[i]);
		CharacterVector genes = geneSets[i];
		vector<string> genes_vec(genes.size());
		for (int k = 0; k < genes.size(); k++) genes_vec[k] = as<string>(genes[k]);
		gs[name] = genes_vec;
	}
	CE.setGeneSets(gs);

	CE.setMinGeneSetSize(minGeneSetSize);
	CE.setMaxGeneSetSize(maxGeneSetSize);
	if (coreGeneSize > 0) CE.setCoreGeneSize(coreGeneSize);
	if (pvalue1Cutoff >= 0.0) {
		CE.setPvalue1Cutoff(pvalue1Cutoff);
	} else {
		CE.suppressPvalue1Cutoff();
	}
	CE.setItrForPvalue2Calculation(itr);
	CE.setRelativeRadius(relativeRadius);

	CE.perform();

	ExpressionModuleSet *ems = CE.getExpressionModuleSet();
	vector<string> ids = ems->getIds();
	size_t n = ids.size();

	CharacterVector out_id(n);
	NumericVector out_pvalue(n), out_pvalue1(n), out_pvalue2(n);
	NumericVector out_absR(n), out_relR(n);
	IntegerVector out_nSeed(n), out_nModule(n);
	CharacterVector out_seedGenes(n), out_moduleGenes(n);

	for (size_t i = 0; i < n; i++) {
		ExpressionModule *em = ems->get(ids[i]);
		EEM *e = em->getEEM();

		vector<string> seedG = e->getSeedGenes();
		vector<string> modG = e->getModuleGenes();

		out_id[i] = ids[i];
		out_pvalue[i] = em->getPvalue();
		out_pvalue1[i] = e->getPvalue1();
		out_pvalue2[i] = e->getPvalue2();
		out_absR[i] = e->getAbsoluteRadius();
		out_relR[i] = e->getRelativeRadius();
		out_nSeed[i] = (int)seedG.size();
		out_nModule[i] = (int)modG.size();
		out_seedGenes[i] = collapse(seedG, ";");
		out_moduleGenes[i] = collapse(modG, ";");
	}

	delete ems;

	DataFrame df = DataFrame::create(
		Named("id") = out_id,
		Named("pvalue") = out_pvalue,
		Named("pvalue1") = out_pvalue1,
		Named("pvalue2") = out_pvalue2,
		Named("nSeedGenes") = out_nSeed,
		Named("nModuleGenes") = out_nModule,
		Named("absoluteRadius") = out_absR,
		Named("relativeRadius") = out_relR,
		Named("seedGenes") = out_seedGenes,
		Named("moduleGenes") = out_moduleGenes,
		Named("stringsAsFactors") = false
	);
	return df;
}
