//
// Dist.cpp
//
#include "Dist.h"
#include "MyFunc.h"
using namespace utility;

Dist::Dist()
	: diagonalElement(0)
{
	diagonalElement= 0;
	b_argCon = false;
}
/*type must be 'e'(euclideanDist), 'c' (pearsonCorrelation),
* or 'C' (pearsonCorrelationForNormarizedList)
*/
//Dist::Dist(const MyMat& m, char type)
void Dist::call_argumentConstruct(const MyMat& m, char type)
{
	diagonalElement= 0;
	b_argCon = true;

	n = m.rowSize();
	name = m.getRowNames();
	int i,j;
	for(i =0; i< n; i++){
		name2index[name[i]] = i;
	}
	M.resize(n);
	for(i=0; i < n; i++){
		M[i].resize(i);
	}

	switch(type){
		/*case('c'):
			{
			MyMat copy_m = m;
			copy_m.normalizeRows();
			for(i = 0; i < n; i++){
				for(j = 0; j < i; j++){
					M[i][j] = MyFunc::pearsonCorrelationForNormarizedList(copy_m.getRow(i),copy_m.getRow(j));
			    }
			}
			setDiagonalElement(1.0);
			}
			break;*/
		case('C'):
			//to be parallelized
			for(i = 0; i < n; i++){
				for(j = 0; j < i; j++){
					M[i][j] = MyFunc::pearsonCorrelationForNormarizedList(m.getRow(i),m.getRow(j));
			    }
			}
			setDiagonalElement(1);
			break;
		/*case('e'):
			for(i = 0; i < n; i++){
				for(j = 0; j < i; j++){
					M[i][j] = MyFunc::euclideanDist(m.getRow(i),m.getRow(j));
			    }
			}
			break;*/
		default:
			//cerr << "type must be 'c','C', or 'e'" << endl;
			cerr << "type must be 'C'"<< endl;
			break;
	}
}

Dist::~Dist()
{
}
//----
//
//----
const vector<string> &Dist::getNames()
{
	return name;
}

void Dist::setDiagonalElement(double d)
{
	diagonalElement = d;	
}

double Dist::get(const string& s, const string& t)
{
	int i = name2index[s];
	int j = name2index[t];
	if(i > j){
		return M[i][j];
	}
	if(j > i){
		return M[j][i];
	}
	return diagonalElement;
}


//--
// operator 2013.11.16
//--
Dist& Dist::operator =(const Dist& dist)
{
	M = dist.M;
	name2index = dist.name2index;	
	n = dist.n;
	name = dist.name;

	diagonalElement = dist.diagonalElement;

	b_argCon = dist.b_argCon;//------------- add 2013.11.20

	return *this;
}

//--
// add method
//--
//void Dist::call_argumentConstruct(const MyMat& m, char type)
Dist::Dist(const MyMat& m, char type)
{
	//Dist(m, type);
	call_argumentConstruct(m, type);
}

void Dist::clear()
{
	size_t nRow = M.size();
	for(size_t i=0; i < nRow; i++) M[i].clear();
	M.clear();

	name2index.clear();
	n = 0;
	name.clear();

	diagonalElement = 0.0;

	b_argCon=false;
}
