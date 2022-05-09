#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
using namespace std;

double D_MAX = 1000000000;
double dist(double Ax, double Ay, double Bx, double By)
{
    return sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
}
struct result;
class matrix
{
private:
	double** data;
	int* Di;
	int* Dj;
	unsigned sz;
	void set_M(vector<pair<double, double>> D);
    void set_M(double** D, int n, int* di = NULL, int* dj = NULL);
	double getMin(int I, bool rows = true);
public:
	matrix();
	matrix(const matrix &A);
	matrix(vector<pair<double, double>> D);
    matrix(double** D, int n, int* di = NULL, int* dj = NULL);
	~matrix();
	double st1();
	void set (int eI, int eJ, double S = -1);
	double getResultSum(vector<pair<int, int>> result);
	void clear();
	bool is_empty() const;
    unsigned size() const;
	void erase(int I, int J);
	double* operator[](int i) const;
	static res Run(res A);
    matrix& operator= (const matrix& D);

};

matrix::matrix()
{
	data = NULL;
	sz = 0;
	Di = NULL;
	Dj = NULL;
}

matrix::matrix(const matrix &A)
{
	this->sz = 0;
	this->set_M(A.data, A.sz, A.Di, A.Dj);
}

matrix::matrix(vector<pair<double, double>> D)
{
	this->sz = 0;
	this->set_M(D);
}

matrix::matrix(double** D, int n, int* di, int* dj)
{
	this->sz = 0;
	this->set_M(D, n, di, dj);
}

matrix::~matrix()
{
	this->clear();
}


double* matrix::operator[](int i) const
{
	return data[i];
}

matrix& matrix::operator= (const matrix& D)
{
	this->clear();
	this->set_M(D.data, D.sz, D.Di, D.Dj);
	return *this;
}
void matrix::set_M(double** D, int n, int* di, int* dj)
{
	this->clear();
	sz = n;
	Di = new int [sz];
	Dj = new int [sz];
	data = new double* [sz];
	for (int i = 0; i < sz; ++i)
	{
		data[i] = new double[sz];
		if (di == NULL)
			Di[i] = i + 1;
		else
			Di[i] = di[i];
		if (dj == NULL)
			Dj[i] = i + 1;
		else
			Dj[i] = dj[i];
		for (int j = 0; j < sz; ++j)
			if (D[i][j] == -1) data[i][j] = D_MAX;
			else data[i][j] = D[i][j];
	}
}
  
void matrix::set_M(vector<pair<double, double>> D)
{
	this->clear();
	sz = D.size();
	Di = new int [sz];
	Dj = new int [sz];
	data = new double* [sz];
	for (int i = 0; i < sz; ++i)
	{
		data[i] = new double[sz];
		Di[i] = i + 1;
		Dj[i] = i + 1;
		for (int j = 0; j < sz; ++j)
			if (i == j) data[i][i] = D_MAX;
			else data[i][j] = dist(D[i].first, D[i].second, D[j].first, D[j].second);
	}
}

void matrix::clear()
{
	if (sz > 0)
	{
		for (int i = 0; i < sz; ++i)
			delete [] data[i];
		delete [] data;
		delete [] Di;
		delete [] Dj;
		data = NULL;
		sz = 0;
	}
}

bool matrix::is_empty() const
{
	return sz == 0;
}


void matrix::set(int eI, int eJ, double S)
{
	for (int i = 0; i < this->sz; ++i)
		for (int j = 0; j < this->sz; ++j)
			if (Di[i] == eI && Dj[j] == eJ)
				if (S == -1)
					this->data[i][j] = D_MAX;
				else
					this->data[i][j] = S;
}

void matrix::erase(int eI, int eJ)
{
	int *dDi, *dDj;
	double** B;
	dDi = new int[this->sz - 1];
	dDj = new int[this->sz - 1];
	B = new double* [this->sz - 1];
	this->set(eJ, eI);
	for (int i = 0; i < this->sz; ++i)
	{
		if (Di[i] != eI)
		{
			dDi[i - (Di[i] > eI)] = this->Di[i];
			B[i - (Di[i] > eI)] = new double [sz - 1];
			for (int j = 0; j < this->sz; ++j)
				if (Dj[j] != eJ)
					B[i - (Di[i] > eI)][j - (Dj[j] > eJ)] = this->data[i][j];
		}
		if (Dj[i] != eJ)
			dDj[i - (Dj[i] > eJ)] = this->Dj[i];
	}
}

bool operator< (const res& A,const res& B)
{
	if (A.H == B.H)
		return A.R.size() < B.R.size();
	return A.H > B.H;
}
