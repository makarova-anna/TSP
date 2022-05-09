#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <string>
#include <fstream>
using namespace std;

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

bool operator< (const res& A,const res& B)
{
	if (A.H == B.H)
		return A.R.size() < B.R.size();
	return A.H > B.H;
}



int main()
{
	int n;
	vector<pair<double, double>> X1;
	double** X2;
	matrix A;

    return 0;
}
