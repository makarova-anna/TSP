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

bool operator< (const res& A,const res& B)
{
	if (A.H == B.H)
		return A.R.size() < B.R.size();
	return A.H > B.H;
}

struct result
{
	double H;
	vector <pair<int, int>> R;
	matrix M;
	res(double H1, vector <pair<int, int>> R1, matrix mf)
	{
		H = H1;
		R = R1;
		M = mf;
	}


	friend bool operator< (const res& A, const res& B);

};

int main()
{

    return 0;
}
