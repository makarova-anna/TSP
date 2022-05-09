#include <iostream>
#include <iomanip>
#include <vector>
#include <queue>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

double D_MAX = 1000000000;

double dist(double Ax, double Ay, double Bx, double By)
{ 
    return sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
}

struct res;

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

    friend ostream& operator<< (ostream& T,const matrix &D);
};

struct res
{
	double H;
	vector <pair<int, int>> R;
	matrix M;
	res()
	{
		H = -1;
	}
	res(double H1, vector <pair<int, int>> R1, matrix mf)
	{
		H = H1;
		R = R1;
		M = mf;
	}
	res(double H1, matrix mf)
	{
		H = H1;
		M = mf;
	}
	res(matrix mf)
	{
		H = 0;
		M = mf;
	}
	friend bool operator< (const res& A, const res& B);
	res& operator= (const res& A)
	{
		this->H = A.H;
		this->R = A.R;
		this->M = A.M;
		return *this;
	}
	bool err()
	{
		return H == -1;
	}
};

bool operator< (const res& A,const res& B)
{
	if (A.H == B.H)
		return A.R.size() < B.R.size();
	return A.H > B.H;
}

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

unsigned matrix::size() const
{
	return this->sz;
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

ostream& operator<< (ostream& out,const matrix &D)
{
    for (int i = 0; i < D.size(); ++i)
		out << D.Dj[i] << ' ';
	out << '\n';
    for (int i = 0; i < D.size(); ++i)
    {
        for(int j = 0; j < D.size(); ++j)
			if (D[i][j] == D_MAX)
				out << " inf  ";
            else
				out << fixed << setprecision(5) << D[i][j] << ' ';
        out << ": " << D.Di[i] << '\n';
    }
    return out;
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
	for (int i = 0; i < this->sz; ++i)
		delete [] this->data[i];
	delete [] this->data;
	delete [] this->Di;
	delete [] this->Dj;
	this->data = B;
	this->Di = dDi;
	this->Dj = dDj;
	this->sz--;
}

double matrix::getMin(int I, bool rows)
{
	double min = D_MAX;
	for (int i = 0; i < this->sz; i++)
		if (rows && min > data[I][i])
			min = data[I][i];
		else if (min > data[i][I] && !rows)
			min = data[i][I];
	return min;
}

double matrix::getResultSum(vector<pair<int, int>> result)
{
	double sum = 0;
	for (auto I : result)
		sum += this->data[I.first - 1][I.second - 1];
	return sum;
}

double matrix::st1()
{
	double sm = 0;
	//cout << *this << '\n';

	for (int i = 0; i < this->sz; i++)
	{
		double min = this->getMin(i, true);
		if (min == D_MAX)
			return D_MAX;
		else if (min != 0)
			for (int j = 0; j < this->sz; j++)
				if(this->data[i][j] != D_MAX) this->data[i][j] -= min;
		sm += min;
		//cout << min << ' ';
	}
	//cout << '\n';
	for (int i = 0; i < this->sz; i++)
	{
		double min = this->getMin(i, false);
		if (min == D_MAX)
			return D_MAX;
		else if (min != 0)
			for (int j = 0; j < this->sz; j++)
				if (this->data[j][i] != D_MAX) this->data[j][i] -= min;
		sm += min;
		//cout << min << ' ';
	}
	//cout << '\n';
	return sm;

}

void Algoritnm(matrix M)
{
	matrix temp = M;
	double H = temp.st1();
	res A(H, temp);
	A = matrix::Run(A);
	if (A.err())
		cout << "Error\n";
	else
	{
		vector<pair<int, int>> result = A.R;
		cout << M.getResultSum(result) << "\n1 ";
		int i;
		for (i = 0; i < result.size(); ++i)
			if (result[i].first == 1)
				break;
		while (result[i].second != 1)
		{
			for (int j = 0; j < result.size(); ++j)
				if (result[i].second == result[j].first)
				{
					i = j;
					cout << result[i].first << ' ';
					break;
				}
		}
		cout << '\n';
	}
}

res matrix::Run(res A) 
{
	static priority_queue <res> Q;
	double H = A.H;
	vector<pair<int, int>> result = A.R;
	matrix M = A.M;

	if (M.size() == 1)
	{
		result.push_back(make_pair(M.Di[0], M.Dj[0]));
		return res(H, result, M);
	}

	M.st1();

	int Max = 0;
	for (int i = 0; i < M.size(); i++)
		for (int j = 0; j < M.size(); j++)
			if (M[i][j] == 0)
			{
				M[i][j] = D_MAX;
				double max = (M.getMin(i, true) == D_MAX || M.getMin(j, false) == D_MAX)? D_MAX: M.getMin(i, true) + M.getMin(j, false);
				if (max > Max) Max = max;
				M[i][j] = 0;
			}

	vector<pair<int, int>> Maxs;
	for (int i = 0; i < M.size(); i++)
		for (int j = 0; j < M.size(); j++)
			if (M[i][j] == 0)
			{
				M[i][j] = D_MAX;
				int max = (M.getMin(i, true) == D_MAX || M.getMin(j, false) == D_MAX)? D_MAX: M.getMin(i, true) + M.getMin(j, false);
				if (max == Max) Maxs.push_back(make_pair(M.Di[i], M.Dj[j]));
				M[i][j] = 0;
			}
	
	if (Maxs.size() == 0)
		return res();

	for (int i = 0; i < Maxs.size(); i++)
	{
		result.push_back(Maxs[i]);
		matrix temp(M);
		temp.erase(Maxs[i].first, Maxs[i].second);
		double H1 = temp.st1();
		Q.push(res(H + H1, result, temp));
		result.pop_back();
		if (Max != D_MAX)
		{
			matrix temp1(M);
			temp1.set(Maxs[i].first, Maxs[i].second);
			Q.push(res(H + Max, result, temp1));
		}
	}
	do
	{
		A = Q.top();
		Q.pop();
		A = matrix::Run(A);
		if (Q.empty())
			return A;
	}
	while (A.err());
	return A;
}

int main(int argc, char ** argv)
{
	ifstream F;
	int n;
	vector<pair<double, double>> X1;
	double** X2;
	matrix A;
	switch(argc)
	{
	case(3):
		if (argv[2][0] != 't')
		{
			cout << "Wrong Arguments\n";
			return 0;
		}
        F.open(argv[1]);
        if (!F.is_open())
        {
            cout << "Can not open file\n";
            return 1;
        }
        F >> n;
		X2 = new double* [n];
        for (int i = 0; i < n; ++i)
        {
			X2[i] = new double [n];
            for (int j = 0; j < n; ++j)
				F >> X2[i][j];
		} 
        F.close();
		A = matrix(X2, n);
        for (int i = 0; i < n; ++i)
			delete [] X2[i];
		delete [] X2;
		break;
	case (2):
		int n;
        F.open(argv[1]);
        if (!F.is_open())
        {
            cout << "Can not open file\n";
            return 1;
        }
        F >> n;
        for (int i = 0; i < n; ++i)
        {
            pair<double, double> N;
            F >> N.first >> N.second;
            X1.push_back(N);
        } 
        F.close();
		A = matrix(X1);
		break;
	case (1):
        cout << "Without File\n";
        cin >> n;
        for (int i = 0; i < n; ++i)
        {
            pair<double, double> N;
            cin >> N.first >> N.second;
            X1.push_back(N);
        }
		A = matrix(X1);
		break;
	default:
		cout << "Wrong Argument\n";
		return 0;
	}
	cout << "start\n";
	Algoritnm(A);
	cout << "finish\n";
	return 0;
}

