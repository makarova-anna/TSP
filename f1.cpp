#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double dist(double Ax, double Ay, double Bx, double By)
{
    return sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
}

class matrix
{
private:
	double** data;
	int* Di;
	int* Dj;
	int D = 10000000;

	void set_M(vector<pair<double, double>> D);
    void set_M(double** D, int n, int* di = NULL, int* dj = NULL);

	double getMin(int I, bool rows = true);
public:
	matrix();
	matrix(const matrix &A);
		matrix(vector<pair<double, double>> D);
    matrix(double** D, int n, int* di = NULL, int* dj = NULL);
	~matrix();

};
int main()
{
    cout << "Hello world!" << endl;
    return 0;
}
