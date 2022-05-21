#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;


double dist(double Ax, double Ay, double Bx, double By)
{ 
    return sqrt(pow(Ax - Bx, 2) + pow(Ay - By, 2));
}


double** toM(vector<pair<double, double>> D)
{
    double** Ans;
    int n = D.size();
	Ans = new double* [n];
	for (int i = 0; i < n; ++i)
	{
		Ans[i] = new double[n];
		for (int j = 0; j < n; ++j)
			if (i == j) Ans[i][i] = -1;
			else Ans[i][j] = dist(D[i].first, D[i].second, D[j].first, D[j].second);
	}
    return Ans;
}

void clear(double** M, int n)
{
    for (int i = 0; i < n; ++i)
        delete [] M[i];
    delete [] M;
}

int Min1(double** M, int n, int I)
{
    int ans = 0;
    for (int i = 1; i < n; ++i)
        if (M[I][i] != -1)
            if (M[I][ans] > M[I][i] || M[I][ans] == -1)
                ans = i;
    if (M[I][ans] == -1)
    {
        cout << "FAlSE MIN";
        exit(1);
    }
    return ans;
}

int Min(double** M, int n, int I, vector <bool> IN)
{
    int ans = 0;
    for (int i = 0; i < n; ++i)
    {
        if (IN[ans] && !IN[i])
            ans = i;
        else if (!IN[i] && M[I][i] != -1)
            if (M[I][ans] > M[I][i] || M[I][ans] == -1)
                ans = i;
    }
    if (M[I][ans] == -1 || IN[ans])
    {
        cout << "FAlSE MIN";
        exit(1);
    }
    return ans;
}

bool In(vector<bool> IN)
{
    for(bool i : IN)
        if (!i)
            return true;
    return false;
}

vector <int> Algorithm(double** M, int n)
{
    int L = 0;
    vector <int> Ans;
    vector <bool> IN (n);
    for (int i = 0; i < n; ++i)
    {
        if (M[L][Min1(M, n, L)] > M[i][Min1(M, n, i)])
            L = i;
        IN[i] = false;
    }
    Ans.push_back(L);
    IN[L] = true;
    while(In(IN))
    {
        L = Min(M, n, L, IN);
        IN[L] = true;
        Ans.push_back(L);
    }
    return Ans;
}

int main(int argc, char ** argv)
{
    fstream F;
        int n;
        vector<pair<double, double>> X;
        double** M;
        switch(argc)
        {
        case (2):
            F.open(argv[1], std::fstream::in);
            if (!F.is_open())
                throw 1;
            F >> n;
            for (int i = 0; i < n; ++i)
            {
                pair<double, double> N;
                F >> N.first >> N.second;
                X.push_back(N);
            } 
            F.close();
            break;
        case (1):
            cout << "Without File\n";
            cin >> n;
            for (int i = 0; i < n; ++i)
            {
                pair<double, double> N;
                cin >> N.first >> N.second;
                X.push_back(N);
            }
            break;
        default:
            throw 2;
        }
        M = toM(X);
        vector <int> Ans = Algorithm(M, n);
        double Sz = 0;
        int L = -1;
        for (int I : Ans)
        {
            if (L != -1)
                Sz += M[L][I];
            L = I;
        }
        Sz += M[*(Ans.end() - 1)][*(Ans.begin())];
        cout << Sz << '\n';
        for (int I : Ans)
            cout << I << ' ';
        clear(M, n);
        return 0;
}
