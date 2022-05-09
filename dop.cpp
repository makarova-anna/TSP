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
