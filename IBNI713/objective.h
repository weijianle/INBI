#include "global.h"


#if(0)
// case 1: no transcation cost
void objectives0(vector<double> &x_var, vector<double>&y_obj, double **ret)
{
	int row = 1000;
	int col = dim;
	// calculate the first objective - mean return

	double sum0 = 0;
	vector<double> rs;
	for (int i = 0; i < row; i++) {
		double sum = 0;
		for (int j = 0; j < col; j++)
			sum += x_var[j] * ret[i][j];
		sum0 += sum;
		rs.push_back(sum);
	}

	double mean = sum0 / row;

	y_obj[0] = -mean;

	// calculate the second objective
	sum0 = 0;
	for (int k = 0; k < row; k++)
		sum0 += (rs[k] - mean)*(rs[k] - mean);
	y_obj[1] = sum0 / row;

}


// case 2: with transcation cost
void objectives2(double *x_var, vector<double>&y_obj, CStock &stock, vector <double> &cost)
{
	int row = 1001;
	int col = dim;
	int day = 365;
	double safe = 0.025;

	// the first objective - the effective return
	// R(x) = sum [x_i*(r_i - rs) - c_i/v0*(1+rs)]
	double sum = 0;
	for (int j = 0; j < col; j++)
		sum += x_var[j] * (day*stock.rate[j] - safe) - cost[j] * (1 + safe) / stock.V0;

	y_obj[0] = -(sum + safe);

	// calculate the second objective
	double sum0 = 0;
	for (int i = 0; i < col; i++)
		for (int j = 0; j < col; j++)
		{
			sum0 += x_var[i] * x_var[j] * stock.cmat[i][j];  // bug found: no addition used
		}
	y_obj[1] = day * sum0;

}
#endif
// case 3: standard model
void objectives1(vector<double> &x_var, vector<double>&y_obj, CStock &stock)
{
	int row = stock.NumStock;
	int col = dim;
	//int day = 365;
	// calculate the first objective - mean return

	double sum = 0;
	for (int j = 0; j < col; j++)
		sum += x_var[j] * stock.mean[j];

	y_obj[0] = -sum;


	// calculate the second objective
	double sum0 = 0;
	for (int i = 0; i < col; i++)
	{
		for (int j = 0; j < col; j++)
		{
			sum0 += x_var[i] * x_var[j] * stock.covar[i][j];  // bug found: no addition used

		}
	}
	y_obj[1] =  sum0;
}