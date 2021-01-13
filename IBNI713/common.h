#ifndef __COMMON_H_
#define __COMMON_H_

#include "global.h"

int n;
void minfastsort(double x[], int idx[], int n, int m)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = i + 1; j < n; j++)
			if (x[i] > x[j])
			{
				double temp = x[i];
				x[i] = x[j];
				x[j] = temp;
				int id = idx[i];
				idx[i] = idx[j];
				idx[j] = id;
			}
	}
}

double distanceArray(double vec1[], double vec2[], int dim)
{
	double sum = 0;
	for (int n = 0; n < dim; n++)
		sum += (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double distanceVector(vector <double> &vec1, vector <double> &vec2)
{
	double sum = 0;
	for (int n = 0; n < num; n++)
		sum += (vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double norm_objective(vector <double> &y)
{
	double sum = 0;
	for (int i = 0; i < num; i++)
		sum = sum + y[i] * y[i];
	return sqrt(sum);
}

double norm_decision(double x[])
{
	double sum = 0;
	for (int i = 0; i < dim; i++)
		sum = sum + x[i] * x[i];
	return sqrt(sum);
}

double norm_sum(vector<double>&vec)
{
	double sum = 0;
	for (int i = 0; i < num; i++)
		sum = sum + vec[i];
	return sum;
}

double innerproduct(vector <double>&vec1, vector <double>&vec2)
{
	double sum = 0;
	for (int i = 0; i < num; i++)
		sum += vec1[i] * vec2[i];
	return sum;
}


double scalarvalue(vector<double> y_obj, vector<double> weight)
{
	double fval;
	// Tchebycheff approach
	if (!strcmp(methodtype, "TCHEBY")) {
		fval = -1.0e+30;
		for (int n = 0; n < num; n++) {
			double temp = weight[n] * (y_obj[n] - low[n]) / (upp[n] - low[n]);
			if (temp > fval) fval = temp;
		}
	}

	if (!strcmp(methodtype, "TCHEBY2")) {

		vector<double> ref;
		double sum = 0;
		for (int n = 0; n < num; n++) {
			sum += upp[n] - low[n];
			ref.push_back(0.0);
			for (int j = 0; j < num; j++)
				ref[n] += weight[j] * simplex[j][n];
		} 
		
		/*cout << "upp: " << upp[0] << " " << upp[1]<<"  ";
		cout << "low: " << low[0] << " " << low[1]<<"  ";
		cout << "simplex[0]: " << simplex[0][0] << " " << simplex[0][1]<< "  ";
		cout << "simplex[1]: " << simplex[1][0] << " " << simplex[1][1]<< "  ";
		cout << "simplex[2]: " << simplex[2][0] << " " << simplex[2][1]<<endl;*/
		
		vector<double> vec; 
		for (n = 0; n < num; n++)
			vec.push_back(sum - (upp[n] - low[n]));

		fval = -1.0e+30;
		for (n = 0; n < num; n++) {
			double temp = vec[n] * (y_obj[n] - ref[n]);
			if (temp > fval) fval = temp;
		}
	}


	// Boundary intersection
	if (!strcmp(methodtype, "BOUNDARY")) {

		double nw = norm_objective(weight);
		vector<double> diff, diff2;
		for (int n = 0; n < num; n++) {
			diff.push_back((y_obj[n] - low[n]) / (upp[n] - low[n]));
		}

		double d1, d2;
		d1 = innerproduct(diff, weight) / nw;
		for (n = 0; n < num; n++)
			diff2.push_back(y_obj[n] - low[n] - d1 * weight[n] / nw);
		d2 = norm_objective(diff2);
		fval = d1 + 1 * d2;
	}


	// Boundary intersection
	if (!strcmp(methodtype, "NBILIKE")) {
		vector<double> ref;
		double sum = 0;
		for (int n = 0; n < num; n++) {
			sum += upp[n] - low[n];
			ref.push_back(0.0);
			for (int j = 0; j < num; j++)
				ref[n] += weight[j] * simplex[j][n];
		}

		vector<double> vec;
		for (n = 0; n < num; n++)
			vec.push_back(sum - (upp[n] - low[n]));

		vector<double> dif1, dif2;
		for (n = 0; n < num; n++)
			dif1.push_back(y_obj[n] - ref[n]);

		double nv = norm_objective(vec);
		double d1, d2;
		d1 = innerproduct(dif1, vec) / nv;

		for (n = 0; n < num; n++)
			dif2.push_back(y_obj[n] - (ref[n] + d1 * vec[n] / nv));
		d2 = norm_objective(dif2);

		fval = d1 + 0.4*d2;
	}

	if (!strcmp(methodtype, "LINEAR")) {
		fval = 0;
		for (int n = 0; n < num; n++) {
			fval += weight[n] * (y_obj[n] - low[n]) / (upp[n] - low[n]);
		}
	}

	return fval;

}

double scalarvalue_sharp(vector<double> y_obj, vector<double> weight, int nb)
{
	double fval;
	if (!strcmp(methodtype, "TCHEBY2")) {

		vector<double> ref;
		double sum = 0;
		if (nb < max_sharpid)
		{
			for (int n = 0; n < num; n++) {
				sum += uppsl[n] - lowsl[n];
				ref.push_back(0.0);
				for (int j = 0,k=0; j < num+1; j+=2,k++)
					ref[n] += weight[k] * simplex[j][n];
			}

			vector<double> vec;
			for (n = 0; n < num; n++)
				vec.push_back(sum - (uppsl[n] - lowsl[n]));

			fval = -1.0e+30;
			for (n = 0; n < num; n++) {
				double temp = vec[n] * (y_obj[n] - ref[n]);
				if (temp > fval) fval = temp;
			}
		}
		else {
			for (int n = 0; n <num ; n++ ) {
				sum += uppsr[n] - lowsr[n];
				ref.push_back(0.0);
				for (int j = 2, k = 0; j >0 ; j --, k++)
					ref[n] += weight[k] * simplex[j][n];
			}
			vector<double> vec;
			for (n = 0; n < num ; n ++ )
				vec.push_back(sum - (uppsr[n] - lowsr[n]));

			fval = -1.0e+30;
			for (n = 0; n < num; n++) {
				double temp = vec[n] * (y_obj[n] - ref[n]);
				if (temp > fval) fval = temp;
			}
		}
		
	}

	return fval;

}

/*template <class T>//MOEAD--NBI+INBI
double hv2(vector <T> &paretofront)
// returns the hypervolume of ps[0 ..] in 2D 
// assumes that ps is sorted improving
{

	double volume = 0;
	int popsize = 51;//NBI使用时为50
	if (paretofront.size() > 0)
		volume = ((yHypervolumeRefPoint_TopRight - paretofront[popsize-1].data.y_obj[0]) *
		(xHypervolumeRefPoint_TopRight - paretofront[popsize-1].data.y_obj[1]));
	for (int i = popsize-2; i >0; i--) // && paretofront[i][0]<xHypervolumeRefPoint_TopRight
		volume += ((xHypervolumeRefPoint_TopRight - paretofront[i].data.y_obj[1]) *
		(paretofront[i+1].data.y_obj[0] - paretofront[i].data.y_obj[0]));

	return volume;

}
/*
template <class T>//NSGA2

double hv2(vector <T> &paretofront)
// returns the hypervolume of ps[0 ..] in 2D 
// assumes that ps is sorted improving
{
	double volume = 0;
	int popsize = 50;
	TIndividual temp;
	for (int i = 0; i < popsize; i++)
	{
		for (int j = i + 1; j < popsize; j++)
			if (paretofront[i].y_obj[0] > paretofront[j].y_obj[0])
			{
				temp = paretofront[i];
				paretofront[i] = paretofront[j];
				paretofront[j] = temp;
			}
	}
	cout << paretofront[0].y_obj[0] << endl;
	if (paretofront.size() > 0)
		volume = ((yHypervolumeRefPoint_TopRight - paretofront[popsize - 1].y_obj[0]) *
		(xHypervolumeRefPoint_TopRight - paretofront[popsize - 1].y_obj[1]));
	for (int i = popsize - 2; i > 0; i--) // && paretofront[i][0]<xHypervolumeRefPoint_TopRight
		volume += ((xHypervolumeRefPoint_TopRight - paretofront[i].y_obj[1]) *
		(paretofront[i + 1].y_obj[0] - paretofront[i].y_obj[0]));

	return volume;

}*/
#endif