#ifndef __TINDIVIDUAL_H_
#define __TINDIVIDUAL_H_

#include "global.h"
#include "objective.h"
#include"portfolio.h"

class TIndividual {
public:
	TIndividual();
	virtual ~TIndividual();

	vector <double> x_var;   // the weight of stocks (including cost)
	vector <int>    n_var;   // the number of stocks   
	vector <double> y_obj;   // the objective values
	double	sharp_index;

	void   rnd_init();               // the initialization of portfolio
	void   obj_eval(CStock &stock);  // the evaluation of objective values
	void   var_norm();               // the normalization of weight
	void   sharp_eval(CStock &stock);//the evaluation of sharp values
	// for transcation cost 
	vector <double> cost;            // the transaction cost
	double R0;                       // the remainder of investment

	void   convert(CStock &stock);   // the conversion of an infeasible solution to integer constraints 
	void   repair(CStock &stock);    // the reinvestment of remainder

	bool   operator<(const TIndividual &ind2);
	bool   operator==(const TIndividual &ind2);
	void   operator=(const TIndividual &ind2);

	void   show_objective();
	void   show_variable();

	void   cardinality(int K);       // the random removal for cardinality constraints
	int    get_dimens();             // the dimensionality of weight (positive weights)

	double sumx, sumy;
	int    rank;
	int    count;

	double density;

	double best;
	double second;
	bool   boo;

};

TIndividual::TIndividual()  // construction function
{
	// allocation of memory for vectors
	for (int i = 0; i < dim; i++) {
		x_var.push_back(0.0);
		n_var.push_back(0);
		cost.push_back(0.0);
	}
	for (int n = 0; n < num; n++)
		y_obj.push_back(0.0);
	rank = 0;
}

TIndividual::~TIndividual() {

}

void TIndividual::rnd_init() {
	// assignment of random number in [0,1] for weight
	for (int n = 0; n < dim; n++)
		x_var[n] = rnd_uni(&rnd_uni_init);
}

void TIndividual::var_norm()
{
	// normalization of weight
	sumx = 0;
	//int flag,t;
	double s=-1.0,sumx_=0.0;
	for (int n = 0; n < dim; n++)
		sumx += x_var[n];	// the sum of weight
	for (int n = 0; n < dim; n++){
        //static int ttt=0;
        //cout<<"test"<<int(a*100)*0.01<<endl;

            x_var[n] = (x_var[n] /sumx);


	}



}

int TIndividual::get_dimens() {
	int dim2 = 0;
	for (int n = 0; n < dim; n++)
		if (x_var[n] != 0) dim2++;
	return dim2;
}

void TIndividual::cardinality(int K) {
	/*int dim2 = get_dimens();
	while (dim2 > K) {
		// randomly choose one stock for removal
		int id = int(rnd_uni(&rnd_uni_init)*dim);
		while (1) {
			if (x_var[id] > 0) {
				x_var[id] = 0;
				dim2--;
				break;
			}
			else {
				id++;
				if (id == dim) id = 0;
			}
		}
	}*/
	var_norm();
}
#if(0)
void TIndividual::convert(CStock &stock) {

	// memo: there is a bug in the repair method. 
	// Qingfu and I discussed this on 08/12/06
	// the variable x should contain the cost
	// the other change is the calculation of objective function

	double sum = 0;
	for (int n = 0; n < dim; n++) {
		double R1 = x_var[n] * stock.V0;
		if (cost_type == 1) {
			//case one: fixed cost		
			// x_var[n]*stock.V0: assume this is the money for investment
			// n_var[n]: the maximal number of stock n
			n_var[n] = int((R1 - fixed_cost) / stock.price[n]);
			// bug found here: nvar might be negative
			if (n_var[n] > 0) {
				cost[n] = fixed_cost;  // if equal or more than one, then pay cost				
			}
			else {
				n_var[n] = 0;
				cost[n] = 0;           // otherwise, don't pay cost
			}
		}

		/*if (cost_type == 2)   //case two: proportional cost
		{
			n_var[n] = int(R1 / (stock.price[n] * (1 + linear_rate)));
			if (n_var[n] > 0) {
				cost[n] = linear_rate * n_var[n] * stock.price[n];
			}
			else {
				n_var[n] = 0;
				cost[n] = 0;
			}
		}

		if (cost_type == 3) {
			int x, y, z;
			x = int((R1 - fixed_cost) / stock.price[n]);
			y = int(fixed_cost / (linear_rate*stock.price[n]));
			z = int(R1 / ((1 + linear_rate)*stock.price[n]));
			if (x > 0 && x <= y) {
				n_var[n] = x;
				cost[n] = fixed_cost;
			}
			else if (z > 0 && z >= y) {
				n_var[n] = z;
				cost[n] = linear_rate * n_var[n] * stock.price[n];
			}
			else {
				n_var[n] = 0;
				cost[n] = 0;
			}
		}

		if (cost_type == 4) {
			n_var[n] = int((R1 - fixed_cost) / (stock.price[n] * (1 + linear_rate)));
			if (n_var[n] > 0) {
				cost[n] = fixed_cost + linear_rate * n_var[n] * stock.price[n];
			}
			else {
				cost[n] = 0;
			}
		}*/

		x_var[n] = (n_var[n] * stock.price[n] + cost[n]) / stock.V0;

		sum += x_var[n];
	}
	R0 = stock.V0*(1 - sum);
	//std::cout<<R0<<endl; getchar();
}

void TIndividual::repair(CStock &stock) {
	// randomly pick an asset
	int p = int(dim*rnd_uni(&rnd_uni_init));
	while (n_var[p] <= 0) {
		p++;   if (p == dim) p = 0;
	}

	// all money to be invested
	double R1 = R0 + n_var[p] * stock.price[p] + cost[p];
	if (cost_type == 1) {	// fixed cost			
		n_var[p] = int((R1 - fixed_cost) / stock.price[p]);
		if (n_var[p] > 0) {
			cost[p] = fixed_cost;
		}
		else {
			n_var[p] = 0;
			cost[p] = 0;
		}
	}

	if (cost_type == 2) {  // proportional cost
		n_var[p] = int(R1 / (stock.price[p] * (1 + linear_rate)));
		if (n_var[p] > 0) {
			cost[p] = n_var[p] * linear_rate*stock.price[p];
		}
		else {
			n_var[p] = 0;
			cost[p] = 0;
		}
	}

	if (cost_type == 3) {
		int x, y, z;
		x = int((R1 - fixed_cost) / stock.price[p]);
		y = int(fixed_cost / (linear_rate*stock.price[p]));
		z = int(R1 / ((1 + linear_rate)*stock.price[p]));
		if (x > 0 && x <= y) {
			n_var[p] = x;
			cost[p] = fixed_cost;
		}
		else if (z > 0 && z >= y) {
			n_var[p] = z;
			cost[p] = linear_rate * n_var[p] * stock.price[p];
		}
		else {
			n_var[p] = 0;
			cost[p] = 0;
		}
	}

	if (cost_type == 4) {
		n_var[p] = int((R1 - fixed_cost) / (stock.price[p] * (1 + linear_rate)));
		if (n_var[p] > 0) {
			cost[p] = fixed_cost + linear_rate * n_var[p] * stock.price[p];
		}
		else {
			cost[p] = 0;
		}
	}

	x_var[p] = (n_var[p] * stock.price[p] + cost[p]) / stock.V0;

	R0 = R1 - n_var[p] * stock.price[p] - cost[p];
}
#endif
void TIndividual::obj_eval(CStock &stock) {
	if (!bTrans) {
		objectives1(x_var, y_obj, stock);      // objective functions in standard mean-variance 
	}
	/*else {
		double *ret = new double[dim];
		for (int i = 0; i < dim; i++)
			ret[i] = n_var[i] * stock.price[i] / stock.V0;
		objectives2(ret, y_obj, stock, cost); // objective functions in transaction cost
		delete[] ret;
	}*/
}

void TIndividual::sharp_eval(CStock &stock) {
	sharp_index = (-y_obj[0] - 0.0005) / y_obj[1];
}

void TIndividual::show_objective()
{
	for (int n = 0; n < num; n++) printf("%f ", y_obj[n]);
	printf(" \n");
}

void TIndividual::show_variable()
{
	for (int n = 0; n < dim; n++)
		printf("%f ", x_var[n]);
	printf("\n");
}

void TIndividual::operator=(const TIndividual &ind2)
{
	x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	n_var = ind2.n_var;
	cost = ind2.cost;
	sumx = ind2.sumx;
	sumy = ind2.sumy;
	rank = ind2.rank;
	count = ind2.count;
	density = ind2.density;
	best = ind2.best;
	second = ind2.second;
	boo = ind2.boo;
	sharp_index = ind2.sharp_index;
}

bool TIndividual::operator<(const TIndividual &ind2)
{
	bool dominated = true;
	for (int n = 0; n < num; n++)
	{
		if (ind2.y_obj[n] < y_obj[n]) return false;
	}
	if (ind2.y_obj == y_obj) return false;
	return dominated;
}


bool TIndividual::operator==(const TIndividual &ind2)
{
	if (ind2.y_obj == y_obj) return true;
	else return false;
}

#endif