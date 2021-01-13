#if(0)
#ifndef __SCALARMOEA_H_
#define __SCALARMOEA_H_


#include "global.h"
#include "individual.h"
#include "common.h"
#include "reproduction.h"
#include "portfolio.h"


class TSUB {
public:
	TIndividual data;
	vector<double> weight;
	vector<int>    neighbor;
	TSUB();
	~TSUB();
	void operator=(const TSUB &sub2);
};


TSUB::TSUB()
{
	for (int n = 0; n < num; n++) weight.push_back(0.0);
}

TSUB::~TSUB()
{
}

void TSUB::operator=(const TSUB &sub2)
{
	data = sub2.data;
	weight = sub2.weight;
	neighbor = sub2.neighbor;
}


void init_weights(vector<TSUB>& population, int ini_seed) {
	popsize = 0;
	vector<int> bits;
	for (int i = 0; i < num; i++) bits.push_back(0);
	while (bits[num - 1] < ini_seed + 1) {
		bits[0]++;
		for (int j = 0; j < num - 1; j++)
			if (bits[j] == ini_seed + 1) {
				bits[j] = 0;
				bits[j + 1]++;
			}
		double sum = 0;
		for (int k = 0; k < num; k++) {
			sum += bits[k];
			//printf("%d ", bits[k]);
		}

		if (sum == ini_seed) {
			TSUB sub;
			for (int s = 0; s < num; s++)
			{
				sub.weight[s] = 1.0*bits[s] / ini_seed;
				//printf("%f ",sub.weight[s]);
			}
			//printf("\n");
			population.push_back(sub);
			popsize++;
		}
	}
}


void init_neighborhood(vector<TSUB> &population)
{
	double *x = new double[population.size()];
	int    *idx = new int[population.size()];
	for (int i = 0; i < popsize; i++) {
		for (int j = 0; j < popsize; j++) {
			x[j] = distanceVector(population[i].weight, population[j].weight);
			idx[j] = j;
		}
		minfastsort(x, idx, popsize, niche);
		for (int k = 0; k < niche; k++)
			population[i].neighbor.push_back(idx[k]);

	}
	delete[] x;
	delete[] idx;

}

void init_population(vector<TSUB> &population, CStock&stock)
{

	for (int i = 0; i < num; i++) {
		upp.push_back(-1.0e+30);
		low.push_back(1.0e+30);
	}

	vector<TIndividual> ini_popu;

	for (int n = 0; n < popsize; n++) {
		population[n].data.rnd_init();
		population[n].data.cardinality(cad);
		/*if (bTrans) {
			population[n].data.convert(stock);
			population[n].data.repair(stock);
		}*/
		population[n].data.obj_eval(stock);
		population[n].data.sharp_eval(stock);

	}

	for (int i = 0; i < num; i++) {
		double min = 1.0e+30;
		int    min_id;
		for (int j = 0; j < popsize; j++) {
			if (population[j].data.y_obj[i] < min) {
				min = population[j].data.y_obj[i];
				min_id = j;
			}
		}
		for (int k = 0; k < num; k++)
			simplex[i][k] = population[min_id].data.y_obj[k];
	}


	for (int i = 0; i < num; i++) {
		for (int k = 0; k < num; k++) {
			if (upp[k] < simplex[i][k]) upp[k] = simplex[i][k];
			if (low[k] > simplex[i][k]) low[k] = simplex[i][k];;
		}
	}
}

void update_reference(TIndividual child) {
	bool updated = false;
	for (int n = 0; n < num; n++) {
		if (low[n] > child.y_obj[n]) {
			updated = true;
			low[n] = child.y_obj[n];
			for (int j = 0; j < num; j++)simplex[n][j] = child.y_obj[j];
		}
	}
	if (updated) {
		for (int i = 0; i < num; i++) {
			for (int k = 0; k < num; k++) {
				if (upp[k] < simplex[i][k]) upp[k] = simplex[i][k];
				if (low[k] > simplex[i][k]) low[k] = simplex[i][k];;
			}
		}
	}
}

void update_neighbors(vector<TSUB> &population, int k, TIndividual child) {
	int count = 0;
	for (int n = 0; n < niche; n++) {
		int nb = population[k].neighbor[n];
		double f1, f2;
		f1 = scalarvalue(population[nb].data.y_obj, population[nb].weight);
		f2 = scalarvalue(child.y_obj, population[nb].weight);
		if (f2 < f1) {
			count++;
			population[nb].data = child;
		}
	}
}

void show_population(vector<TSUB> population) {
	for (int n = 0; n < popsize; n++)
		population[n].data.show_objective();
}

void save_population(vector<TSUB> population, char *filename) {
	fstream fout;
	fout.open(filename, std::ios::out);
	for (int n = 0; n < popsize; n++)
	{
		for (int k = 0; k < num; k++)
		{
			fout << population[n].data.y_obj[k] << "  ";
			//cout << population[n].data.y_obj[k] << "  ";
		}
		//cout << population[n].data.sharp_index;
		fout << population[n].data.sharp_index;
		fout << "\n";
		//cout << "\n";
	}
	fout.close();
}


void dmoea(int ini_seed, int max_gen, CStock &stock, int run, char *root)
{
	strcpy(methodtype, "TCHEBY2");
	vector <TSUB>  population;

	simplex = new double*[num];
	for (int i = 0; i < num; i++)
	{
		simplex[i] = new double[num];
	}


	init_weights(population, ini_seed - 1);
	init_neighborhood(population);
	init_population(population, stock);


	char filename[1024];

	for (int gen = 2; gen <= max_gen; gen++)
	{
		for (int k = 0; k < popsize; k++)
		{

			int p1, p2, p3, n1, n2, n3;
			p1 = int(rnd_uni(&rnd_uni_init)*niche);
			n1 = population[k].neighbor[p1];
			p2 = int(rnd_uni(&rnd_uni_init)*niche);
			n2 = population[k].neighbor[p2];
			p3 = int(rnd_uni(&rnd_uni_init)*niche);
			n3 = population[k].neighbor[p3];

			TIndividual child;
			diffevolution(population[k].data.x_var, population[n2].data.x_var, population[n3].data.x_var, child.x_var);
			realmutation(child, 0.01, 0, 1);
			child.cardinality(cad);
			/*if (bTrans) {
				child.convert(stock);
				child.repair(stock);
			}*/
			child.obj_eval(stock);
			child.sharp_eval(stock);
			update_reference(child);
			update_neighbors(population, k, child);
		}
		
		if ( gen % 100 == 0) {
			
			printf("MOEA/D ----> gen = %d\n", gen);

			sprintf(filename, "./pof1/MOEAD_%d_%d_%d.txt",dim, gen, run);
			save_population(population, filename);

		}
	}
	
	for (int i = 0; i < num; i++) delete[]simplex[i];
	delete[]simplex;

	population.clear();

	upp.clear(); low.clear();  // bug found: memory for global variables should be released.


}

#endif
#endif