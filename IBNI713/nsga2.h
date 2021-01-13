#if(0)
#ifndef __NSGA2_H_

#define __NSGA2_H_
#include "global.h"
#include "individual.h"
#include "common.h"
#include "reproduction.h"

int i, k;
void init_population(vector<TIndividual> &population, CStock &stock) {
	for (int n = 0; n < popsize; n++) {
		TIndividual ind;
		ind.rnd_init();
		ind.cardinality(cad);
		/*if (bTrans) {
			ind.convert(stock);
			ind.repair(stock);
		}*/
		ind.obj_eval(stock);
		ind.sharp_eval(stock);
		population.push_back(ind);
	}
}


void fast_popurank(vector<TIndividual> &population) {
	int size = population.size();
	int** cset;

	int* rank = new int[size];

	cset = new int*[size];

	for (int i = 0; i < size; i++)
		cset[i] = new int[size];

	for (i = 0; i < size; i++)
	{
		rank[i] = 0;
		population[i].rank = -1;
		population[i].count = 0;

	}

	for (int k = 0; k < size; k++)
		for (int j = 0; j < size; j++)
		{
			if (k != j)
			{
				if (population[j] < population[k] && !(population[j] == population[k])) rank[k]++;

				if (population[k] < population[j] && !(population[j] == population[k]))
				{
					population[k].count++;
					int m = population[k].count - 1;
					cset[k][m] = j;
				}
			}
		}

	/*
	for(k=0; k<size;k++)
	{
		population[k].Show();

		std::cout<<"rank = "<<rank[k]<<"\n";

		for(int j=0; j<population[k].count; j++)
		{
			int id = cset[k][j];
			population[id].Show();
		}
		std::cout<<"\n";
	}//*/

	int curr_rank = 0;
	while (1)
	{
		int stop_count = 0;
		int* rank2 = new int[size];
		for (int k = 0; k < size; k++)
			rank2[k] = rank[k];

		for (k = 0; k < size; k++)
		{
			if (population[k].rank == -1 && rank[k] == 0)
			{
				population[k].rank = curr_rank;
				for (int j = 0; j < population[k].count; j++)
				{
					int id = cset[k][j];
					rank2[id]--;
					stop_count++;
				}
			}
		}

		for (k = 0; k < size; k++)
			rank[k] = rank2[k];

		delete[] rank2;
		curr_rank++;
		if (stop_count == 0) break;
	}

	delete[] rank;

	for (i = 0; i < size; i++)
		delete cset[i];
	delete[] cset;

}

void slow_popurank(vector<TIndividual> &population)
{
	int size = population.size();
	for (int i = 0; i < size; i++)
		population[i].rank = -1;
	int rank = 0;
	while (1) {
		int count = 0;
		for (int i = 0; i < size; i++) {
			if (population[i].rank == -1) {
				population[i].rank = rank;
				for (int j = 0; j < size; j++) {
					if (population[j].rank == -1 && population[j] < population[i] && !(population[j] == population[i])) {
						population[i].rank = -1;
						break;
					}
				}
			}
			else
				count++;
		}
		if (count == size) break;
		rank++;
	}
}

void densitysharing(vector<TIndividual> &population, vector<TIndividual> &offspring) {
	population.clear();
	int size = offspring.size();
	int rank = 0;
	while (1) {
		int count = 0;
		for (int i = 0; i < size; i++)
			if (offspring[i].rank == rank)
				count++;

		int size2 = population.size() + count;
		if (size2 > popsize) {
			break;
		}

		for (i = 0; i < size; i++)
			if (offspring[i].rank == rank)
				population.push_back(offspring[i]);
		rank++;
		if (population.size() > popsize) break;
	}

	if (population.size() < popsize) {
		vector<TIndividual> list;
		for (int i = 0; i < size; i++)
			if (offspring[i].rank == rank)
				list.push_back(offspring[i]);
		int s2 = list.size();
		double *density = new double[s2];
		int    *idx = new int[s2];
		for (i = 0; i < s2; i++) {
			list[i].density = 0;
			for (int j = 0; j < num; j++) {
				double left = 1.0e+30, right = 1.0e+30;
				for (int k = 0; k < s2; k++) {
					if (list[k].y_obj[j] < list[i].y_obj[j] && list[i].y_obj[j] - list[k].y_obj[j] < left)  left = list[i].y_obj[j] - list[k].y_obj[j];
					if (list[k].y_obj[j] > list[i].y_obj[j] && list[k].y_obj[j] - list[i].y_obj[j] < right) right = list[k].y_obj[j] - list[i].y_obj[j];
				}
				list[i].density += left + right;
			}
			density[i] = -list[i].density;
			idx[i] = i;
		}

		int s3 = popsize - population.size();

		minfastsort(density, idx, s2, s3);

		for (i = 0; i < s3; i++)
			population.push_back(list[idx[i]]);

		delete[] density;
		delete[] idx;
	}
	offspring.clear();
}


void show_population(vector<TIndividual> &population) {
	for (int i = 0; i < population.size(); i++)
		population[i].show_objective();
}

void tour_selection(vector<TIndividual> &population, TIndividual &parent) {
	int p1, p2;
	p1 = int(rnd_uni(&rnd_uni_init)*popsize);
	p2 = int(rnd_uni(&rnd_uni_init)*popsize);
	if (population[p1].rank < population[p2].rank)
		parent = population[p1];
	else
		parent = population[p2];
}

void save_population(vector<TIndividual> population, char *filename) {
	fstream fout;
	fout.open(filename, std::ios::out);
	for (int n = 0; n < popsize; n++) {
		for (int k = 0; k < num; k++)
		{
			fout << population[n].y_obj[k] << "  ";
			//cout << population[n].y_obj[k] << "  ";
		}
		//cout << population[n].sharp_index;
		fout << population[n].sharp_index;
		fout << "\n";
		//cout << "\n";
	}
	fout.close();

}

void fill_offspring(vector<TIndividual> &offspring, TIndividual child) {
	bool flag = true;
	int  size = offspring.size();
	for (int i = 0; i < size; i++) {
		if (child == offspring[i]) {
			flag = false;
			break;
		}
	}
	if (flag) offspring.push_back(child);
}


void nsga2(int psize, int max_gen, CStock &stock, int run, char* root)
{

	char filename[1024];

	popsize = psize;;

	vector<TIndividual> population, offspring;

	init_population(population, stock);
	fast_popurank(population);
	sprintf(filename, "C:\\Users\\weijianle\\Desktop\\NBI\\INBI\\%s\\NSGA2_%d_%d_%d.txt", root, dim, cad, run);
	//save_population(population, filename);

	for (int gen = 2; gen <= max_gen; gen++) {
		for (int i = 0; i < popsize; i++) {
			TIndividual parent1, parent2, parent3, child;
			tour_selection(population, parent1);
			tour_selection(population, parent2);
			tour_selection(population, parent3);
			diffevolution(parent1.x_var, parent2.x_var, parent3.x_var, child.x_var);
			realmutation(child, 0.01, 0, 1);
			child.cardinality(cad);
			/*if (bTrans) {
				child.convert(stock);
				child.repair(stock);
			}*/
			child.obj_eval(stock);
			child.sharp_eval(stock);
			fill_offspring(offspring, child);
			fill_offspring(offspring, population[i]);
		}

		fast_popurank(offspring);
		densitysharing(population, offspring);

		/*if (gen % 100 == 0) {
			printf("NSGA-II ---> gen = %d\n", gen);
			sprintf(filename, "C:\\Users\\weijianle\\Desktop\\NBI\\INBI\\%s\\NSGA2_%d_%d_%d_%d.txt", root, dim, cad, gen, run);
			save_population(population, filename);
		}*/
		if (gen % 100 == 0) {

			printf("NSGA-II ---> gen = %d\n", gen);

			sprintf(filename, "./pof/nsga2_%d_%d_%d.txt",  dim,  gen, run);
			save_population(population, filename);

		}
	}
	
	
}

#endif
#endif