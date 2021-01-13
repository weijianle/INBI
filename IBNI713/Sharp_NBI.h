

#if(0)
#ifndef __SCALARMOEA_H_
#define __SCALARMOEA_H_


#include "global.h"
#include "individual.h"
#include "common.h"
#include "reproduction.h"
#include "portfolio.h"



class SHARP {
public:
	TIndividual data;
	vector<double> weight;
	vector<int>    neighbor;

	SHARP();
	~SHARP();
	void operator=(const SHARP &sub2);
};


SHARP::SHARP()
{
	for (int n = 0; n < num; n++) weight.push_back(0.0);
}

SHARP::~SHARP()
{
}

void SHARP::operator=(const SHARP &sub2)
{
	data = sub2.data;
	weight = sub2.weight;
	neighbor = sub2.neighbor;
}


void init_weights(vector<SHARP>& population, int ini_seed) {
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

			SHARP sub;

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


void init_neighborhood(vector<SHARP> &population)
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

void init_population(vector<SHARP> &population, CStock&stock)
{

	for (int i = 0; i < num; i++) {//��ʼ���˲ο���
		uppsl.push_back(-1.0e+30);
		lowsl.push_back(1.0e+30);
		uppsr.push_back(-1.0e+30);
		lowsr.push_back(1.0e+30);
	}

	vector<TIndividual> ini_popu;

	for (int n = 0; n < popsize; n++) {//��ʼ����Ⱥ����
		population[n].data.rnd_init();
		population[n].data.cardinality(cad);
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
		if(i==0)
			population[0].data = population[min_id].data;
		else
			population[popsize-1].data = population[min_id].data;

	}

	for (int i = 0; i < popsize; i++) {
		if (population[i].data.sharp_index > max_sharp) {
			max_sharp = population[i].data.sharp_index;
			max_sharpid = i;
			//cout << max_sharp<<endl;
		}
	}
	simplex[num][0] = population[max_sharpid].data.y_obj[0];
	simplex[num][1] = population[max_sharpid].data.y_obj[1];
	population[24].data = population[max_sharpid].data;
	max_sharpid = 25;
	//cout << simplex[num][0] << " " << simplex[num][1] << endl;

	for (int i = 0; i < num + 1; i += 2) {
		for (int k = 0; k < num; k++) {

			if (uppsl[k] < simplex[i][k]) uppsl[k] = simplex[i][k];
			if (lowsl[k] > simplex[i][k]) lowsl[k] = simplex[i][k];;
		}
	}

	for (int i = 2; i >= 1; i--) {
		for (int k = 0; k < num; k++) {
			if (uppsr[k] < simplex[i][k]) uppsr[k] = simplex[i][k];
			if (lowsr[k] > simplex[i][k]) lowsr[k] = simplex[i][k];
		}
	}


}
void update_reference(TIndividual child) {
	bool updated = false;
	//bool updated_r = false;
	for (int n = 0; n < num; n++) {
		if ((low[n] > child.y_obj[n])) {
			updated = true;
			low[n] = child.y_obj[n];
			for (int j = 0; j < num; j++)
				simplex[n][j] = child.y_obj[j];
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
double distance_calc(int index1, int index2, vector<SHARP>& population)
{
	double distance = 0;
	int i;
	for (i = index1; i < index2 - 7; i+=5)
	{
		distance += sqrt(pow(population[i].data.y_obj[0] - population[i + 5].data.y_obj[0], 2) + pow(population[i].data.y_obj[1] - population[i + 5].data.y_obj[1], 2));
	}
	distance += sqrt(pow(population[index2-1].data.y_obj[0] - population[i].data.y_obj[0], 2) + pow(population[index2-1].data.y_obj[1] - population[i].data.y_obj[1], 2));
	return distance;
}

void update_reference_sp(TIndividual child, vector<SHARP>& population, int gen) {
	if ((max_sharp < child.sharp_index))
	{

		max_sharp = child.sharp_index;
		population[max_sharpid - 1].data = child;

		for (int j = 0; j < num; j++)simplex[num][j] = child.y_obj[j];
		uppsl[0] = simplex[0][0] > simplex[2][0] ? simplex[0][0] : simplex[2][0];
		uppsl[1] = simplex[0][1] > simplex[2][1] ? simplex[0][1] : simplex[2][1];
		lowsl[0] = simplex[0][0] < simplex[2][0] ? simplex[0][0] : simplex[2][0];
		lowsl[1] = simplex[0][1] < simplex[2][1] ? simplex[0][1] : simplex[2][1];
		uppsr[0] = simplex[1][0] > simplex[2][0] ? simplex[1][0] : simplex[2][0];
		uppsr[1] = simplex[1][1] > simplex[2][1] ? simplex[1][1] : simplex[2][1];
		lowsr[0] = simplex[1][0] < simplex[2][0] ? simplex[1][0] : simplex[2][0];
		lowsr[1] = simplex[1][1] < simplex[2][1] ? simplex[1][1] : simplex[2][1];

	}
	if (lowsl[0] > child.y_obj[0])
	{

		population[0].data = child;
		simplex[0][0] = child.y_obj[0]; simplex[0][1] = child.y_obj[1];
		uppsl[0] = simplex[0][0] > simplex[2][0] ? simplex[0][0] : simplex[2][0];
		uppsl[1] = simplex[0][1] > simplex[2][1] ? simplex[0][1] : simplex[2][1];
		lowsl[0] = simplex[0][0] < simplex[2][0] ? simplex[0][0] : simplex[2][0];
		lowsl[1] = simplex[0][1] < simplex[2][1] ? simplex[0][1] : simplex[2][1];
	}
	else if (lowsr[1] > child.y_obj[1])
	{

		population[popsize - 1].data = child;

		simplex[1][0] = child.y_obj[0]; simplex[1][1] = child.y_obj[1];
		uppsr[0] = simplex[1][0] > simplex[2][0] ? simplex[1][0] : simplex[2][0];
		uppsr[1] = simplex[1][1] > simplex[2][1] ? simplex[1][1] : simplex[2][1];
		lowsr[0] = simplex[1][0] < simplex[2][0] ? simplex[1][0] : simplex[2][0];
		lowsr[1] = simplex[1][1] < simplex[2][1] ? simplex[1][1] : simplex[2][1];

	}
}
void update_weight(vector<SHARP>& population, int ini_seed)
{
	int popsize_0 = 0;
	vector<int> bits;
	for (int i = 0; i < num; i++) bits.push_back(0);
	while (bits[num - 1] < max_sharpid) {
		bits[0]++;
		for (int j = 0; j < num - 1; j++)
			if (bits[j] == max_sharpid) {
				bits[j] = 0;
				bits[j + 1]++;
			}
		double sum = 0;
		for (int k = 0; k < num; k++) {
			sum += bits[k];
			//printf("%d ", bits[k]);
		}
		if (sum == max_sharpid - 1) {
			SHARP sub;
			for (int s = 0; s < num; s++)
			{
				sub.weight[s] = 1.0*bits[s] / (max_sharpid - 1);
				//printf("%f ",sub.weight[s]);
			}
			//printf("\n");
			population[popsize_0].weight[0] = sub.weight[0];
			population[popsize_0].weight[1] = sub.weight[1];
			popsize_0++;
		}
	}
	vector<int> bitsr;
	for (int i = 0; i < num; i++) bitsr.push_back(0);
	while (bitsr[num - 1] < (ini_seed + 1 - max_sharpid)) {
		bitsr[0]++;
		for (int j = 0; j < num - 1; j++)
			if (bitsr[j] == (ini_seed + 1 - max_sharpid)) {
				bitsr[j] = 0;
				bitsr[j + 1]++;
			}
		double sum = 0;
		for (int k = 0; k < num; k++) {
			sum += bitsr[k];
			//printf("%d ", bits[k]);
		}
		if (sum == (ini_seed - max_sharpid)) {
			SHARP subr;
			for (int s = 0; s < num; s++)
			{
				subr.weight[s] = 1.0*bitsr[s] / (ini_seed - max_sharpid);
				//printf("%f ",subr.weight[s]);
			}
			//printf("\n");
			population[popsize_0].weight[0] = subr.weight[0];
			population[popsize_0].weight[1] = subr.weight[1];
			popsize_0++;
			//printf("%d", popsize);

		}
	}
}
void update_neighbors(vector<SHARP> &population, int k, TIndividual child) {
	int count = 0;
	for (int n = 0; n < niche; n++) {
		int nb = population[k].neighbor[n];
		double f1, f2;
		f1 = scalarvalue_sharp(population[nb].data.y_obj, population[nb].weight, nb);
		f2 = scalarvalue_sharp(child.y_obj, population[nb].weight, nb);
		if (f2 < f1) {
			count++;
			population[nb].data = child;

		}
	}
}

void show_population(vector<SHARP> population) {
	for (int n = 0; n < popsize; n++)
		population[n].data.show_objective();
}

void save_population(vector<SHARP> population, char *filename) {
	fstream fout;
	fout.open(filename, std::ios::out);
	for (int n = 0; n < popsize; n++)
	{
		if (n != max_sharpid)
		{
			for (int k = 0; k < num; k++)
			{
				fout << population[n].data.y_obj[k] << "  ";
				//cout << population[n].data.y_obj[k] << "  ";
			}
			//cout << population[n].data.sharp_index << "\n";
			fout << population[n].data.sharp_index << "\n";
		}

	}

	fout.close();
}


void dmoea_sharp(int ini_seed, int max_gen, CStock &stock, int run, char *root)
{
	int r = 5, tip =10;
	


	ini_seed++;
	strcpy(methodtype, "TCHEBY2");
	vector <SHARP>  population;

	simplex = new double*[num + 1];
	for (int i = 0; i < num + 1; i++)
	{
		simplex[i] = new double[num];
	}


	init_weights(population, ini_seed - 1);
	init_neighborhood(population);

	init_population(population, stock);
	update_weight(population, ini_seed - 1);


	char filename[1024];

	for (int gen = 2; gen <= max_gen; gen++)
	{
		int times = 0;
		int tag = 0;
		int temp = 0;


		for (int k = 0; k < popsize; k++)
		{

			if ((gen % tip == 0 && gen != max_gen)) {
				
				if (times < (3 * r ))
					tag = temp;
				else
					tag = int(rnd_uni(&rnd_uni_init)*(popsize - 1));
			}
			else
				tag = k;

			int p1, p2, p3, n1, n2, n3;
			p1 = int(rnd_uni(&rnd_uni_init)*niche);
			n1 = population[tag].neighbor[p1];
			p2 = int(rnd_uni(&rnd_uni_init)*niche);
			n2 = population[tag].neighbor[p2];
			p3 = int(rnd_uni(&rnd_uni_init)*niche);
			n3 = population[tag].neighbor[p3];

			TIndividual child;
			diffevolution(population[tag].data.x_var, population[n2].data.x_var, population[n3].data.x_var, child.x_var);
			realmutation(child, 0.01, 0, 1);
			child.cardinality(cad);

			child.obj_eval(stock);
			child.sharp_eval(stock);
			//cout << child.y_obj[0];

			update_reference_sp(child, population, gen);

			if (gen % 200 == 1&&k==0)
			{

				double dis1, dis2;
				
				dis1 = distance_calc(0, max_sharpid,population);
				dis2 = distance_calc(max_sharpid-1, popsize,population);
				
				max_sharpid = popsize*dis1 / (dis1 + dis2)+1;
				

				//sqrt(pow((simplex[1][0] - simplex[2][0]), 2.0) + pow(simplex[1][1] - simplex[2][1], 2.0))))+1;	

				update_weight(population, ini_seed - 1);



			}




			update_neighbors(population, tag, child);
			if ((gen % tip == 0 && gen != max_gen))
			{
				if (times < (r))
				{
					temp = 0;
					times++;
				}
				if (times < (r * 2) )
				{
					temp = max_sharpid - 1;
					times++;
				}
				if (times < (r * 3))

				{
					temp = popsize - 1;
					times++;
				}
			}

		}

		if ((gen) % 100 == 0) {
			//r--;
			//tip = 10 + 10*(gen / 100 + 1);
			//tip = 100-sqrt(100*100-pow(25*gen/200,2));
			printf("sharp-MOEA/D_sharp ----> gen = %d.maxsharp:", gen);
			cout << max_sharp << endl;
			sprintf(filename, "./pof/Sharp%d_%d_%d.txt", dim, gen, run);
			save_population(population, filename);

		}
	}

	for (int i = 0; i < num + 1; i++) {if(NULL!=simplex[i]){delete[]simplex[i];simplex[i]=NULL;}}
    if(simplex!=NULL){delete[]simplex;simplex=NULL;}

	population.clear();

	upp.clear(); low.clear();  // bug found: memory for global variables should be released.
}

#endif
#endif
