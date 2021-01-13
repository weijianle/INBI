#if(1)
#ifndef __SCALARMOEA_H_
#define __SCALARMOEA_H_

#include<random>
#include "global.h"
#include "individual.h"
#include "common.h"
#include "reproduction.h"
#include "portfolio.h"
#include<ctime>
#include<cassert>

using namespace std;
double sigmoid(double value);
double softplus(double value);
double relu(double value);
double leaky_relu(double value);
vector<double>m;
vector<double>v;

class SHARP {
public:
    TIndividual data;
    vector<double> weight;
    vector<int>    neighbor;
    vector<double> sigma;
    vector<double> chromosome;
    //vector<double> old_sigma;
    vector<double> old_chromosome;

    SHARP();
    ~SHARP();
    void ini_sigma();
    void ini_chromosome();
    void updata_sigma();
    void updata_chromosome();
    void repair();
    void operator=(const SHARP &sub2);
};
class NEAnet {
public:
    vector< vector<double> > b1;
    vector< double > b2;
    vector<double> b3;
    vector< vector<double> > w1;
    vector< double > w2;
    vector< vector<double> > w3;
    vector< bool >   hide;
    //vector<bool>   hide2;

    NEAnet();
    ~NEAnet();
    void decode(SHARP);
    void ini_net(SHARP);

    SHARP net(CStock,SHARP);
};




SHARP NEAnet::net(CStock stock, SHARP population) {

    int input = stock.NumStock;
    vector< vector<double> > inputlayer;

    inputlayer.push_back(m);
    inputlayer.push_back(v);
    //inputlayer.push_back(mv);
    vector< vector<double> > hiddenlayer;
    for (int i = 0; i < input; i++) {
        double temp=0;
        vector<double> arr;
        for (int j =0;j<hidesize;j++) {
            double l[2]={inputlayer[0][i],inputlayer[1][i]};
            double r[2]={w1[0][j],w1[1][j]};
            temp = l[0]*r[0]+l[1]*r[1]+b1[0][j];
            arr.push_back(leaky_relu(temp));
        }


        hiddenlayer.push_back(arr);
        //arr.clear();
    }

    vector< vector<double> > hiddenlayer1;
    for (int i = 0; i < input; i++) {
        double temp=0;
        vector<double> arr;
        for (int j =0;j<hidesize2;j++) {
            for(int m=0;m<hidesize;m++){
                temp+=hiddenlayer[i][m]*w3[m][j];
            }
            temp+=b3[j];
            arr.push_back(leaky_relu(temp));
        }


        hiddenlayer1.push_back(arr);
        //arr.clear();
    }

    vector<double>outputlayer;
    double sum=10000000;int flag=0;
    for (int i = 0; i < input; i++) {
        double temp=0;

        for (int j = 0; j < hidesize2; j++) {

            temp += w2[j] * hiddenlayer1[i][j];
            //cout<<w2[j]<<" ";
        }
        temp += b2[0];

        if(relu(temp)<exp(-10) ){
            flag++;
        }
        outputlayer.push_back((temp));
    }
    if(flag==input)
        outputlayer[0]=0.5;
    for (int i = 0; i < input; i++) {

        population.data.x_var[i] = (relu(outputlayer[i]));
        //cout<<population.data.x_var[i]<<" ";

    }
    //cout<<endl;

    return population;

}
NEAnet::NEAnet() {

}
NEAnet::~NEAnet() {

}
SHARP::SHARP()
{
    for (int n = 0; n < num; n++) weight.push_back(0.0);
}

SHARP::~SHARP()
{
}
void SHARP::repair()
{
    int numofpositive=0;
    int minnum = 3;
    for (int i = lenOfChromosome - hidesize; i < lenOfChromosome; i++) {
        if (chromosome[i] >= 0)
            numofpositive++;
    }
    if (numofpositive < minnum)
    {
        srand(time(0));
        int repair_num = minnum - numofpositive;
        while (repair_num > 0) {
            int k = rand() % hidesize;
            if (chromosome[lenOfChromosome-1-k] < 0) {
                chromosome[lenOfChromosome - 1 - k] = chromosome[lenOfChromosome - 1 - k] * -1;
                repair_num--;
            }
        }

    }



}
void NEAnet::ini_net(SHARP sub) {
    vector<double> w,ww,www;
    vector<double> b;

    for (int i = 0; i < lenOfChromosome; i++) {

        if(i < (2 * hidesize + hidesize)) {
            if (i < 2 * hidesize) {
                if (i<hidesize ) {
                    w.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(w);
                }else if (i<2*hidesize ){
                    ww.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(ww);
                }else{
                    www.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(www);
                }
            } else {
                b.push_back(sub.chromosome[i]);
                if (i == 2 * hidesize + hidesize - 1) {
                    b1.push_back(b);
                    b1.push_back(b);
                }
            }
        }
        else if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2))
        {
            if(i<(2* hidesize + hidesize) + (hidesize2*hidesize))
            {
                for(int k=0;k<hidesize;k++){
                    vector<double> w;
                    for(int j=0;j<hidesize2;j++)
                        w.push_back(sub.chromosome[i]);
                    w3.push_back(w);

                }
            }else{
                //cout<<"here+"<<i;
                //getchar();
                b3.push_back(sub.chromosome[i]);
                //cout<<b2[0]<<" ";getchar();

            }
        }else if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2)+ (1*hidesize2 + 1))
        {
            if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2)+hidesize2)
            {
                w2.push_back(sub.chromosome[i]);
            }else{
                //cout<<"here+"<<i;
                //getchar();
                b2.push_back(sub.chromosome[i]);
                //cout<<b2[0]<<" ";getchar();

            }
        }else {
            bool bl;
            if (sub.chromosome[i] >= 0)
                bl = 1;
            else
                bl = 0;
            hide.push_back(bl);
        }
    }
}
void NEAnet::decode(SHARP sub) {

    vector<double> w,ww,www;
    vector<double> b;

    for (int i = 0; i < lenOfChromosome; i++) {

        if(i < (2 * hidesize + hidesize)) {
            if (i < 2 * hidesize) {
                if (i<hidesize ) {
                    w.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(w);
                }else if (i<2*hidesize ){
                    ww.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(ww);
                }else{
                    www.push_back(sub.chromosome[i]);
                    if((i+1)%hidesize==0)
                        w1.push_back(www);
                }
            } else {
                b.push_back(sub.chromosome[i]);
                if (i == 2 * hidesize + hidesize - 1) {
                    b1.push_back(b);
                    b1.push_back(b);
                }
            }
        }
        else if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2))
        {
            if(i<(2* hidesize + hidesize) + (hidesize2*hidesize))
            {
                for(int k=0;k<hidesize;k++){
                    vector<double> w;
                    for(int j=0;j<hidesize2;j++)
                        w.push_back(sub.chromosome[i]);
                    w3.push_back(w);

                }
            }else{
                //cout<<"here+"<<i;
                //getchar();
                b3.push_back(sub.chromosome[i]);
                //cout<<b2[0]<<" ";getchar();

            }
        }else if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2)+ (1*hidesize2 + 1))
        {
            if(i<(2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2)+hidesize2)
            {
                w2.push_back(sub.chromosome[i]);
            }else{
                //cout<<"here+"<<i;
                //getchar();
                b2.push_back(sub.chromosome[i]);
                //cout<<b2[0]<<" ";getchar();

            }
        }else {
            bool bl;
            if (sub.chromosome[i] >= 0)
                bl = 1;
            else
                bl = 0;
            hide.push_back(bl);
        }
    }

}
double sigmoid(double value)
{
    return 1.0 / (1 + exp(-1.0 * value));
}

double relu(double value)
{
    return value>=0?value:0;
}
double leaky_relu(double value)
{
    return value>=0?value:0.1*value;
}
double softplus(double value){
    return log(1.0+1.0/(1+exp(value*-1)));
}
void SHARP::operator=(const SHARP &sub2)
{
    data = sub2.data;
    weight = sub2.weight;
    neighbor = sub2.neighbor;
    sigma = sub2.sigma;
    chromosome =sub2. chromosome;
    old_chromosome =sub2. old_chromosome;
}

void SHARP::ini_sigma() {
    lenOfChromosome;
    double sigm;
    sigm = sqrt(1.0 / (3 * lenOfChromosome));
    for (int i=0; i < lenOfChromosome; i++) {
        sigma.push_back(sigm);
    }

}
default_random_engine  e;
uniform_real_distribution<double> u(0,1);
void SHARP::ini_chromosome() {
    double chrom;;
    double  r;

    for (int i=0; i < lenOfChromosome; i++) {
        chrom=u(e);

        chromosome.push_back(chrom);
    }
}
void SHARP::updata_sigma() {

    double eta0 = 1.0 / sqrt(2 * lenOfChromosome);
    double eta1 = 1.0 / sqrt(2 * sqrt(lenOfChromosome));
    double rand0 = gaussrand();

    for (int i = 0; i < lenOfChromosome;i++) {
        double oldsigm = sigma[i];
        sigma[i] = sigma[i] * exp(eta0*rand0+eta1*gaussrand());
        if (sigma[i] < exp(-10))
            sigma[i] = oldsigm;
    }

}
void SHARP::updata_chromosome() {
    double r;
    if(rnd_uni(&rnd_uni_init) <= 0.05){
        for (int i = 0; i < lenOfChromosome; i++)
        {
            if(rnd_uni(&rnd_uni_init) <= 0.5){
                //do{
                    r=gaussrand();
                //}while(r<-1||r>1);
                chromosome[i] += r*sigma[i];
            }
        }
    }
    else{
        a=2;
        while(a>0.1||a<-0.1)a=gaussrand();
        for (int i = 0; i < lenOfChromosome; i++)
        {
            r=gaussrand();
            chromosome[i]=chromosome[i]+(chromosome[i]-old_chromosome[i])*a;
            if(rnd_uni(&rnd_uni_init) <= 0.1){
                chromosome[i] += sigma[i]*r;
            }
        }
    }


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
        }
        if (sum == ini_seed) {
            SHARP sub;
            for (int s = 0; s < num; s++)
            {
                sub.weight[s] = 1.0*bits[s] / ini_seed;
            }
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

void init_population(vector<SHARP> &population, CStock&stock,vector<NEAnet> NET)
{

    for (int i = 0; i < num; i++) {//初始化端参考点
        uppsl.push_back(-1.0e+30);
        lowsl.push_back(1.0e+30);
        uppsr.push_back(-1.0e+30);
        lowsr.push_back(1.0e+30);
    }

    vector<TIndividual> ini_popu;

    for (int n = 0; n < popsize; n++) {//初始化种群个体
        population[n]=NET[n].net(stock, population[n]);
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

void update_reference_sp(SHARP child, vector<SHARP>& population, int gen) {
    if ((max_sharp < child.data.sharp_index))
    {

        max_sharp = child.data.sharp_index;
        population[max_sharpid - 1].data = child.data;

            population[max_sharpid - 1].sigma = child.sigma;
            population[max_sharpid - 1].chromosome = child.chromosome;
            population[max_sharpid - 1].old_chromosome = child.old_chromosome;



        for (int j = 0; j < num; j++)simplex[num][j] = child.data.y_obj[j];
        uppsl[0] = simplex[0][0] > simplex[2][0] ? simplex[0][0] : simplex[2][0];
        uppsl[1] = simplex[0][1] > simplex[2][1] ? simplex[0][1] : simplex[2][1];
        lowsl[0] = simplex[0][0] < simplex[2][0] ? simplex[0][0] : simplex[2][0];
        lowsl[1] = simplex[0][1] < simplex[2][1] ? simplex[0][1] : simplex[2][1];
        uppsr[0] = simplex[1][0] > simplex[2][0] ? simplex[1][0] : simplex[2][0];
        uppsr[1] = simplex[1][1] > simplex[2][1] ? simplex[1][1] : simplex[2][1];
        lowsr[0] = simplex[1][0] < simplex[2][0] ? simplex[1][0] : simplex[2][0];
        lowsr[1] = simplex[1][1] < simplex[2][1] ? simplex[1][1] : simplex[2][1];

    }
    if (lowsl[0] > child.data.y_obj[0])
    {

        population[0].data = child.data;

            population[max_sharpid - 1].sigma = child.sigma;
            population[max_sharpid - 1].chromosome = child.chromosome;
            population[max_sharpid - 1].old_chromosome = child.old_chromosome;

        simplex[0][0] = child.data.y_obj[0]; simplex[0][1] = child.data.y_obj[1];
        uppsl[0] = simplex[0][0] > simplex[2][0] ? simplex[0][0] : simplex[2][0];
        uppsl[1] = simplex[0][1] > simplex[2][1] ? simplex[0][1] : simplex[2][1];
        lowsl[0] = simplex[0][0] < simplex[2][0] ? simplex[0][0] : simplex[2][0];
        lowsl[1] = simplex[0][1] < simplex[2][1] ? simplex[0][1] : simplex[2][1];
    }
    else if (lowsr[1] > child.data.y_obj[1])
    {

        population[popsize - 1].data = child.data;

            population[max_sharpid - 1].sigma = child.sigma;
            population[max_sharpid - 1].chromosome = child.chromosome;
            population[max_sharpid - 1].old_chromosome = child.old_chromosome;


        simplex[1][0] = child.data.y_obj[0]; simplex[1][1] = child.data.y_obj[1];
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
void update_neighbors(vector<SHARP> &population, int k, SHARP child,int gen) {
    int count = 0;
    for (int n = 0; n < niche; n++) {
        int nb = population[k].neighbor[n];
        double f1, f2;
        f1 = scalarvalue_sharp(population[nb].data.y_obj, population[nb].weight, nb);
        f2 = scalarvalue_sharp(child.data.y_obj, population[nb].weight, nb);
        if (f2 < f1) {
            count++;
            population[nb].data = child.data;

                population[nb].sigma = child.sigma;
                population[nb].old_chromosome=population[nb].chromosome;
                population[nb].chromosome = child.chromosome;



            updatetimes++;


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

            /*for(int k=0;k<lenOfChromosome;k++)
                printf("%4.3lf ",population[n].chromosome[k]);
            cout<<endl;*/
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
SHARP proliferation(CStock stock,SHARP population) {

    SHARP child = population;

    child.updata_sigma();

    child.updata_chromosome();

    child.repair();
    NEAnet NETproliferation;

    NETproliferation.decode(child);
    child=NETproliferation.net(stock, child);
    return child;
}

void dmoea_sharpnea(int ini_seed, int max_gen, CStock &stock, int run, char *root)
{

    for (int i = 0; i < dim; i++) {
        m.push_back((stock.mean[i]));
        v.push_back((stock.var[i]));
        //mv.push_back((stock.mean[i])/(stock.var[i]));
    }

    int r = 17, tip =10;

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
    vector<NEAnet> NET(ini_seed);
    for (int i = 0; i < ini_seed; i++)
    {
        population[i].ini_sigma();
        population[i].ini_chromosome();
        population[i].old_chromosome=population[i].chromosome;
        NET[i].ini_net(population[i]);

    }

    init_population(population, stock, NET);

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
                {
                    tag = temp;
                }
                else
                    tag = int(rnd_uni(&rnd_uni_init)*popsize);
            }
            else
                tag = k;
            int p1,  n1,p2,  n2,p3,  n3;
            p1 = int(rnd_uni(&rnd_uni_init)*niche);
            n1 = population[tag].neighbor[p1];
            p2 = int(rnd_uni(&rnd_uni_init)*niche);
            n2 = population[tag].neighbor[p2];
            p3 = int(rnd_uni(&rnd_uni_init)*niche);
            n3 = population[tag].neighbor[p3];




            SHARP child;
            child=population[tag];
            if(gen>600){
                diffevolution(population[tag].data.x_var, population[n2].data.x_var, population[n3].data.x_var, child.data.x_var);
                realmutation(child.data, 0.01, 0, 1);
            }
            else
            //diffevolution(population[tag].chromosome, population[n2].chromosome, population[n3].chromosome, child.chromosome);
            {
                child=proliferation(stock,child);
            }




            child.data.cardinality(cad);

            child.data.obj_eval(stock);
            child.data.sharp_eval(stock);

            update_reference_sp(child, population, gen);

            if (gen % 200 == 1&&k==0)
            {

                double dis1, dis2;

                dis1 = distance_calc(0, max_sharpid,population);
                dis2 = distance_calc(max_sharpid-1, popsize,population);

                max_sharpid = popsize*dis1 / (dis1 + dis2)+1;
                update_weight(population, ini_seed - 1);
            }

            update_neighbors(population, tag, child,gen);
            if ((gen % tip == 0 && gen != max_gen))
            {
                if (times < (1))
                {
                    temp = 0;
                    times++;
                }
                if (times < (26) )
                {
                    temp = max_sharpid - 1;
                    times++;
                }
                if (times < (51))

                {
                    temp = popsize - 1;
                    times++;
                }
            }

        }

        if ((gen) % 100 == 0) {
            //show_population(population);
            delay-=0.01;
            //tip = 10 + 10*(gen / 100 + 1);
            //tip = 100-sqrt(100*100-pow(25*gen/200,2));
            //if(gen==10000){
            printf("sharpnea-MOEA/D_sharp ----> gen = %d.maxsharp:", gen);
            cout << max_sharp ;
            cout<<"  更新次数："<<updatetimes<<endl;
            //}
            sprintf(filename, "./pof/SharpNEA%d_%d_%d.txt", dim, gen, run);
            save_population(population, filename);

        }
    }

    for (int i = 0; i < num + 1; i++) delete[]simplex[i];
    delete[]simplex;
    m.clear();
    v.clear();
    //mv.clear();
    population.clear();
    NET.clear();
    upp.clear(); low.clear();  // bug found: memory for global variables should be released.
}

#endif
#endif
