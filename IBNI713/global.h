#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
//#include "engine.h"
//#include <windows.h>
#define _CRT_SECURE_NO_WARNINGS

using namespace std;



vector <double> low, upp ,uppsl ,lowsl ,uppsr,lowsr;//

int
        dim,						   // dimensionality
num,                           // number of objectives
popsize,                       // population size
niche,                         // neighborhood size
cad,                           // maximal number of stocks in portfolio

cost_type;                     // type of transaction cost
int    id_cx = 20,                    // indexes used in mutation and crossover
id_mu = 20,
        max_sharpid,
        hidesize,
        hidesize1,
        hidesize2,

        lenOfChromosome;

double a=2;
int    seed = 448;                   // seed for random number
int    updatetimes;
long   rnd_uni_init;

char   methodtype[1024] = "TCHEBY";   // type of decompostion

bool   bTrans,                       // decide if transaction cost is considered=
flagsp = false;
double fixed_cost,                    // fixed transaction cost
linear_rate,                   // linear transaction cost                  
init_wealth,           // initial wealth
max_sharp,
        delay=0.5,
        HV;

double xHypervolumeRefPoint_TopRight = 0.005;
double yHypervolumeRefPoint_TopRight = 0;

double **simplex;
//double simplex[10][10];


int    loop;

#include "random.h"

#endif