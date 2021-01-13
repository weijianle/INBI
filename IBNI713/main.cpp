// ****************************************************************************
// MOEA/D: A multi-objective evolutionary algorithm based on decomposition    *
// for portfolio optimization                                                 *
// Author: HUI LI                                                             *
// Department of Computer Science, University of Essex                        *
// Date: November 2006                                                        *
// Copyright @ 2006 by HUI LI                                                 *
// ****************************************************************************
#include<iostream>
//#include "common.h"               // implement the basic operations in mathematics
#include "portfolio.h"            // define the class for portfolio data
#include "dmoea.h"           // implement the MOEA/D
#include"nsga2.h"				// implement the nsga2
#include"Sharp_NBI.h"
#include"SharpNEA.h"
//#include "nsga2.h"                // implement the NSGA-II

int main()
{

	int gen, pop, max_run;

	num = 2;         // the number of objectives--Ŀ�����

	char root[100];   sprintf(root, "pof");//ǰ���ļ���

	init_wealth = 1.0e+05;     // initial wealth for investment--Ͷ�ʼ�ֵ

	//dim = 50;						// number of stocks--��Ʊ��
	//cad = dim;					 // cardinality constraints--����Լ��

	bTrans = false;				 //false;
	cost_type = 1;
	linear_rate = 0.01;			 // proportion cost --����������--1%
	fixed_cost = 50;			  // fixed cost--�̶�������

	max_run = 30;				//������д���
	gen = 2000;			//��������
	pop = 50;					//��Ⱥ��С
	niche = pop/2;

	seed = (seed + 131) % 177;
	rnd_uni_init = -(long)seed;

	char dataset[5][20] = { "./port1.txt", "./port2.txt", "./port3.txt", "./port4.txt", "./port5.txt" };

    CStock stock;
	for (int i = 0; i <5; i++) {

        stock.Loadtradehistory(dataset[i]);		// load the daily return data--��������

		for (int j = 0; j <30;j++)
		{
		    updatetimes=0;
            delay=1.0;
			loop = j;
			cout << dataset[i] <<" : run = "<<loop<< endl;
			dim = stock.NumStock;
            hidesize=10;

            hidesize2=2;

            cout << dim << endl;

            lenOfChromosome = (2* hidesize + hidesize) + (hidesize2*hidesize+hidesize2)+ (1*hidesize2 + 1) + hidesize;
			cad = dim;
			max_sharp = -10;
            //nsga2(pop, gen, stock, loop, root);
            //dmoea(pop, gen, stock, loop, root);
            //dmoea_sharpnea(pop, gen, stock, loop, root);
			dmoea_sharpnea(pop, gen, stock, loop, root);
            //dmoea_sharp(pop, gen, stock, loop, root);

		}

	}

return 0;
}

