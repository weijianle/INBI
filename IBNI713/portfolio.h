
#ifndef __CSTOCK_H_
#define __CSTOCK_H_

#include "global.h"

// class for tsp test instances
class CStock {
public:
	int     NumStock;			//�ʲ�����
    double   mean[250];	//��ֵ
    double   var[250];	//����

	double covar[250][250];//����
	void   Loadtradehistory(const char* pch);

	CStock();
	~CStock();
};

CStock::CStock() {

}

CStock::~CStock() {

}

void CStock::Loadtradehistory(const char* pch) {
	int i=0, j=0;
	double m=0.0, v=0.0, cv=0.0, index1=0.0, index2=0.0;

	std::ifstream readf(pch);
	readf >> NumStock;
    double max_mean=0,max_var=0;
	for (i = 0; i < NumStock; i++) {
		readf >> m;
		readf >> v;
		mean[i]=m;
		var[i]=v;

	}


	for (i = 0; i < NumStock; i++)
		for (j = i; j < NumStock; j++)
		{
			readf >> index1;
			readf >> index2;
			readf >> cv;
			covar[i][j] = cv * var[i] * var[j];
			covar[j][i] = covar[i][j];
		}
	readf.close();
}

#endif