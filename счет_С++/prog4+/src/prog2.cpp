/*
 * prog2.cpp
 * Программа счета, когда в качестве расчетного слоя принимается траектория С^- характеристика
 *  Created on: 3 мая 2020 г.
 *      Author: 1
 */
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include <vector>
#include "mx.hpp"
using namespace std;

const int maxLen2 = 200000;
double p[maxLen2][4], Layer1[maxLen2][4], Layer2[maxLen2][4];   // t,x,R,L

struct prn1 {
	int qp;    // количество точек, которые сделали вклад в значения параметров газа
	double t;
	double x;
	double u;
	double c;
};
vector<prn1> netFinal;

void crossLayerTf(double l1[maxLen2][4], int count)
{
	double lambda, u0, u1, c0, c1;
	prn1 p;
	for (int i=0; i<count;i++)
		if ((l1[i][0]<=tf)and(l1[i+1][0]>=tf))
			{
			lambda = (tf-l1[i][0])/(l1[i+1][0]-l1[i][0]);
			p.t = l1[i][0] + lambda*(l1[i+1][0]-l1[i][0]);
			p.x = l1[i][1] + lambda*(l1[i+1][1]-l1[i][1]);
			u0 = (l1[i][2]+l1[i][3])/2;
			c0 = (l1[i][2]-l1[i][3])*(gamma-1)/4;
			u1 = (l1[i+1][2]+l1[i+1][3])/2;
			c1 = (l1[i+1][2]-l1[i+1][3])*(gamma-1)/4;
			p.u = u0 + lambda*(u1-u0);
			p.c = c0 + lambda*(c1-c0);
			if (p.u < 0.1) netFinal.push_back(p);  // добавил фильтр отбрасывающий точки, примыкающие к (t_f, r_f)
			break;
			}
}

void L2ToL1_OutToFile(double l1[maxLen2][4], double l2[maxLen2][4], const int count)
{
	for (int i=0; i <= count; i++)
		{
		l1[i][0] = l2[i][0];
		l1[i][1] = l2[i][1];
		l1[i][2] = l2[i][2];
		l1[i][3] = l2[i][3];
//		double u = (l1[i][2]+l1[i][3])/2;
//		double c = (l1[i][2]-l1[i][3])*(gamma-1)/4;
//		output << numOfLayer << ";" << i << ";" << l1[i][0] <<  ";" << Layer1[i][1] << ";" << u << ";" <<
//				c << ";" << u+c << ";" << u-c << ";" << endl;
		}
}

//  далее основная процедура ======================================================================
int main_prog2(){
//	double p[maxLen2][4], Layer1[maxLen2][4], Layer2[maxLen2][4];   // t,x,R,L
	cout << fixed << setprecision(16);
	readPiston("piston5.csv", p); // прочитали траекторию поршня

	double timeEndRegular, rEndRegular;
	int nr = 25000;
	double dr = (rw-rs)/nr;
	int i;

	// инициализируем первый слой
	int numPointOnLayer = 0;
	int numOfLayer = 0;
	int maxNumPointOnLayer1;
	Layer1[numPointOnLayer][0] = ts;
	Layer1[numPointOnLayer][1] = rs;
	Layer1[numPointOnLayer][2] = 2*c0/(gamma-1);
	Layer1[numPointOnLayer][3] = -2*c0/(gamma-1);
	cout << "построили слой номер = " << numOfLayer << endl;
// строим следующие слои до линии r=rw
	cout << "Строим область 1: номер слоя 1, номер на поршне = " << numOnPiston << endl;
	while (Layer1[0][1] < rw)
		{
		numOfLayer ++;
		Layer2[0][0] = ts;
		Layer2[0][1] = rs + numOfLayer*dr;
		Layer2[0][2] = 2*c0/(gamma-1);
		Layer2[0][3] = -2*c0/(gamma-1);
		for (i=1; i < numOfLayer*2; i++)
			{
			newPoint1(Layer2[i-1], Layer1[i-1], Layer2[i]);
			}
		while (newPoint3_prog2(p[numOnPiston], p[numOnPiston-1], Layer2[numOfLayer*2-1], Layer2[numOfLayer*2]) != 1)
			{
			numOnPiston -=1;
			}

		L2ToL1_OutToFile(Layer1, Layer2, numOfLayer*2);
		}
		cout << "построили область 1: номер слоя = " << numOfLayer << " номер на поршне = " << numOnPiston << endl;

	maxNumPointOnLayer1 = numOfLayer*2;
// строим следующие слои вдоль r=rw пока минус характеристики еще пересекают траекторию поршня
	cout << "Строим область 2: номер слоя = " << numOfLayer << " номер на поршне = " << numOnPiston << endl;
	while (numOnPiston > 1)
		{
		newPoint2(Layer1[1], Layer2[0]);
		for (i=2; i <= maxNumPointOnLayer1; i++)
			{
			newPoint1(Layer2[i-2], Layer1[i], Layer2[i-1]);
			}
		while (newPoint3_prog2(p[numOnPiston], p[numOnPiston-1], Layer2[maxNumPointOnLayer1-1], Layer2[maxNumPointOnLayer1]) != 1)
			{
			numOnPiston -=1;
			if (numOnPiston == 0)
				{
				timeEndRegular = Layer2[maxNumPointOnLayer1-1][0];
				rEndRegular = Layer2[maxNumPointOnLayer1-1][1];
				break;
				}
			}
		L2ToL1_OutToFile(Layer1, Layer2, maxNumPointOnLayer1);
		numOfLayer ++;
//		cout << "область 2: построили слой номер = " << numOfLayer <<
//			", номер на поршне = " << numOnPiston << endl;
		}
	cout << "Построили область 2: номер слоя = " << numOfLayer << " номер на поршне = " << numOnPiston << endl;
// строим следующие слои вдоль линии r=rw когда траектория поршня закончилась
	cout << "Строим область 3: номер слоя = " << numOfLayer << " число точек на слое = " << maxNumPointOnLayer1 << endl;
	while (maxNumPointOnLayer1 > 1)
		{
		newPoint2(Layer1[1], Layer2[0]);
		for (i=2; i <= maxNumPointOnLayer1; i++)
			{
			newPoint1(Layer2[i-2], Layer1[i], Layer2[i-1]);
			}
		maxNumPointOnLayer1 -= 1;
		if ((Layer2[0][0]<=tf)and(Layer2[maxNumPointOnLayer1][0]>=tf))
			{
//			cout << "слой номер = " << numOfLayer << ", пересекли линию t == tf" << endl;
			crossLayerTf(Layer2, maxNumPointOnLayer1);
			}
		else{
			// не пересекли t=tf, проверяем что продолжение С^- характеристики прилетит в точку (tf,rf)
			double k = (Layer2[maxNumPointOnLayer1][2]+Layer2[maxNumPointOnLayer1][3])/2;
			k = k - (gamma-1)*(Layer2[maxNumPointOnLayer1][2]-Layer2[maxNumPointOnLayer1][3])/4;
//			cout << "характеристика приходит в точку " << Layer2[maxNumPointOnLayer1][1] + k*(tf-Layer2[maxNumPointOnLayer1][0]) << endl;
			}

		L2ToL1_OutToFile(Layer1, Layer2, maxNumPointOnLayer1);
		numOfLayer ++;
//		cout << "область 3: построили слой номер = " << numOfLayer << ", число точек на слое = " << maxNumPointOnLayer1 << endl;
		}
	cout << "Построили область 3: номер слоя = " << numOfLayer << " число точек на слое = " << maxNumPointOnLayer1 << endl;
	cout << "Расчет завершен " << endl;
//	output.close();

// ================================================================================================
	ofstream outNetFinal;
	outNetFinal.open("netFinal.csv");
	outNetFinal << fixed << setprecision(16);
	outNetFinal << "t;x;u;c;rho" << endl;
	for (auto pnt:netFinal)
		{
		cout << "t = " << pnt.t << ", x = " << pnt.x << ", u = " << pnt.u << ", c = " << pnt.c << ", rho = " << pow(pnt.c, 2/(gamma-1)) << endl;
		outNetFinal << pnt.t << ";" << pnt.x << ";" << pnt.u << ";" << pnt.c << ";" << pow(pnt.c, 2/(gamma-1)) << endl;
		}
	outNetFinal.close();
// ================================================================================================
	ofstream tabl1;
	tabl1.open("tabl1.csv");
	tabl1 << fixed << setprecision(16);
	tabl1 << "r_s;t_f;r_f;maxNum;u_max;c_max;rho_max" << endl;
	tabl1 << rs << ";" << tf << ";" << rf << ";" << maxNum << ";" << u_max << ";" << c_max << ";" << rho_max << endl;
	tabl1.close();
// ================================================================================================
	ofstream tabl2;
	tabl2.open("tabl2.csv");
	tabl2 << fixed << setprecision(160);
	tabl2 << "gamma;nu;rho_;m;n;КоличествоCлоев" << endl;
	tabl2 << gamma << ";" << nu << ";" << rho_ << ";" << m << ";" << nr << ";" << numOfLayer << endl;
	tabl2.close();
// ================================================================================================
	ofstream tabl3;
	tabl3.open("tabl3.csv");
	tabl3 << fixed << setprecision(16);
	tabl3 << "dm;t_E1;r_E1;r_E2;T%;d_reg%" << endl;
	tabl3 << dm(rs, rho0, netFinal[0].x, pow((netFinal[0].c+netFinal[netFinal.size()-1].c)/2, 2/(gamma-1)), rw)
			<< ";" << timeEndRegular << ";" << rEndRegular << ";" << netFinal[0].x << ";" << 100*timeEndRegular/tf
			<< ";" << 100*(rw-netFinal[0].x)/(rw-rf) << endl;
	tabl3.close();
// ================================================================================================
	cout << "=====================================================================================" << endl;
	cout << "Погрешность масс менее: " << dm(rs, rho0, netFinal[0].x,
			pow((netFinal[0].c+netFinal[netFinal.size()-1].c)/2, 2/(gamma-1)), rw) << "%" <<endl;
	cout << "последний момент времени когда сетка была регулярной (точка Е1): " << timeEndRegular << endl;
	cout << "последнее значение переменной r когда сетка была регулярной (точка Е1): " << rEndRegular << endl;
	cout << "значение переменной r когда сетка снова стала регулярной (точка на линии t = tf, точка Е2): " << netFinal[0].x << endl;
	cout << 100*timeEndRegular/tf << " % времени сжатия сетка строилась регулярно (T_reg)" << endl;
//	cout << "расстояние по области где сетка не строится: " << netFinal[0].x-rEndRegular << endl;
//	cout << 100*(netFinal[0].x-rEndRegular)/(rw-rf) << " %, соотношение расстояния отсутствия регулярности и расстояния между финишом поршня и неподвижной стенкой" << endl;
//	cout << 100*(netFinal[0].x-rEndRegular)/(rw-rs) << " %, соотношение расстояния отсутствия регулярности и расстояния между стартом поршня и неподвижной стенкой" << endl;
	cout << 100*(rw-netFinal[0].x)/(rw-rf) << " % расстояния в момент сжатия где построена регулярная сетка, d_reg" << endl;

	cout << "ts = " << ts << endl;
	cout << "rs = " << rs << endl;
	cout << "tf = " << tf << endl;
	cout << "rf = " << rf << endl;
	cout << "rw = " << rw << endl;
	cout << "c_ = " << pow(rho_,(gamma-1)/2) << endl;

	cout << "==  Расчет окончен  =================================================================" << endl;
	return 0;
}


