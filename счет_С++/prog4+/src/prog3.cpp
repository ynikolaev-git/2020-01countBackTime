/*
 * prog3.cpp
 *
 *  Created on: 5 мая 2020 г.
 *      Author: 1
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <iomanip>
#include <cmath>
#include "mx.hpp"
using namespace std;
const int maxLenLayer3 = 20000;

ofstream output;

// ========================================================================================
// Функция копирует в массив 1 count значений из массива 2
void outputLayer(double l1[maxLenLayer3][4], const int count, const int nl)
{
	double u, c;
	for (auto i=0; i<=count; i++)
		{
		u = (l1[i][2]+l1[i][3])/2;
		c = (l1[i][2]-l1[i][3])*(gamma-1)/4;
		output << nl << ";" << i << ";" << l1[i][0] <<  ";" << l1[i][1] << ";" << u << ";" <<
				c << ";" << u+c << ";" << u-c << ";" << endl;
		}
}

// ========================================================================================
// Функция копирует в массив 1 count значений из массива 2
void Layer2toLayer1(double l2[maxLenLayer3][4], const int count2, double l1[maxLenLayer3][4], int& count1)
{
	for (auto i=0; i<=count2; i++)
		{
		l1[i][0] = l2[i][0];
		l1[i][1] = l2[i][1];
		l1[i][2] = l2[i][2];
		l1[i][3] = l2[i][3];
		}
	count1 = count2;
}

// ========================================================================================
// Функция поиска минимального шага по времени (формулы в тетради)
double initTL(double layer[maxLenLayer3][4],   // последний построенный слой
	int tempmaxNumPointOnLayer1,        // число точек на слое
	double rw)                          // координата неподвижного поршня-стенки
{
	double tempu = (layer[1][2]+layer[1][3])/2;
	double tempc = (layer[1][2]-layer[1][3])*(gamma-1)/4;
	double mint = layer[1][0]+(rw-layer[1][1])/(tempu+tempc);  // время нулевой точки
return mint;
}

// ========================================================================================
int findmaxNumPointOnLayer2(double l1[maxLenLayer3][4],   // последний построенный слой
	double piston[maxLenLayer3][4],                       // траектория поршня
	int num1,                                      // число точек на последнем построенном слое
	double tL,                                     // время нового слоя
	int& cnop)                                     // номер отрезка достигнутый на траектории поршня
{
	int i = 0;
	while ((piston[i][0]>=tL)and(piston[i+1][0]>=tL))
			i++;
	// нашли отрезок поршня, соответствующий времени слоя
	cnop = i;

	double lambda = (tL-piston[i+1][0])/(piston[i][0]-piston[i+1][0]);
	double x = piston[i+1][1]+lambda*(piston[i][1]-piston[i+1][1]);
	// нашли пространственную координату поршня

	int N = num1-1;
	while ((x>l1[N+1][1])and(x>=l1[N][1]))
			N--;
	// нашли отрезок на последнем построенном слое, соответствующий координате х
return N+1;
}

// ======================
int main_prog3(){
	double p[maxLenLayer3][4];
	readPiston("piston10000_2.csv", p); // прочитали траекторию поршня

	double Layer1[maxLenLayer3][4], Layer2[maxLenLayer3][4];   // t,x,R,L
	int nr = 15000;
	double dr = (rw-rs)/nr;
	double drf = dr/100;
	int i = 0;

	double tOfCurrentLayer = ts;
	output.open("net3.csv");
	output << fixed << setprecision(10);
	output << "numOfLayer;numPointOnLayer;t;x;u;c;u+c;u-c;" << endl;

// строим первый слой
	int numOfLayer = 0;
	int maxNumPointOnLayer1;
	i = 0;
	double rCurrentPoint = rw;
	do  // заполняем слой сетки находящийся под слоем сжатого газа
		{
		Layer1[i][0] = tOfCurrentLayer;
		Layer1[i][1] = rCurrentPoint;
		Layer1[i][2] = 2*c0 / (gamma-1);
		Layer1[i][3] = - 2*c0 / (gamma-1);
		i ++;
		rCurrentPoint -= drf;
		}
	while (rCurrentPoint >= rf);
	while (rCurrentPoint-dr > rs)  // заполняем слой сетки находящийся под слоем траекторией поршня
		{
		Layer1[i][0] = tOfCurrentLayer;
		Layer1[i][1] = rCurrentPoint;
		Layer1[i][2] = 2*c0 / (gamma-1);
		Layer1[i][3] = - 2*c0 / (gamma-1);
		i ++;
		rCurrentPoint -= dr;
		}
	Layer1[i][0] = tOfCurrentLayer;
	Layer1[i][1] = rs;
	Layer1[i][2] = 2*c0 / (gamma-1);
	Layer1[i][3] = - 2*c0 / (gamma-1);

//	Layer2toLayer1(Layer2, maxNumPointOnLayer2, Layer1, maxNumPointOnLayer1);
	maxNumPointOnLayer1 = i;
//	outputLayer(Layer1, maxNumPointOnLayer1, numOfLayer);
	cout << "Построили слой номер = " << numOfLayer << "; t =" << tOfCurrentLayer << endl;
	numOfLayer++;

	int maxNumPointOnLayer2;
	int currentNumOnPiston;
	bool flag = 1;
// === начало основного цикла построения сетки
	while (Layer1[0][0]<tf)
		{
		if (flag)
			{
			tOfCurrentLayer = initTL(Layer1, maxNumPointOnLayer1, rw);
			tOfCurrentLayer -= (tOfCurrentLayer-Layer1[0][0])/10;
			newPoint2_prog3(Layer1, Layer2, tOfCurrentLayer);		// посчитали нулевую точку

			maxNumPointOnLayer2 = findmaxNumPointOnLayer2(Layer1, p, maxNumPointOnLayer1, tOfCurrentLayer, currentNumOnPiston);
			// посчитали maxNumPointOnLayer2
			for (i=1; i <= maxNumPointOnLayer2-1; i++)
				{
				flag *= newPoint1_prog3(Layer1, Layer2, tOfCurrentLayer, i);
				}
			newPoint3_prog3(Layer1, Layer2, p, tOfCurrentLayer, maxNumPointOnLayer1);
			// посчитали точку пересечения с поршнем
			}
		else
			{
			flag = 1;
			tOfCurrentLayer -= (tOfCurrentLayer-Layer1[0][0])/10;
			newPoint2_prog3(Layer1, Layer2, tOfCurrentLayer);		// посчитали нулевую точку

			maxNumPointOnLayer2 = findmaxNumPointOnLayer2(Layer1, p, maxNumPointOnLayer1, tOfCurrentLayer, currentNumOnPiston);
			// посчитали maxNumPointOnLayer2
			for (i=1; i <= maxNumPointOnLayer2-1; i++)
				{
				flag *= newPoint1_prog3(Layer1, Layer2, tOfCurrentLayer, i);
				}
			newPoint3_prog3(Layer1, Layer2, p, tOfCurrentLayer, maxNumPointOnLayer1);
			// посчитали точку пересечения с поршнем
			}
		if (flag)
			{
			if (Layer2[0][0]>=tf) outputLayer(Layer1, maxNumPointOnLayer1, numOfLayer);
			Layer2toLayer1(Layer2, maxNumPointOnLayer2, Layer1, maxNumPointOnLayer1);
//			outputLayer(Layer1, maxNumPointOnLayer2, numOfLayer);
			cout << "Построили слой номер = " << numOfLayer << "; t = " << tOfCurrentLayer <<
					"; точек на слое = " << maxNumPointOnLayer2 <<
					"; номер точки на поршне = "<< currentNumOnPiston << endl;
			numOfLayer++;
			}
		}
// === конец основного цикла
	outputLayer(Layer1, maxNumPointOnLayer2, numOfLayer);
	output.close();

	cout << "Расчет завершен " << endl;
	cout << "ts = " << ts << endl;
	cout << "rs = " << rs << endl;
	cout << "tf = " << tf << endl;
	cout << "rf = " << rf << endl;
	cout << "rw = " << rw << endl;
	cout << "c_ = " << pow(rho_,(gamma-1)/2) << endl;

	return 0;
}




