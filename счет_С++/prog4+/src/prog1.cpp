//============================================================================
// Name        : prog1.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : start 17-04-2020
// Description : end   27-04-2020
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
using namespace std;

#include "mx.hpp"

//extern double rho0, c0, gamma, rho_, m, ts;
//extern double c0, gamma, rho_, m, ts;
//extern double rs, tf, rf, pi;
//extern int nu;
//extern const int len_piston;
double nr = 100;  // параметр шага по пространственной переменной, определяет точность расчетов
double nt = 100;
int flag = 1;
int doit = 1;
const int n = 100;
int numOnLayer;
int numLayer;

double layer1[n][4];  // t,r,R,L == t,x,R,L
double layer2[n][4];

double calc_rw_old() {
	double res;
	if (nu == 0) {
        res = rs + m/rho0;
	}
	if (nu == 1) {
        res = pow(rs*rs + m/(pi*rho0), 1.0/2.0);
	}
	if (nu == 2) {
        res = pow(rs*rs*rs + (3*m)/(4*pi*rho0), 1.0/3.0);
	}
	return res;
}

//rw = calc_rw();
double dr = (rw-rs)/nr;
double dt = 0.0001;
double layerEnd_time;
double u, c;
int i = -1;



//=======================================================================
//void newPoint2(double* p1, double* p2) {
// Пересечение плюс характеристики с прямой r=rw
// Из р1 выпускаем С+, в р2 новая точка
// считаем первое приближение
//    double a = Cp(p1[2],p1[3]);
//    p2[1] = rw;
//    p2[0] = p1[0]+(p2[1]-p1[1])/a;
//    p2[2] = p1[2]+(p2[0]-p1[0])*Fp(p1[1],p1[2],p1[3]);
//    p2[3] = -p2[2];
// уточняем полученные значения
//    a = (Cp(p1[2],p1[3])+Cp(p2[2],p2[3]))/2;
//    p2[0] = p1[0]+(rw-p1[1])/a;
//    p2[2] = p1[2]+(p2[0]-p1[0])*(Fp(p1[1],p1[2],p1[3])+Fp(p2[1],p2[2],p2[3]))/2;
//    p2[3] = -p2[2];
//}


//=======================================================================
bool newPoint3(double* p1, double* p2, double* p3) {
// в р1, р2 отрезок траектории поршня
// в р3 точка сетки из которой ищем пересечение с траекторией поршня. Если пересечение нашли (пересечение может быть на другом отрезке траеткории),
// то p3 сохраняем новую точку и функция возврщает true, иначе p3 не меняем и функция возвращает false.
    bool flag3 = 0;
    double u=(p3[2]+p3[3])/2 - (gamma-1)*(p3[2]-p3[3])/4;
   double a = p2[1]-p1[1]+u*(p1[0]-p2[0]);
   double l, t;
    if (a == 0) {
        l = 10;
    }
    else {
    	l = (p3[1]+u*p1[0]-u*p3[0]-p1[1])/a;
    	t = p1[0] + l*(p2[0]-p1[0]);
    }
    if ((0<=l)and(l<=1)and(p3[0]<=t)) {
        p3[0] = p1[0] + l*(p2[0]-p1[0]);
        p3[1] = p1[1] + l*(p2[1]-p1[1]);
        p3[2] = p1[2] + l*(p2[2]-p1[2]);
        p3[3] = p1[3] + l*(p2[3]-p1[3]);
        flag3 = 1;
    }
    return flag3;
}

//=======================================================================
bool check(int m) {
//	return 1;
	if (m == 0) return 1;
	if (layer2[m][1] >= layer2[m-1][1]){
		cout << "обрезали слой, количество точек стало = " << m << endl;
		return 0;  // регулярность сетки нарушена
	}
	else return 1;  // регулярность не нарушена
}


//=======================================================================
int main1() {
	double p[10000][4];
	ofstream output;
	output.open("net.csv");
	output << "numOfLayer;numOnLayer;t;x;u;c;u+c;u-c;" << endl;

// читаем траекторию поршня
	readPiston("piston32.csv", p);
// прочитали траекторию поршня
// ==================================================================
//  Далее считаем слой номер 1
	i = 0;
	flag = 1;
	while (rw-i*dr > rs) {
		layer1[i][0] = ts;
		layer1[i][1] = rw-i*dr;
		layer1[i][2] = 2*c0 / (gamma-1);
		layer1[i][3] = - 2*c0 / (gamma-1);
		u = (layer1[i][2]+layer1[i][3])/2;
		c = (layer1[i][2]-layer1[i][3])*(gamma-1)/4;
		output << "1" << ";" << i << ";" << layer1[i][0] <<  ";" << layer1[i][1] << ";" << u << ";" <<
				c << ";" << u+c << ";" << u-c << ";" << endl;
		i ++;
	}
	layer1[i][0] = ts;
	layer1[i][1] = rs;
	layer1[i][2] = 2*c0 / (gamma-1);
	layer1[i][3] = - 2*c0 / (gamma-1);
	u = (layer1[i][2]+layer1[i][3])/2;
	c = (layer1[i][2]-layer1[i][3])*(gamma-1)/4;
	output << "1" << ";" << i << ";" << layer1[i][0] <<  ";" << layer1[i][1] << ";" << u << ";" <<
			c << ";" << u+c << ";" << u-c << ";" << endl;
	layerEnd_time = layer1[i][0];
	numLayer = 1;
//	R_pred = layer1[0][2];
	cout << "построили слой номер 1" << endl;
	cout << "на слое точек (счет начинается с 0) = " << i+1 << endl;
// ==================================================================
// Далее строим сетку
	numOnLayer = i;  // +1 это число точек на нечетном слое, нумерация с 0
	flag = 2;
	bool withOutPiston = 0;
//	flag == 1, строим слой с правой точкой на стенке и левой на траектории сжимающего поршня (нечентные номера)
//	flag == 2, строим слой без точек на поршнях (четные номера)
//	flag == 3, строим слой с правой точкой на стенке, траектория сжимающего поршня закончилась
	doit = 1;
	while (doit) {
		if (flag == 1) {
			newPoint2(layer1[0], layer2[0]);
			for (i=0; i <= numOnLayer-2; i++) {
				newPoint1(layer1[i], layer1[i+1], layer2[i+1]);
				if (check(i+1) != 1){
					numOnLayer = i;
					withOutPiston = 1;
				}
			}
			if (withOutPiston == 0) {
				while (newPoint3(p[numOnPiston], p[numOnPiston-1], layer1[i]) != 1) {
					numOnPiston -=1;
					}
				layer2[i+1][0] = layer1[i][0];
				layer2[i+1][1] = layer1[i][1];
				layer2[i+1][2] = layer1[i][2];
				layer2[i+1][3] = layer1[i][3];
				if (check(i+1) != 1){  // check
					numOnLayer = i;
					withOutPiston = 1;
				}
				}
			for (i=0; i <= numOnLayer; i++) {
				layer1[i][0] = layer2[i][0];
				layer1[i][1] = layer2[i][1];
				layer1[i][2] = layer2[i][2];
				layer1[i][3] = layer2[i][3];
				u = (layer1[i][2]+layer1[i][3])/2;
				c = (layer1[i][2]-layer1[i][3])*(gamma-1)/4;
				output << numLayer+1 << ";" << i << ";" << layer1[i][0] <<  ";" << layer1[i][1] << ";" << u << ";" <<
						c << ";" << u+c << ";" << u-c << ";" << endl;
				}
			flag = 2;
			numLayer ++;
			cout << "построили слой номер = " << numLayer << endl;
			layerEnd_time = layer1[i][0];
		}
		else if (flag == 2) {
			for (i=0; i <= numOnLayer-1; i++) {
				newPoint1(layer1[i], layer1[i+1], layer2[i]);
				if (check(i) != 1){  // check
					numOnLayer = i;
					withOutPiston = 1;
//					flag = 3;
				}
				}
//			corrLayer2();
			for (i=0; i <= numOnLayer-1; i++) {
				layer1[i][0] = layer2[i][0];
				layer1[i][1] = layer2[i][1];
				layer1[i][2] = layer2[i][2];
				layer1[i][3] = layer2[i][3];
				u = (layer1[i][2]+layer1[i][3])/2;
				c = (layer1[i][2]-layer1[i][3])*(gamma-1)/4;
				output << numLayer+1  << ";" << i << ";" << layer1[i][0] <<  ";" << layer1[i][1] << ";" << u << ";" <<
						c << ";" << u+c << ";" << u-c << ";" << endl;
				}
			flag = 1;
			numLayer ++;
			cout << "построили слой номер = " << numLayer << endl;
			layerEnd_time = layer1[i][0];
		}
		for (i=0; i <= numOnLayer; i++) {
			if (layer1[i][0] < tf) {
				doit = 1;
//					doit = 0;
				i = numOnLayer;
			} else {
				doit = 0;
			}
		}
		if (numOnLayer <=1) doit = 0;

	}
	cout << "Расчет успешно завершен" << endl;
	cout << "numLayer = " << numLayer << endl;
	output.close();
	cout << "с* = " << pow(rho_, (gamma-1)/2.0) << endl;
	return 0;
}

