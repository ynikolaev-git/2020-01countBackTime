/*
 * mx.cpp
 *
 *  Created on: 3 мая 2020 г.
 *      Author: 1
 */
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <cmath>
using namespace std;

//#include "mx.hpp"

double rho0 = 1;
double c0 = 1;
//double gamma = 1.4;
double gamma = 1.6666666666666667;
double rho_ = 10000;
double m = 10;
int nu = 1;
double ts, tf, rf, rs, rw, maxNum, u_max, c_max, rho_max;
double pi = 3.1415926535;

double pistonEnd;
int numOnPiston;
const int maxLenPiston = 200000;
const int maxLenLayer = 200000;

double calc_rw() {
	double res;
	if (nu == 0) {
//        res = rs + m/rho0;
        res = rf + m/rho_;
	}
	if (nu == 1) {
//        res = pow(rs*rs + m/(pi*rho0), 1.0/2.0);
        res = pow(rf*rf + m/(pi*rho_), 1.0/2.0);
	}
	if (nu == 2) {
//        res = pow(rs*rs*rs + (3*m)/(4*pi*rho0), 1.0/3.0);
        res = pow(rf*rf*rf + (3*m)/(4*pi*rho_), 1.0/3.0);
	}
	return res;
}

// double rw = calc_rw();

//=======================================================================
// наклон положительной характеристики
double Cp(double R,double L) {
    return R*(gamma+1)/4+L*(3-gamma)/4;
}

//=======================================================================
// наклон отрицательной характеристики
double Cm(double R,double L){
    return L*(gamma+1)/4+R*(3-gamma)/4;
}

//=======================================================================
double Fp(double x,
		double R,
		double L) {
    return (-nu)*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
double Fm(double x,
		double R,
		double L) {
    return nu*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
void newPoint1(double p1[4],
		double p2[4],
		double p3[4]) {
// Пересечение плюс и минус характеристик в р
// Из р1 выпускаем С-, из р2 С+
// считаем первое приближение
    double a = Cm(p1[2],p1[3]);
    double b = Cp(p2[2],p2[3]);
    p3[0] = (a*p1[0]-b*p2[0]+p2[1]-p1[1])/(a-b);
    p3[1] = p1[1]+(p2[1]-p1[1]+b*(p1[0]-p2[0]))*a/(a-b);
    p3[2] = p2[2]+(p3[0]-p2[0])*Fp(p2[1],p2[2],p2[3]);
    p3[3] = p1[3]+(p3[0]-p1[0])*Fm(p1[1],p1[2],p1[3]);
// уточняем полученные значения
    a = (Cm(p1[2],p1[3])+Cm(p3[2],p3[3]))/2;
    b = (Cp(p2[2],p2[3])+Cp(p3[2],p3[3]))/2;
    p3[0] = (a*p1[0]-b*p2[0]+p2[1]-p1[1])/(a-b);
    p3[1] = p1[1]+(p2[1]-p1[1]+b*(p1[0]-p2[0]))*a/(a-b);
    p3[2] = p2[2]+(p3[0]-p2[0])*(Fp(p2[1],p2[2],p2[3])+Fp(p3[1],p3[2],p3[3]))/2;
    p3[3] = p1[3]+(p3[0]-p1[0])*(Fm(p1[1],p1[2],p1[3])+Fm(p3[1],p3[2],p3[3]))/2;
}

//=======================================================================
bool newPoint1_prog3(double l1[maxLenLayer][4],  // массив слоя 1, данные используем для расчета
		double l2[maxLenLayer][4],               // массив слоя 2, здесь рассчитываем точку
		double tLayer,                           // время на слое
		int numPoint) {                          // номер рассчитываемой точки
    double dt = tLayer - l1[0][0];
    double lambda, u0, c0, um1, cm1, up1, cp1;
    double pp[4], pm[4], pnp1[4];
    bool flg = 1;
    um1 = (l1[numPoint-1][2]+l1[numPoint-1][3])/2;
    cm1 = (l1[numPoint-1][2]-l1[numPoint-1][3])*(gamma-1)/4;
    up1 = (l1[numPoint+1][2]+l1[numPoint+1][3])/2;
    cp1 = (l1[numPoint+1][2]-l1[numPoint+1][3])*(gamma-1)/4;
    u0 = (l1[numPoint][2]+l1[numPoint][3])/2;
    c0 = (l1[numPoint][2]-l1[numPoint][3])*(gamma-1)/4;

    // считаем lambda+ и p+
    lambda = l1[numPoint][1]-l1[numPoint+1][1] - dt*(up1+cp1);
    lambda = lambda/(l1[numPoint][1]-l1[numPoint+1][1]+dt*(u0+c0-up1-cp1));
    flg = ((lambda>=0)and(lambda<=1));
    pp[0] = l1[0][0];
    pp[1] = l1[numPoint+1][1]+lambda*(l1[numPoint][1]-l1[numPoint+1][1]);
    pp[2] = l1[numPoint+1][2]+lambda*(l1[numPoint][2]-l1[numPoint+1][2]);
    pp[3] = l1[numPoint+1][3]+lambda*(l1[numPoint][3]-l1[numPoint+1][3]);

	// считаем lambda- и pm
    lambda = -dt*(u0-c0);
    lambda = lambda/(l1[numPoint-1][1]-l1[numPoint][1]+dt*(um1-cm1-u0+c0));
    if ((lambda>=0)and(lambda<=1))
		{
		pm[0] = l1[0][0];
		pm[1] = l1[numPoint][1]+lambda*(l1[numPoint-1][1]-l1[numPoint][1]);
		pm[2] = l1[numPoint][2]+lambda*(l1[numPoint-1][2]-l1[numPoint][2]);
		pm[3] = l1[numPoint][3]+lambda*(l1[numPoint-1][3]-l1[numPoint][3]);
		}
	else
		{
	    lambda = l1[numPoint][1]-l1[numPoint+1][1]-dt*(up1-cp1);
	    lambda = lambda/(l1[numPoint][1]-l1[numPoint+1][1]+dt*(u0-c0-up1+cp1));
		pm[0] = l1[0][0];
		pm[1] = l1[numPoint+1][1]+lambda*(l1[numPoint][1]-l1[numPoint+1][1]);
		pm[2] = l1[numPoint+1][2]+lambda*(l1[numPoint][2]-l1[numPoint+1][2]);
		pm[3] = l1[numPoint+1][3]+lambda*(l1[numPoint][3]-l1[numPoint+1][3]);
		}
    flg = ((flg)and((lambda>=0)and(lambda<=1)));

	// считаем первое приближение
    double a = Cm(pm[2],pm[3]);
    double b = Cp(pp[2],pp[3]);
    pnp1[0] = (a*pm[0]-b*pp[0]+pp[1]-pm[1])/(a-b);
    pnp1[1] = pm[1]+(pp[1]-pm[1]+b*(pm[0]-pp[0]))*a/(a-b);
    pnp1[2] = pp[2]+(pnp1[0]-pp[0])*Fp(pp[1],pp[2],pp[3]);
    pnp1[3] = pm[3]+(pnp1[0]-pm[0])*Fm(pm[1],pm[2],pm[3]);
    double unp1 = (pnp1[2]+pnp1[3])/2;
    double cnp1 = (pnp1[2]-pnp1[3])*(gamma-1)/4;

    l2[numPoint][0] = pnp1[0];
    l2[numPoint][1] = pnp1[1];
    l2[numPoint][2] = pnp1[2];
    l2[numPoint][3] = pnp1[3];

// =====  делаем пересчет =========================
    // уточняем лямбда+ и pp
    lambda = 2*(l1[numPoint][1]-l1[numPoint+1][1])-dt*(up1+cp1+unp1+cnp1);
    lambda = lambda/(2*(l1[numPoint][1]-l1[numPoint+1][1]) + dt*(u0-c0-up1+cp1));
    flg = ((lambda>=0)and(lambda<=1));
    pp[0] = l1[0][0];
    pp[1] = l1[numPoint+1][1]+lambda*(l1[numPoint][1]-l1[numPoint+1][1]);
    pp[2] = l1[numPoint+1][2]+lambda*(l1[numPoint][2]-l1[numPoint+1][2]);
    pp[3] = l1[numPoint+1][3]+lambda*(l1[numPoint][3]-l1[numPoint+1][3]);
	// утоняем lambda- и pm
    lambda = -dt*(u0-c0+unp1-cnp1);
    lambda = lambda/(2*(l1[numPoint-1][1]-l1[numPoint][1])+dt*(um1-cm1-u0+c0));
    if ((lambda>=0)and(lambda<=1))
		{
		pm[0] = l1[0][0];
		pm[1] = l1[numPoint][1]+lambda*(l1[numPoint-1][1]-l1[numPoint][1]);
		pm[2] = l1[numPoint][2]+lambda*(l1[numPoint-1][2]-l1[numPoint][2]);
		pm[3] = l1[numPoint][3]+lambda*(l1[numPoint-1][3]-l1[numPoint][3]);
		}
	else
		{
		lambda = 2*(l1[numPoint][1]-l1[numPoint+1][1])-dt*(up1-cp1+unp1-cnp1);
	    lambda = lambda/(2*(l1[numPoint][1]-l1[numPoint+1][1])+dt*(u0-c0-up1+cp1));
		pm[0] = l1[0][0];
		pm[1] = l1[numPoint+1][1]+lambda*(l1[numPoint][1]-l1[numPoint+1][1]);
		pm[2] = l1[numPoint+1][2]+lambda*(l1[numPoint][2]-l1[numPoint+1][2]);
		pm[3] = l1[numPoint+1][3]+lambda*(l1[numPoint][3]-l1[numPoint+1][3]);
		}
    flg = ((flg)and((lambda>=0)and(lambda<=1)));
    l2[numPoint][2] = pp[2]+(pnp1[0]-pp[0])*(Fp(pp[1],pp[2],pp[3])+Fp(pnp1[1],pnp1[2],pnp1[3]))/2;
    l2[numPoint][3] = pm[3]+(pnp1[0]-pm[0])*(Fm(pm[1],pm[2],pm[3])+Fm(pnp1[1],pnp1[2],pnp1[3]))/2;
return flg;
}


//=======================================================================
void newPoint2(double* p1, double* p2) {
// Пересечение плюс характеристики с прямой r=rw
// Из р1 выпускаем С+, в р2 новая точка
// считаем первое приближение
    double a = Cp(p1[2],p1[3]);
    p2[1] = rw;
    p2[0] = p1[0]+(p2[1]-p1[1])/a;
    p2[2] = p1[2]+(p2[0]-p1[0])*Fp(p1[1],p1[2],p1[3]);
    p2[3] = -p2[2];
// уточняем полученные значения
    a = (Cp(p1[2],p1[3])+Cp(p2[2],p2[3]))/2;
    p2[0] = p1[0]+(rw-p1[1])/a;
    p2[2] = p1[2]+(p2[0]-p1[0])*(Fp(p1[1],p1[2],p1[3])+Fp(p2[1],p2[2],p2[3]))/2;
    p2[3] = -p2[2];
}

//=======================================================================
void newPoint2_prog3(double l1[maxLenLayer][4], double l2[maxLenLayer][4], double tLayer) {
// Пересечение плюс характеристики с прямой r=rw
    double pp[4];
// считаем
    double dt = tLayer - l1[0][0];
    double lambda, u, c, up1, cp1;
    u = (l1[0][2]+l1[0][3])/2;
    c = (l1[0][2]-l1[0][3])*(gamma-1)/4;
    up1 = (l1[1][2]+l1[1][3])/2;
    cp1 = (l1[1][2]-l1[1][3])*(gamma-1)/4;

    // считаем lambda+ и p+
    lambda = rw-l1[1][1]-dt*(up1+cp1);
    lambda = lambda/(l1[0][1]-l1[1][1]+dt*(u+c-up1-cp1));
    pp[0] = l1[0][0];
    pp[1] = l1[1][1]+lambda*(l1[0][1]-l1[1][1]);
    pp[2] = l1[1][2]+lambda*(l1[0][2]-l1[1][2]);
    pp[3] = l1[1][3]+lambda*(l1[0][3]-l1[1][3]);

	// считаем первое приближение
	double a = Cp(pp[2],pp[3]);
    l2[0][1] = rw;
    l2[0][0] = pp[0]+(l2[0][1]-pp[1])/a;
    l2[0][2] = pp[2]+(l2[0][0]-pp[0])*Fp(pp[1],pp[2],pp[3]);
    l2[0][3] = -l2[0][2];
}

//=======================================================================
void readPiston(
		const char* fName,
		double p[maxLenPiston][4])
{
	ifstream input(fName);
	if (input) {
		cout << fName << " is opened" << endl;
	} else {
		cout << fName << " is not opened" << endl;
	}
	string rd;
	getline(input, rd);  // прочитали строку с названиями колонок
	int i = -1;
	while (getline(input, rd, ';')) {
		i ++;
		getline(input, rd, ';');
		p[i][0] = atof(rd.c_str());
		getline(input, rd, ';');
		p[i][1] = atof(rd.c_str());
		getline(input, rd, ';');
		p[i][2] = atof(rd.c_str());
		getline(input, rd);
		p[i][3] = atof(rd.c_str());
	}
	numOnPiston = i;
	pistonEnd = p[0][0];
	ts = p[i][0];
	rs = p[i][1];
	tf = p[0][0];
	rf = p[0][1];
	u_max = (p[0][2]+p[0][3])/2;
	c_max = (p[0][2]-p[0][3])*(gamma-1)/4;
	rho_max = pow(c_max, 2/(gamma-1));
	maxNum = i;
	rw = calc_rw();
	cout << "прочитали поршень из файла" << endl;

}

//=======================================================================
bool newPoint3_prog2(double p1[4], double p2[4], double p3[4], double p4[4]) {
// в р1, р2 отрезок траектории поршня
// в р3 точка сетки из которой ищем пересечение с траекторией поршня. Если пересечение нашли (пересечение может быть на другом отрезке траеткории),
// то в p4 сохраняем новую точку и функция возврщает true, иначе функция возвращает false.
    bool flag3 = 0;
    double u=(p3[2]+p3[3])/2 - (gamma-1)*(p3[2]-p3[3])/4;
    double a = p2[1]-p1[1]+u*(p1[0]-p2[0]);
    double l;
    double unew;
    if (a == 0) {
        l = 10;
    }
    else {
    	l = (p3[1]+u*p1[0]-u*p3[0]-p1[1])/a;
    }
    if ((0<=l)and(l<=1)) {
//    if (((0<=l)and(l<=1))or(numOnPiston == 100)) {
        p4[0] = p1[0] + l*(p2[0]-p1[0]);
        p4[1] = p1[1] + l*(p2[1]-p1[1]);
//        p4[2] = p1[2] + l*(p2[2]-p1[2]);
//        p4[3] = p1[3] + l*(p2[3]-p1[3]);
        unew = (p1[2]+p1[3])/2 + l*((p2[2]+p2[3])/2-(p1[2]+p1[3])/2);
        p4[3] = p3[3] + (p4[0]-p3[0])*Fm(p3[1],p3[2],p3[3]);
        p4[2] = 2*unew-p4[3];
        flag3 = 1;
    }
    return flag3;
}

//=======================================================================
void newPoint3_prog3(double l1[maxLenLayer][4],  // последний построенный слой
		double l2[maxLenLayer][4],              // слой, который строим
		double piston[maxLenPiston][4],         // траектория поршня
		double tL,                              // значение времени на строящемся поршне
		int count) {                            // число точек на последнем построенном слое
	int i = 0;
	while ((piston[i][0]>=tL)and(piston[i+1][0]>=tL))
			i++;
	// нашли отрезок поршня, соответствующий времени слоя

	double lambda = (tL-piston[i+1][0])/(piston[i][0]-piston[i+1][0]);
	double x = piston[i+1][1]+lambda*(piston[i][1]-piston[i+1][1]);
	double u = (piston[i+1][2]+piston[i+1][3])/2;
	u = u + lambda*((piston[i][2]+piston[i][3])/2 - (piston[i+1][2]+piston[i+1][3])/2);
	// нашли пространственную координату поршня и значение скорости газа в ней (совпадает со скоростью поршня)

	int N = count-1;
	while ((x>l1[N+1][1])and(x>l1[N][1]))
			N--;
	// нашли отрезок на последнем построенном слое, соответствующий координате х

	// считаем lambda- и pm
    double dt = tL - l1[0][0];
	double rr, ur, cr, rl, ul, cl;
    double pm[4], pnp1[4];
	for (int k=count; k>N-1; k--)
		{
		rl = l1[k][1];
	    ul = (l1[k][2]+l1[k][3])/2;
	    cl = (l1[k][2]-l1[k][3])*(gamma-1)/4;
		rr = l1[k-1][1];
	    ur = (l1[k-1][2]+l1[k-1][3])/2;
	    cr = (l1[k-1][2]-l1[k-1][3])*(gamma-1)/4;
	    lambda = x-rl-dt*(ul-cl);
	    lambda = lambda/(rr-rl + dt*(ur-ul-cr+cl));
	    pm[0] = l1[0][0];
	    pm[1] = rl+lambda*(rr-rl);
	    pm[2] = l1[k][2]+lambda*(l1[k-1][2]-l1[k][2]);
	    pm[3] = l1[k][3]+lambda*(l1[k-1][3]-l1[k][3]);
	    if ((lambda>=0)and(lambda<=1))
	    	{
//	    	cout << "нашли правильное лямбда !!!" << endl;
	    	break;
	    	}
		}
	if ((lambda<=0)or(lambda>=1))
    	{
    	cout << "неправильное лямбда !!!" << endl;
    	}

	// считаем первое приближение
	double a = Cm(pm[2],pm[3]);

    pnp1[1] = x;
    pnp1[0] = pm[0]+(pnp1[1]-pm[1])/a;
    pnp1[3] = pm[3]+(pnp1[0]-pm[0])*Fm(pm[1],pm[2],pm[3]);
    pnp1[2] = 2*u-pnp1[3];
    double unp1 = (pnp1[2]+pnp1[3])/2;
    double cnp1 = (pnp1[2]-pnp1[3])*(gamma-1)/4;

    l2[N+1][0] = pnp1[0];
    l2[N+1][1] = pnp1[1];
    l2[N+1][2] = pnp1[2];
    l2[N+1][3] = pnp1[3];

//  ==== пересчитываем ========
	for (int k=count; k>N-1; k--)
		{
		rl = l1[k][1];
	    ul = (l1[k][2]+l1[k][3])/2;
	    cl = (l1[k][2]-l1[k][3])*(gamma-1)/4;
		rr = l1[k-1][1];
	    ur = (l1[k-1][2]+l1[k-1][3])/2;
	    cr = (l1[k-1][2]-l1[k-1][3])*(gamma-1)/4;
	    lambda = 2*(x-rl)-dt*(ul-cl+unp1-cnp1);
	    lambda = lambda/(2*(rr-rl) + dt*(ur-ul-cr+cl));
	    pm[0] = l1[0][0];
	    pm[1] = rl+lambda*(rr-rl);
	    pm[2] = l1[k][2]+lambda*(l1[k-1][2]-l1[k][2]);
	    pm[3] = l1[k][3]+lambda*(l1[k-1][3]-l1[k][3]);
	    if ((lambda>=0)and(lambda<=1))
	    	{
//	    	cout << "нашли правильное лямбда !!!" << endl;
	    	break;
	    	}
		}

	if ((lambda<=0)or(lambda>=1))
    	{
    	cout << "неправильное лямбда !!!" << endl;
    	}
    pnp1[3] = pm[3]+(pnp1[0]-pm[0])*(Fm(pm[1],pm[2],pm[3]) + Fm(pnp1[1],pnp1[2],pnp1[3]))/2;
    pnp1[2] = 2*u-pnp1[3];

    l2[N+1][0] = pnp1[0];
    l2[N+1][1] = pnp1[1];
    l2[N+1][2] = pnp1[2];
    l2[N+1][3] = pnp1[3];


//    return N+1;  // возвращаем наибольший номер точек на новом построенном слое
}

// ======================================================================================
// расчет погрешности масс
double dm(double rs, double rho0, double rfnew, double rho_, double rw)
{
	double m0, m_;
    if (nu == 0)
    	{
        m0 = (rw - rs)*rho0;
        m_ = (rw - rfnew)*rho_;
    	}
    else if (nu == 1)
		{
        m0 = 2*pi*(rw*rw - rs*rs)*rho0;
        m_ = 2*pi*(rw*rw - rfnew*rfnew)*rho_;
		}
    else if (nu == 2)
		{
        m0 = (4/3)*pi*(rw*rw*rw - rs*rs*rs)*rho0;
        m_ = (4/3)*pi*(rw*rw*rw - rfnew*rfnew*rfnew)*rho_;
		}
	return 100*abs(m0-m_)/m;
}
