/*
 * mx.h
 *
 *  Created on: 3 мая 2020 г.
 *      Author: 1
 */

#ifndef MX_HPP_
#define MX_HPP_

extern double rho0;
extern double c0;
extern double gamma;
extern double rho_;
extern double m;
extern double pi;
extern double ts, tf, rf, rs, rw, maxNum, u_max, c_max, rho_max;

extern int nu;
extern const int maxLenPiston;
extern const int maxLenLayer;
extern double pistonEnd;
extern int numOnPiston;

extern double Cp(double R,double L);
extern double Cm(double R,double L);
extern double Fp(double x, double R, double L);
extern double Fm(double x, double R, double L);
//extern void newPoint1(double* p1, double* p2, double* p3);
extern void newPoint1(double p1[4], double p2[4], double p3[4]);
extern void newPoint2(double* p1, double* p2);
//extern bool newPoint3_prog2(double* p1, double* p2, double* p3,  double* p4);
extern bool newPoint3_prog2(double p1[4], double p2[4], double p3[4],  double p4[4]);
extern bool newPoint1_prog3(double l1[10000][4], double l2[10000][4], double tLayer, int numPoint);
extern void newPoint2_prog3(double l1[10000][4], double l2[10000][4], double tLayer);
extern void newPoint3_prog3(double l1[10000][4], double l2[10000][4], double piston[10000][4], double tL, int count);
extern double calc_rw();
extern void readPiston(const char* fName, double p[200000][4]);
double dm(double rs, double rho0, double rfnew, double rho_, double rw);

#endif /* MX_HPP_ */
