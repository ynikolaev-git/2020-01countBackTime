//============================================================================
// Name        : prog3.cpp
// Author      : Nikolaev Yuri
// Version     : 01
// Copyright   :
// Description : ������� �� ������ ��������� ����� � �������� ����������� �������.
// ��������� ��� �������� ��������� �����, �������� �������� ������� ���������� ��� �������� ��������� ��������� ������
// � ������� ������� ������� � ���������� ������.
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const int maxSizeOfLayer = 1000000;
double gamma;
int nu;
double rho_;
int n0;      // ����� ����� �� 0-����
int n1;      // ����� ����� �� ��������� ��������� � ������ ������
int n2;      // ����� ����� �� ����-������������� ������� 3
double m;      // ����� ���������� ����
double rho0;    // �������������� ���������
double tf;
double r1;
double pi = 3.1415926535;
int cC = 0;                           // ����� ����� ������� ����� ����� � ����������� ���.�����
bool msg = 1;
double rw;
double r_;
double c_;
double t012;
double dt0;

struct pointOfNet {
	int s,c;
	double t,x,R,L;
};

vector<pointOfNet> Layer1(maxSizeOfLayer);
vector<pointOfNet> Layer2(maxSizeOfLayer);
vector<pointOfNet> piston;
pointOfNet p1, p2, p, pp;


//=======================================================================
// ������ ���������� r ������������ ������
double calc_rw() {
	double res;
    if (nu == 0) res = r1 + m/rho0;
    else if (nu == 1) res = pow(r1*r1 + m/(pi*rho0) , 1.0/2.0);
    else if (nu == 2) res = pow(r1*r1*r1 + (3*m)/(4*pi*rho0) , 1.0/3.0);
    return res;
}

//=======================================================================
// ������ ���������� r_
double calc_r_() {
	double res;
    if (nu == 0) res = rw - m/rho_;
    else if (nu == 1) res = pow(rw*rw - m/(rho_*pi) , 1.0/2.0);
    else if (nu == 2) res = pow(rw*rw*rw - (3*m)/(4*rho_*pi) , 1.0/3.0);
    return res;
}

//=======================================================================
// ������ ������������� ��������������
double Cp(double R, double L) {
    return R*(gamma+1)/4+L*(3-gamma)/4;
}

//=======================================================================
// ������ ������������� ��������������
double Cm(double R, double L) {
    return L*(gamma+1)/4+R*(3-gamma)/4;
}

//=======================================================================
double Fp(double x, double R, double L) {
    return (-nu)*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
double Fm(double x, double R, double L) {
    return nu*(R*R-L*L)*(gamma-1)/(8*x);
}

//=======================================================================
pointOfNet newPoint1(pointOfNet p1, pointOfNet p2) {
// ����������� ���� � ����� ������������� � �
// �� �1 ��������� �-, �� �2 �+
// ������� ������ �����������
	pointOfNet p;
	double a = Cm(p1.R,p1.L);
    double b = Cp(p2.R,p2.L);
    p.t = (a*p1.t-b*p2.t+p2.x-p1.x)/(a-b);
    p.x = p1.x+(p2.x-p1.x+b*(p1.t-p2.t))*a/(a-b);
    p.R = p2.R+(p.t-p2.t)*Fp(p2.x,p2.R,p2.L);
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
// �������� ���������� ��������
    a = (Cm(p1.R,p1.L)+Cm(p.R,p.L))/2;
    b = (Cp(p2.R,p2.L)+Cp(p.R,p.L))/2;
    p.t = (a*p1.t-b*p2.t+p2.x-p1.x)/(a-b);
    p.x = p1.x+(p2.x-p1.x+b*(p1.t-p2.t))*a/(a-b);
    p.R = p2.R+(p.t-p2.t)*(Fp(p2.x,p2.R,p2.L)+Fp(p.x,p.R,p.L))/2;
    p.L = p1.L+(p.t-p1.t)*(Fm(p1.x,p1.R,p1.L)+Fm(p.x,p.R,p.L))/2;
    return p;
}
//=======================================================================
pointOfNet newPoint2(const pointOfNet p1) {
// ����������� ����� �������������� � ������ r=rw
// ������� ������ �����������
	pointOfNet p;
    double a = Cm(p1.R,p1.L);
    p.x = rw;
    p.t = p1.t+(p.x-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
    p.R = -p.L;
// �������� ���������� ��������
    a = (Cm(p1.R,p1.L)+Cm(p.R,p.L))/2;
    p.t = p1.t+(rw-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*(Fm(p1.x,p1.R,p1.L)+Fm(p.x,p.R,p.L))/2;
    p.R = -p.L;
    return p;
}

//==========================================================================
bool cross(pointOfNet p1, pointOfNet p2, pointOfNet& p) {
// �1, �2 ���� ������������������ ����� ����� ������� �������� ��������������.
// � � ����� ���������� �� ������� ���������� �������� ������. ���� flag=true,
// �� � � ����� ����� ����������, ����� � �� ����������}
    bool flagI = 0;
    double l;
    double u=(p.R+p.L)/2;
    double a = p2.x-p1.x+u*(p1.t-p2.t);
    if (a!=0) l= (p.x+u*p1.t-u*p.t-p1.x)/a;
    else l = 10;
    if ((0<=l)and(l<=1)) {
        p.t = p1.t + l*(p2.t-p1.t);
        p.x = p1.x + l*(p2.x-p1.x);
        p.L = p1.L + l*(p2.L-p1.L);
        p.R = p1.R + l*(p2.R-p1.R);
//          print('(0<=l<=1)')
        flagI = 1;
    }
    return flagI;
}

//=======================================================================
// ������ ����������� ����
double dm() {
    double rp_ = piston[0].x;
    double rp0 = piston[piston.size()-1].x;
    double m0, m_;
    if (nu == 0) {
        m0 = (rw - rp0)*rho0;
        m_ = (rw - rp_)*rho_;
    	}
    else if (nu == 1) {
        m0 = 2*pi*(rw*rw - rp0*rp0)*rho0;
        m_ = 2*pi*(rw*rw - rp_*rp_)*rho_;
    	}
    else if (nu == 2) {
        m0 = (4.0/3.0)*pi*(rw*rw*rw - rp0*rp0*rp0)*rho0;
        m_ = (4.0/3.0)*pi*(rw*rw*rw - rp_*rp_*rp_)*rho_;
    	}
    return 100 * abs(m0-m_)/m;
}
//=======================================================================
// ������� ������ ������
double Ap() {
    double res0, res1, res = 0;
    int i, k = piston.size();
    if (nu == 0) {
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= (res1+res0)*(piston[i+1].x-piston[i].x)/(2*gamma);
        	}
    	}
    else if (nu == 1) {
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= pi*(res1+res0)*(piston[i+1].x*piston[i+1].x-piston[i].x*piston[i].x)/(2*gamma);
        	}
    	}
    else if (nu == 2)
    	{
        for (i=0; i<k-1;i++) {
            res1 = pow(((piston[i+1].R-piston[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res0 = pow(((piston[i].R-piston[i].L)/4)*(gamma-1), 2*gamma/(gamma-1));
            res -= pi*(4.0/3.0)*(res1+res0)*(piston[i+1].x*piston[i+1].x*piston[i+1].x
                           -piston[i].x*piston[i].x*piston[i].x)/(2*gamma);
        	}
    	}
    return res;
}

//=======================================================================
// ��������� � ������ �����
void printStart() {
    cout << "     ==== �������� ���� ====" << endl;
//    print('������� ����� � ' + str(l-1))
    cout << "���������� �������� ���� = " << gamma << endl;
    cout << "���������� ��������� ������� = " << nu << endl;
    cout << "���������� �� ������ ��������� �� ����������� ������ = " << rw << endl;
    cout << "���������� �� ������ ��������� �� ������ � ������ ������ = " << r_ << endl;
    cout << "��������� ������� ���� = " << rho_ << endl;
    cout << "�������� ����� � ������ ���� = " << c_ << endl;
}

//=======================================================================
// ��������� � ����� �����
void printEnd() {
    cout << "                     ==== ��������� ���� ====" << endl;
    cout << "����� ������ (���������������� ����������)= " << piston[piston.size()-1].x << endl;
    cout << "����� ������ (�����)= " << piston[piston.size()-1].t << endl;
    cout << "����������� ����= " << dm() << endl;
    cout << "����� ������ (���������������� ����������)= " << piston[0].x << endl;
    cout << "����� ������ (�����)= " << piston[0].t << endl;
    cout << "������ ������ = " << Ap() << endl;
//    print('=== ������������ ����� ��������, ������= ', int(time.time() - start_time))
}

//=======================================================================
// ������ �����
void mainFunction() {
	int i, j;
	double u;
//    minOnL = {};
    int numL = 0;    // ����� ���� (������������� ��������������)
    int numOnL = 0;  // ����� �� ���� (������ � ����� � ������������ - ������ �� ������ � ������ ������)

//##################### ������ ������� ���� (�� �^- ��������������)
    int k = 1;
    double t, x, R, L;
    t = tf;
    R = c_*2 / (gamma-1);
    L = -R;
    while (t > t012) {
    	x = r_-c_*(t-tf);
        Layer1[numOnL] = {0, numOnL, t, x, R, L};
        if (k < 100) t -= dt0*0.01;
        else t -= dt0;
        k += 1;
        numOnL += 1;
    	}
    x = r_-c_*(t012-tf);
    Layer1[numOnL] = {0, numOnL, t012, x, R, L};
//    minOnL[0] = 0;
    cout << "�������� 0-����, ������������ ����� ����� �� ����= " << numOnL << endl;

//############ ����� �������� ����� � ������� 1, 2
    double c = c_;
    cout << "c_ = " << c_ << endl;
    double dc1 = c_/n1;
    for (j=0; j<n1; j++) {      // j - ���������� ����
        c = c - dc1;
        numL += 1;
        u = 2*(c_-c)/(gamma-1);
        p = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        Layer2[0] = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        for (i=0; i<numOnL + j; i++) {      // i - ���������� ����� �� �����
            p1 = Layer2[i];
            p2 = Layer1[i+1];
            p = newPoint1(p1, p2);
            p.s = j+1;
            p.c = i+1;
            //Layer2[i+1] = {j+1, i+1, p.t, p.x, p.R, p.L};
            Layer2[i+1] = p;
            cC += 1;
        	}
        // ������ ����������� � ������  rw
        p = newPoint2(p);
//        Layer2[i+2] = {j+1, i+2, p.t, p.x, p.R, p.L};
        p.s = j+1;
        p.c = i+1;
        Layer2[i+1] = p;
        cC += 1;
        Layer1 = Layer2;
        rho0 = (p.R - p.L)*(gamma-1)/4;
        cout << "��������� ���� ����� " << j+1 << endl;
//        cout << "�������� ��������� c= " << rho0 << endl;
        if (rho0 <=1) {
            cout << "�������� ��������� c= " << rho0 << endl;
            break;
        	}
    	}
    numL = j + 1;    // ��������� ����������� ����
    numOnL = i + 1;    // ��������� ����� �� ������������� ����
    cout << "��������� ������� 1: ����� ����= " << numL << ", ������������ ����� ����� �� ����= " << numOnL << endl;


// ����� ������ ����� � ������� 2
    int p_numL = numL;
    int p_numOnL = 0;
    pp = Layer2[0];
    piston.push_back(pp);
    double t0 = p.t;
//    print ('����� � ����� D =', str(t0))
//    print ('��������� r � ����� D =', str(x))
//    #RLxt2 = RLxt[RLxt['s']==numL].reset_index(drop=True)
    double dt2 = (rw-r1)/n2;
    int cnt = 0;
    i = 0;
    while (p_numOnL < numOnL) { // ���������� ���� ���������� � �������� �������� �����
        Layer2[numOnL] = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
        j = numOnL;
//        minOnL[numL+1+i] = j
    //    Layer1 = RLxt2[RLxt2['s']==numL+i]
        p = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
        while (j > cnt) {   // ������ ���� � ������� ����� �������
            p1 = p;
            p2 = Layer1[j-1];
            p = newPoint1(p1, p2);
            p.s = numL+1+i;
            p.c = j-1;
            Layer2[j-1] = p;
            cC += 1;
//            minOnL[numL+1+i] = j-1;
            j -= 1;
    		}
        cout << "��������� ���� ����� " << numL+1+i << " �� ���������� ������ " << numOnL - j << " �����" << endl;
//       �������� ����������� ������ ������� ������ ������ �����
        bool flagDown = 1;
        bool flag;
        while ((flagDown==1)and(p_numOnL < numOnL)) {
            p1 = Layer2[p_numOnL+1];
            p2 = Layer1[p_numOnL+1];
            flag = 0;
            flag = cross(p1, p2, pp);
            if (flag==1) {
                p_numOnL += 1;
                pp.s = p1.s;
                pp.c = p_numOnL;
                piston.push_back(pp);
//                cout << "t= " << pp.t << "; x= " << pp.x << endl;
            	}
            else flagDown = 0;
        	}
//       �������� ����������� ����� ������� �������� ������
        if (p_numOnL == numOnL) {
            cout << "��������� ������ ���������� ������" << endl;
            cout << "������ ������ = " << piston.size() << endl;
            break;
        	}
        p1 = Layer2[p_numOnL];
        p2 = Layer2[p_numOnL+1];
        flag = 0;
        flag = cross(p1,  p2, pp);
        if (flag==1) {
            p_numL += 1;
            piston.push_back(pp);
            cnt = p_numOnL;
        	}
        else {
            cout << "��������� ����������, �������� ��������� �����" << endl;
            //print (p_numL, p_numOnL, t, x, R, L)
            break;
        	}
        Layer1 = Layer2;
        i += 1;
    	}
    numL += i;

    for (i=0;i<piston.size();++i)
    	piston[i].t -= piston[piston.size()-1].t;

    cout << "��������� ������� 2: ����� ����= " << numL << ", ������������ ����� ����� �� ����= " << numOnL << endl;
    cout << "����� �����= " << cC << endl;
}




//############################################
//    ### ������ ���� ��������� ����� ###
//############################################
int main() {
	// ������ excel-����
	//	wb = openpyxl.load_workbook('examples.xlsx')
	gamma = 1.6666667;
	nu = 2;
	rho_ = 1000;
	n0 = 10000;      // ����� ����� �� 0-����
	n1 = 10000;      // ����� ����� �� ��������� ��������� � ������ ������
	n2 = 10000;      // ����� ����� �� ����-������������� ������� 3
	m = 10;      // ����� ���������� ����
	rho0 = 1;    // �������������� ���������
	tf = 1;
	r1 = 1;
	rw = calc_rw();
	r_ = calc_r_();
	c_ = pow(rho_, (gamma-1)/2.0);
	t012 = tf - (rw - r_) / c_;
	dt0 = (tf - t012) / n0;

    if (msg==1) printStart();
    mainFunction();
    if (msg==1) printEnd();

    double u_pred, u, u_next, c_pred, c, c_next, t_s, t_f;
    t_s = piston[piston.size()-1].t;
    t_f = piston[0].t;
    for (int i=piston.size()-1;i>0;i--)
    	{
    	u_pred = (piston[i-1].R+piston[i-1].L)/2;
    	u = (piston[i].R+piston[i].L)/2;
    	u_next = (piston[i+1].R+piston[i+1].L)/2;
    	c_pred = (gamma-1)*(piston[i-1].R-piston[i-1].L)/4;
    	c = (gamma-1)*(piston[i].R-piston[i].L)/4;
    	c_next = (gamma-1)*(piston[i+1].R-piston[i+1].L)/4;
    	// ���� �������
    	if ((u_pred>u)&(u_next>u))  // ����� ������� ��������
    		cout << "����� min ��������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred>c)&(c_next>c))  // ����� ������� ���������
    		cout << "����� min ���������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	// ���� ��������
    	if ((u_pred<u)&(u_next<u))  // ����� �������� ��������
    		cout << "����� max ��������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred<c)&(c_next<c))  // ����� �������� ���������
    		cout << "����� max ���������: t=" << piston[i].t << ", r=" << piston[i].x <<
			", ����� ����� ������= " << i <<
			", ������ r ������ (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
    		", ������ t ������ (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	}
    cout << "������������ ����� ����� �� ���������� ������ = " << piston.size()-1 << endl;


//	for (int l=2; l<2; l++) {
//	    sheet = wb['---']   # �������� ����
//	    p = pd.DataFrame(columns = ['p_numL', 'p_numOnL', 't', 'x', 'R', 'L'])
//	    start_time = time.time()
//	    gamma = float(sheet['B'+str(l)].value) # ���������� ��������
//	    nu = int(sheet['C'+str(l)].value)      # ��� ���������
//	    rho_ = float(sheet['E'+str(l)].value)  # ��������� ������� ����
//	    n0 = int(sheet['F'+str(l)].value)      # ����� ����� �� 0-����
//	    n1 = int(sheet['G'+str(l)].value)      # ����� ����� �� ��������� ��������� � ������ ������
//	    n2 = int(sheet['H'+str(l)].value)      # ����� ����� �� ����-������������� ������� 3
//	    m = int(sheet['I'+str(l)].value)       # ����� ���������� ����
//	    c_ = pow(rho_, (gamma-1) / 2)
//	    rw = calc_rw();
//	    r_ = calc_r_();
//	    t012 = tf - (rw - r_) / c_
//	    dt0 = (tf - t012) / n0
//	    cC = 0                           # ����� ����� ������� ����� ����� � ����������� ���.�����
//	    mainFunction();
//	    if (msg==1) printEnd();
//	    print('�������='+str(l), ' & ', 'nu='+str(nu), ' & ', 'gamma='+str(gamma),
//	          ' & ', 'rho_='+str(rho_), ' & ', 'm='+str(m), '&', 'dm='+str(dm()), ' & - \\')
//	    print('- & ', p.loc[0].t, ' & ', p.loc[0].x, ' & ', Ap(), ' & - & - & - \\')
//		}

return 0;
}
