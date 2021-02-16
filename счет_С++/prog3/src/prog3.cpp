//============================================================================
// Name        : prog3.cpp
// Author      : Nikolaev Yuri
// Version     : 01
// Copyright   :
// Description : Перенос из Питона программы счета в обратном направлении времени.
// Программа для быстрого массового счета, плавного перебора входных параметров для изучения характера поведения поршня
// в моменты времени близкие к финальному сжатию.
//============================================================================

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

const int maxSizeOfLayer = 1000000;
double gamma;
int nu;
double rho_;
int n0;      // число точек на 0-слое
int n1;      // число точек на ступеньке плотности в момент сжатия
int n2;      // число точек на плюс-харктеристике области 3
double m;      // масса сжимаемого газа
double rho0;    // первоначальная плотность
double tf;
double r1;
double pi = 3.1415926535;
int cC = 0;                           // здесь будет сколько всего ячеек в построенной хар.сетке
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
// расчет координаты r неподвижного поршня
double calc_rw() {
	double res;
    if (nu == 0) res = r1 + m/rho0;
    else if (nu == 1) res = pow(r1*r1 + m/(pi*rho0) , 1.0/2.0);
    else if (nu == 2) res = pow(r1*r1*r1 + (3*m)/(4*pi*rho0) , 1.0/3.0);
    return res;
}

//=======================================================================
// расчет координаты r_
double calc_r_() {
	double res;
    if (nu == 0) res = rw - m/rho_;
    else if (nu == 1) res = pow(rw*rw - m/(rho_*pi) , 1.0/2.0);
    else if (nu == 2) res = pow(rw*rw*rw - (3*m)/(4*rho_*pi) , 1.0/3.0);
    return res;
}

//=======================================================================
// наклон положительной характеристики
double Cp(double R, double L) {
    return R*(gamma+1)/4+L*(3-gamma)/4;
}

//=======================================================================
// наклон отрицательной характеристики
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
// Пересечение плюс и минус характеристик в р
// Из р1 выпускаем С-, из р2 С+
// считаем первое приближение
	pointOfNet p;
	double a = Cm(p1.R,p1.L);
    double b = Cp(p2.R,p2.L);
    p.t = (a*p1.t-b*p2.t+p2.x-p1.x)/(a-b);
    p.x = p1.x+(p2.x-p1.x+b*(p1.t-p2.t))*a/(a-b);
    p.R = p2.R+(p.t-p2.t)*Fp(p2.x,p2.R,p2.L);
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
// уточняем полученные значения
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
// Пересечение минус характеристики с прямой r=rw
// считаем первое приближение
	pointOfNet p;
    double a = Cm(p1.R,p1.L);
    p.x = rw;
    p.t = p1.t+(p.x-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*Fm(p1.x,p1.R,p1.L);
    p.R = -p.L;
// уточняем полученные значения
    a = (Cm(p1.R,p1.L)+Cm(p.R,p.L))/2;
    p.t = p1.t+(rw-p1.x)/a;
    p.L = p1.L+(p.t-p1.t)*(Fm(p1.x,p1.R,p1.L)+Fm(p.x,p.R,p.L))/2;
    p.R = -p.L;
    return p;
}

//==========================================================================
bool cross(pointOfNet p1, pointOfNet p2, pointOfNet& p) {
// р1, р2 узлы характеристической сетки через которые проходит характеристика.
// В р точка траектории из которой продолжаем движение поршня. Если flag=true,
// то в р новая точка траектории, иначе р не изменяется}
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
// расчет погрешности масс
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
// считаем работу поршня
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
// сообщение в начале счета
void printStart() {
    cout << "     ==== Начинаем счет ====" << endl;
//    print('Вариант счета № ' + str(l-1))
    cout << "Показатель адиабаты газа = " << gamma << endl;
    cout << "Показатель симметрии течений = " << nu << endl;
    cout << "Расстояние от центра симметрии до неподвижной стенки = " << rw << endl;
    cout << "Расстояние от центра симметрии до поршня в момент сжатия = " << r_ << endl;
    cout << "Плотность сжатого газа = " << rho_ << endl;
    cout << "Скорость звука в сжатом газе = " << c_ << endl;
}

//=======================================================================
// сообщение в конце счета
void printEnd() {
    cout << "                     ==== Завершили счет ====" << endl;
    cout << "Старт поршня (пространственная координата)= " << piston[piston.size()-1].x << endl;
    cout << "Старт поршня (время)= " << piston[piston.size()-1].t << endl;
    cout << "Погрешность масс= " << dm() << endl;
    cout << "Финал поршня (пространственная координата)= " << piston[0].x << endl;
    cout << "Финал поршня (время)= " << piston[0].t << endl;
    cout << "Работа поршня = " << Ap() << endl;
//    print('=== Длительность счета варианта, секунд= ', int(time.time() - start_time))
}

//=======================================================================
// строим сетку
void mainFunction() {
	int i, j;
	double u;
//    minOnL = {};
    int numL = 0;    // номер слоя (отрицательная характеристика)
    int numOnL = 0;  // номер на слое (начало в точке с особенностью - разрыв на поршне в момент сжатия)

//##################### строим нулевой слой (на С^- характеристике)
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
    cout << "построен 0-слой, максимальный номер точки на слое= " << numOnL << endl;

//############ далее построим сетку в области 1, 2
    double c = c_;
    cout << "c_ = " << c_ << endl;
    double dc1 = c_/n1;
    for (j=0; j<n1; j++) {      // j - перебираем слои
        c = c - dc1;
        numL += 1;
        u = 2*(c_-c)/(gamma-1);
        p = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        Layer2[0] = {numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)};
        for (i=0; i<numOnL + j; i++) {      // i - перебираем точки на слоях
            p1 = Layer2[i];
            p2 = Layer1[i+1];
            p = newPoint1(p1, p2);
            p.s = j+1;
            p.c = i+1;
            //Layer2[i+1] = {j+1, i+1, p.t, p.x, p.R, p.L};
            Layer2[i+1] = p;
            cC += 1;
        	}
        // строим пересечение с линией  rw
        p = newPoint2(p);
//        Layer2[i+2] = {j+1, i+2, p.t, p.x, p.R, p.L};
        p.s = j+1;
        p.c = i+1;
        Layer2[i+1] = p;
        cC += 1;
        Layer1 = Layer2;
        rho0 = (p.R - p.L)*(gamma-1)/4;
        cout << "построили слой номер " << j+1 << endl;
//        cout << "Достигли плотности c= " << rho0 << endl;
        if (rho0 <=1) {
            cout << "Достигли плотности c= " << rho0 << endl;
            break;
        	}
    	}
    numL = j + 1;    // последний построенный слой
    numOnL = i + 1;    // последняя точка на отрицательном слое
    cout << "построена область 1: номер слоя= " << numL << ", максимальный номер точки на слое= " << numOnL << endl;


// далее строим сетку в области 2
    int p_numL = numL;
    int p_numOnL = 0;
    pp = Layer2[0];
    piston.push_back(pp);
    double t0 = p.t;
//    print ('Время в точке D =', str(t0))
//    print ('Координта r в точке D =', str(x))
//    #RLxt2 = RLxt[RLxt['s']==numL].reset_index(drop=True)
    double dt2 = (rw-r1)/n2;
    int cnt = 0;
    i = 0;
    while (p_numOnL < numOnL) { // перебираем слои приходящие в исходный несжатый покой
        Layer2[numOnL] = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
        j = numOnL;
//        minOnL[numL+1+i] = j
    //    Layer1 = RLxt2[RLxt2['s']==numL+i]
        p = {numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)};
        while (j > cnt) {   // строим слой в сторону роста времени
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
        cout << "построили слой номер " << numL+1+i << " до траектории поршня " << numOnL - j << " узлов" << endl;
//       проверим пересечения нижних боковых сторон ячейки сетки
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
//       проверим пересечение левой боковой стороной ячейки
        if (p_numOnL == numOnL) {
            cout << "Завершили расчет траектории поршня" << endl;
            cout << "Размер поршня = " << piston.size() << endl;
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
            cout << "Аварийное завершение, измените параметры счета" << endl;
            //print (p_numL, p_numOnL, t, x, R, L)
            break;
        	}
        Layer1 = Layer2;
        i += 1;
    	}
    numL += i;

    for (i=0;i<piston.size();++i)
    	piston[i].t -= piston[piston.size()-1].t;

    cout << "построена область 2: номер слоя= " << numL << ", максимальный номер точки на слое= " << numOnL << endl;
    cout << "Число ячеек= " << cC << endl;
}




//############################################
//    ### дальше сама программа счета ###
//############################################
int main() {
	// читаем excel-файл
	//	wb = openpyxl.load_workbook('examples.xlsx')
	gamma = 1.6666667;
	nu = 2;
	rho_ = 1000;
	n0 = 10000;      // число точек на 0-слое
	n1 = 10000;      // число точек на ступеньке плотности в момент сжатия
	n2 = 10000;      // число точек на плюс-харктеристике области 3
	m = 10;      // масса сжимаемого газа
	rho0 = 1;    // первоначальная плотность
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
    	// ищем минимум
    	if ((u_pred>u)&(u_next>u))  // нашли минимум скорости
    		cout << "точка min скорости: t=" << piston[i].t << ", r=" << piston[i].x <<
			", номер точки поршня= " << i <<
			", прошло r сжатия (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", прошло t сжатия (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred>c)&(c_next>c))  // нашли минимум плотности
    		cout << "точка min плотности: t=" << piston[i].t << ", r=" << piston[i].x <<
			", номер точки поршня= " << i <<
			", прошло r сжатия (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", прошло t сжатия (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	// ищем максимум
    	if ((u_pred<u)&(u_next<u))  // нашли максимум скорости
    		cout << "точка max скорости: t=" << piston[i].t << ", r=" << piston[i].x <<
			", номер точки поршня= " << i <<
			", прошло r сжатия (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
			", прошло t сжатия (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	if ((c_pred<c)&(c_next<c))  // нашли максимум плотности
    		cout << "точка max плотности: t=" << piston[i].t << ", r=" << piston[i].x <<
			", номер точки поршня= " << i <<
			", прошло r сжатия (%) = " << 100*(piston[i].x-piston[piston.size()-1].x)/(piston[0].x-piston[piston.size()-1].x) <<
    		", прошло t сжатия (%) = " << 100*(piston[i].t-t_s)/(t_f-t_s) << endl;
    	}
    cout << "Максимальный номер точки на траектории поршня = " << piston.size()-1 << endl;


//	for (int l=2; l<2; l++) {
//	    sheet = wb['---']   # получаем лист
//	    p = pd.DataFrame(columns = ['p_numL', 'p_numOnL', 't', 'x', 'R', 'L'])
//	    start_time = time.time()
//	    gamma = float(sheet['B'+str(l)].value) # показатель адиабаты
//	    nu = int(sheet['C'+str(l)].value)      # вид симметрии
//	    rho_ = float(sheet['E'+str(l)].value)  # плотность сжатого газа
//	    n0 = int(sheet['F'+str(l)].value)      # число точек на 0-слое
//	    n1 = int(sheet['G'+str(l)].value)      # число точек на ступеньке плотности в момент сжатия
//	    n2 = int(sheet['H'+str(l)].value)      # число точек на плюс-харктеристике области 3
//	    m = int(sheet['I'+str(l)].value)       # масса сжимаемого газа
//	    c_ = pow(rho_, (gamma-1) / 2)
//	    rw = calc_rw();
//	    r_ = calc_r_();
//	    t012 = tf - (rw - r_) / c_
//	    dt0 = (tf - t012) / n0
//	    cC = 0                           # здесь будет сколько всего ячеек в построенной хар.сетке
//	    mainFunction();
//	    if (msg==1) printEnd();
//	    print('вариант='+str(l), ' & ', 'nu='+str(nu), ' & ', 'gamma='+str(gamma),
//	          ' & ', 'rho_='+str(rho_), ' & ', 'm='+str(m), '&', 'dm='+str(dm()), ' & - \\')
//	    print('- & ', p.loc[0].t, ' & ', p.loc[0].x, ' & ', Ap(), ' & - & - & - \\')
//		}

return 0;
}
