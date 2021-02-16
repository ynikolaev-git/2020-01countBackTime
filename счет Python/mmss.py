# -*- coding: utf-8 -*-
"""
Created on Sat Dec 14 21:08:52 2019

@author: 1
"""
import math
import pandas as pd
pd.set_option('precision', 16)
import time
import openpyxl

msg = 1      # Флажок - выводить сообщения в начале и в конце счета
#=======================================================================
# расчет координаты r неподвижного поршня
#def calc_rw(nu, m, r1, rho0):
def calc_rw():
    if nu == 0:
        res = r1 + m/rho0
    elif nu == 1:
        res = math.pow(r1*r1 + m/(math.pi*rho0) , 1 / 2)
    elif nu == 2:
        res = math.pow(r1*r1*r1 + (3*m)/(4*math.pi*rho0) , 1 / 3)
    return res;

#=======================================================================
# расчет координаты r_
def calc_r_():
    if nu == 0:
        res = rw - m/rho_
    elif nu == 1:
        res = math.pow(rw*rw - m/(rho_*math.pi) , 1 / 2)
    elif nu == 2:
        res = math.pow(rw*rw*rw - (3*m)/(4*rho_*math.pi) , 1 / 3)
    return res;

#=======================================================================
# наклон положительной характеристики
def Cp(R,L):
    return R*(gamma+1)/4+L*(3-gamma)/4;

#=======================================================================
# наклон отрицательной характеристики
def Cm(R,L):
    return L*(gamma+1)/4+R*(3-gamma)/4;

#=======================================================================
def Fp(x,R,L):
    return (-nu)*(R*R-L*L)*(gamma-1)/(8*x)

#=======================================================================
def Fm(x,R,L):
    return nu*(R*R-L*L)*(gamma-1)/(8*x)

#=======================================================================
def newPoint1(t1,x1,R1,L1, t2,x2,R2,L2):
# Пересечение плюс и минус характеристик в р 
# Из р1 выпускаем С-, из р2 С+ 
# считаем первое приближение
    a = Cm(R1,L1)
    b = Cp(R2,L2)
    t = (a*t1-b*t2+x2-x1)/(a-b)
    x = x1+(x2-x1+b*(t1-t2))*a/(a-b)
    R = R2+(t-t2)*Fp(x2,R2,L2)
    L = L1+(t-t1)*Fm(x1,R1,L1)
#{ уточняем полученные значения }
    a = (Cm(R1,L1)+Cm(R,L))/2
    b = (Cp(R2,L2)+Cp(R,L))/2
    t = (a*t1-b*t2+x2-x1)/(a-b)
    x = x1+(x2-x1+b*(t1-t2))*a/(a-b);
    R = R2+(t-t2)*(Fp(x2,R2,L2)+Fp(x,R,L))/2;
    L = L1+(t-t1)*(Fm(x1,R1,L1)+Fm(x,R,L))/2;
    return (t,x,R,L)

#=======================================================================
def newPoint2(t1, x1, R1, L1):
# Пересечение минус характеристики с прямой r=rw
# считаем первое приближение }
    a = Cm(R1,L1)
    x = rw
    t = t1+(x-x1)/a
    L = L1+(t-t1)*Fm(x1,R1,L1)
    R = -L
#{ уточняем полученные значения }
    a = (Cm(R1,L1)+Cm(R,L))/2
    t = t1+(rw-x1)/a
    L = L1+(t-t1)*(Fm(x1,R1,L1)+Fm(x,R,L))/2
    R = -L
    return (t,x,R,L)
    
#{==========================================================================}
def cross(t1, x1, R1, L1,  t2, x2, R2, L2, t, x, R, L):
#р1, р2 узлы характеристической сетки через которые проходит характеристика.
# В р точка траектории из которой продолжаем движение поршня. Если flag=true,
# то в р новая точка траектории, иначе р не изменяется}
    flagI = 0
    u=(R+L)/2
    a = x2-x1+u*(t1-t2)
    if a!=0:
        l = (x+u*t1-u*t-x1)/a
    else:
        l = 10
    if (0<=l<=1):
        t = t1 + l*(t2-t1)
        x = x1 + l*(x2-x1)
        L = L1 + l*(L2-L1)
        R = R1 + l*(R2-R1)
#        print('(0<=l<=1)')
        flagI = 1
    return (t,x,R,L,flagI)

#=======================================================================
# расчет погрешности масс
def dm():
    rp_ = p.loc[0].x
    rp0 = p.loc[len(p)-1].x
    if nu == 0:
        m0 = (rw - rp0)*rho0
        m_ = (rw - rp_)*rho_
    elif nu == 1:
        m0 = 2*math.pi*(rw*rw - rp0*rp0)*rho0
        m_ = 2*math.pi*(rw*rw - rp_*rp_)*rho_
    elif nu == 2:
        m0 = (4/3)*math.pi*(rw*rw*rw - rp0*rp0*rp0)*rho0
        m_ = (4/3)*math.pi*(rw*rw*rw - rp_*rp_*rp_)*rho_
    return 100 * abs(m0-m_)/m

#=======================================================================
# считаем работу поршня
def Ap():
    res = 0
    k = len(p)
    if nu == 0:
        for i in range (0,k-1):
            res1 = math.pow(((p.loc[i+1].R-p.loc[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res0 = math.pow(((p.loc[i].R-p.loc[i].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res -= (res1+res0)*(p.loc[i+1].x-p.loc[i].x)/(2*gamma)
    elif nu == 1:
        for i in range (0,k-1):
            res1 = math.pow(((p.loc[i+1].R-p.loc[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res0 = math.pow(((p.loc[i].R-p.loc[i].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res -= math.pi*(res1+res0)*(p.loc[i+1].x*p.loc[i+1].x-p.loc[i].x*p.loc[i].x)/(2*gamma)
    elif nu == 2:
        for i in range (0,k-1):
            res1 = math.pow(((p.loc[i+1].R-p.loc[i+1].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res0 = math.pow(((p.loc[i].R-p.loc[i].L)/4)*(gamma-1), 2*gamma/(gamma-1))
            res -= math.pi*(4/3)*(res1+res0)*(p.loc[i+1].x*p.loc[i+1].x*p.loc[i+1].x
                           -p.loc[i].x*p.loc[i].x*p.loc[i].x)/(2*gamma)
    return res

#=======================================================================
# сообщение в начале счета
def printStart():
    print('     ==== Начинаем счет ====')
    print('Вариант счета № ' + str(l-1))
    print('Показатель адиабаты газа = ' + str(gamma))
    print('Показатель симметрии течений = ' + str(int(nu)))
    print('Расстояние от центра симметрии до неподвижной стенки = ' + str(rw))
    print('Расстояние от центра симметрии до поршня в момент сжатия = ' + str(r_))
    print('Плотность сжатого газа = ' + str(rho_))
    print('Скорость звука в сжатом газе = ' + str(c_))

#=======================================================================
# сообщение в конце счета
def printEnd():
    print('                     ==== Завершили счет ====')
    print('Старт поршня (пространственная координата)= ', p.loc[len(p)-1].x)
    print('Старт поршня (время)= ', p.loc[len(p)-1].t)
    print('Погрешность масс= ', dm())
    print('Финал поршня (время)= ', p.loc[0].t)
    print('Финал поршня (пространственная координата)= ', p.loc[0].x)
    print('Работа поршня = ', Ap())
    print('=== Длительность счета варианта, секунд= ', int(time.time() - start_time))

#=======================================================================
# строим сетку
def mainFunction():
    global cC
    if msg:
        printStart()
    #RLxt = pd.DataFrame(columns = ['s', 'c', 't', 'x', 'R', 'L'])
    Layer1 = pd.DataFrame(columns = ['s', 'c', 't', 'x', 'R', 'L'])
    Layer2 = pd.DataFrame(columns = ['s', 'c', 't', 'x', 'R', 'L'])
    minOnL = {}
    t = tf
    numL = 0  # номер слоя (отрицательная характеристика)
    numOnL = 0  # номер на слое (начало в точке с особенностью - разрыв на поршне в момент стажития)

##################### строим нулевой слой (строку)
    k = 1
    R = c_*2 / (gamma-1)
    L = -R
    while (t > t012):
        Layer1.loc[numOnL] = [0, numOnL, t, r_-c_*(t-tf), R, L]
        if k < 100:
            t -= dt0*0.01
        else:
            t -= dt0
        k += 1
        numOnL += 1
    Layer1.loc[numOnL] = [0, numOnL, t012, r_-c_*(t012-tf), R, L]
    minOnL[0] = 0
    print('построен 0-слой, максимальный номер точки на слое= ', numOnL)

############ далее построим сетку в области 1
    c = c_
    print('c_ = ', c_)
    dc1 = c_/n1
    minOnL[numL] = 0
    for j in range (0,n1):      # j - перебираем слои
        c = c - dc1
        numL += 1
        u = 2*(c_-c)/(gamma-1)
        st, ct, t, x, R, L = numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)
        Layer2.loc[0] = [numL, 0, tf, r_, u + c*2 / (gamma-1), u - c*2 / (gamma-1)]
    ##    RLxt.loc[len(RLxt)] = [s, c, t, x, R, L]
        for i in range(0, numOnL + j):      # i - перебираем точки на слоях
            t1, x1, R1, L1 = t, x, R, L
            st, ct, t2, x2, R2, L2 = Layer1.loc[(Layer1['c']==i+1)].iloc[0].tolist()
            t, x, R, L = newPoint1(t1, x1, R1, L1, t2, x2, R2, L2)
            Layer2.loc[i+1] = [j+1, i+1, t, x, R, L]
            cC += 1
    ##        RLxt.loc[len(RLxt)] = [j+1, i+1, t, x, R, L]
        # строим переечение с линией  rw
        t, x, R, L = newPoint2(t, x, R, L)
        Layer2.loc[i+2] = [j+1, i+2, t, x, R, L]
        cC += 1
    ##    RLxt.loc[len(RLxt)] = [j+1, i+2, t, x, R, L]
        Layer1 = Layer2.reset_index(drop=True)
        minOnL[j+1] = 0
        print('построили слой номер ', j+1)
        rho0 = (R - L)*(gamma-1)/4
        if rho0 <=1:
            print('Достигли плотности c= ', rho0)
            break
    numL = j + 1    # последняя построенный слой
    numOnL = i + 2    # последняя точка на отрицательном слое
    print('построена область 1:', ' номер слоя= ', numL, ', максимальный номер точки на слое= ', numOnL)


# далее строим сетку в области 2
    p_numL = numL
    p_numOnL = 0
    st, ct, tp, xp, Rp, Lp = Layer2.loc[Layer2['c']==0].iloc[0].tolist()
    p.loc[len(p)] = [p_numL, p_numOnL, tp, xp, Rp, Lp]
    t0 = t
    print ('Время в точке D =', str(t0))
    print ('Координта r в точке D =', str(x))
    #RLxt2 = RLxt[RLxt['s']==numL].reset_index(drop=True)
    dt2 = (rw-r1)/n2
    cnt = 0
    i = 0
    while p_numOnL < numOnL: # перебираем слои приходящие в исходный несжатый покой
        Layer2.loc[numOnL] = [numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)]
        j = numOnL
        minOnL[numL+1+i] = j
    #    Layer1 = RLxt2[RLxt2['s']==numL+i]
        st, ct, t, x, R, L = numL+1+i, numOnL, t0 - dt2*(i+1), rw - dt2*(i+1), 2 / (gamma-1), 2 / (1-gamma)
        while (j > cnt):   # строим слой в сторону роста времени
            t1, x1, R1, L1 = t, x, R, L
            st, ct, t2, x2, R2, L2 = Layer1.loc[Layer1['c']==j-1].iloc[0].tolist()
            t, x, R, L = newPoint1(t1, x1, R1, L1, t2, x2, R2, L2)
            Layer2.loc[j-1] = [numL+1+i, j-1, t, x, R, L]
            cC += 1
            minOnL[numL+1+i] = j-1
            j -= 1
        print('построили слой номер ', numL+1+i, ' до траектории поршня', numOnL - j,'узлов')
#       проверим пересечения нижних боковый сторон ячейки сетки
        flagDown = 1
        while (flagDown)and(p_numOnL < numOnL):
            st, ct, t1, x1, R1, L1 = Layer2.loc[Layer2['c']==p_numOnL+1].iloc[0].tolist()
            st, ct, t2, x2, R2, L2 = Layer1.loc[Layer1['c']==p_numOnL+1].iloc[0].tolist()
            flag = 0
            tp, xp, Rp, Lp, flag = cross(t1, x1, R1, L1,  t2, x2, R2, L2, tp, xp, Rp, Lp)
            if flag:
                p_numOnL += 1
                p.loc[len(p)] = [p_numL, p_numOnL, tp, xp, Rp, Lp]
            else:
                flagDown = 0
#       проверим пересечение левой боковой стороной ячейки
        if p_numOnL == numOnL:
            print('Завершили расчет траектории поршня')
            break
        st, ct, t1, x1, R1, L1 = Layer2.loc[Layer2['c']==p_numOnL].iloc[0].tolist()
        st, ct, t2, x2, R2, L2 = Layer2.loc[Layer2['c']==p_numOnL+1].iloc[0].tolist()
        flag = 0
        tp, xp, Rp, Lp, flag = cross(t1, x1, R1, L1,  t2, x2, R2, L2, tp, xp, Rp, Lp)
        if flag:
            p_numL += 1
            p.loc[len(p)] = [p_numL, p_numOnL, tp, xp, Rp, Lp]
            cnt = p_numOnL
        else:
            print('Аварийное завершение, измените параметры счета')
            print (p_numL, p_numOnL, t, x, R, L)
            break
        Layer1 = Layer2.reset_index(drop=True)
        Layer2[:] = 0
        i += 1
    numL += i
    
    p.t -= p.loc[len(p)-1].t
    
    print('построена область 2:', ' номер слоя= ', numL, ', максимальный номер точки на слое= ', numOnL)
    print('Число ячеек= ', cC)
    

############################################
    ### дальше сама программа счета ###
############################################
rho0 = 1       #первлначальная плотность
tf = 1
r1 = 1  # 

# читаем excel-файл
wb = openpyxl.load_workbook('examples.xlsx')

for l in range(2,3):
    sheet = wb['---']   # получаем лист
    p = pd.DataFrame(columns = ['p_numL', 'p_numOnL', 't', 'x', 'R', 'L'])
    start_time = time.time()
    gamma = float(sheet['B'+str(l)].value) # показатель адиабаты
    nu = int(sheet['C'+str(l)].value)      # вид симметрии
    rho_ = float(sheet['E'+str(l)].value)  # плотность сжатого газа
    n0 = int(sheet['F'+str(l)].value)      # число точек на 0-слое
    n1 = int(sheet['G'+str(l)].value)      # число точек на ступеньке плотности в момент сжатия
    n2 = int(sheet['H'+str(l)].value)      # число точек на плюс-харктеристике области 3
    m = int(sheet['I'+str(l)].value)       # масса сжимаемого газа
    c_ = pow(rho_, (gamma-1) / 2)
    rw = calc_rw();
    r_ = calc_r_();
    t012 = tf - (rw - r_) / c_
    dt0 = (tf - t012) / n0
    cC = 0                           # здесь будет сколько всего ячеек в построенной хар.сетке
    mainFunction()
    if msg:
        printEnd()
    print('вариант='+str(l), ' & ', 'nu='+str(nu), ' & ', 'gamma='+str(gamma),
          ' & ', 'rho_='+str(rho_), ' & ', 'm='+str(m), '&', 'dm='+str(dm()), ' & - \\')
    print('- & ', p.loc[0].t, ' & ', p.loc[0].x, ' & ', Ap(), ' & - & - & - \\')
    sheet['L'+str(l)] = dm()
    sheet['M'+str(l)] = Ap()
    sheet['N'+str(l)] = p.loc[0].t
    sheet['O'+str(l)] = p.loc[0].x
    sheet['P'+str(l)] = '-'
    sheet['Q'+str(l)] = '-'
    sheet['R'+str(l)] = cC
    wb.remove(wb[str(l-1)])
    wb.create_sheet(title = str(l-1), index = l)
    sheet = wb[str(l-1)]       # получаем лист
    sheet.append(['№ point','t', 'r', 'R', 'L'])
    for j in range(0,len(p)):
        sheet.append([j, p.loc[j].t, p.loc[j].x, p.loc[j].R, p.loc[j].L])
    sheet = wb['forArticle']   # получаем лист
    sheet['A1'] = '№ варианта'
    sheet['B1'] = 'В таблицу 1'
    sheet['C1'] = 'В таблицу 2'
    sheet['A'+str(l)] = l-1
    sheet['B'+str(l)] = ('вариант='+str(l-1)
        +' & '+'nu='+str(nu)
        +' & '+'gamma='+str(gamma)
        +' & '+'rho_='+str(rho_)
        +' & '+'m='+str(m)
        +'&'+'dm='+str(dm())
        +' & - \\')
    sheet['C'+str(l)] = ('вариант='+str(l-1)
        +' & '+'tf='+str(p.loc[0].t)
        +' & '+'rf='+str(p.loc[0].x)
        +' & '+'Aopt='+str(Ap())
        +' & - & - & - \\')

# записываем файл
wb.save('examples.xlsx')

#import sys
#sys.exit(0)
#exit
    
