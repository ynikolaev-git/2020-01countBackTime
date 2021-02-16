# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 07:09:33 2019
@author: 1
"""

############################################
    ### дальше сама программа счета ###
############################################

#mainFunction()
#print(p.to_string())
#sys.exit(0)
#exit
from mmss import *
cntCells = 0   # здесь будет сколько всего ячеек в построенной хар.сетке
#import mmss


import openpyxl
# читаем excel-файл
wb = openpyxl.load_workbook('examples.xlsx')
sheet = wb['---']   # получаем лист

for i in range(2,3):
    gamma = float(sheet['B'+str(i)].value) # показатель адиабаты
    nu = int(sheet['C'+str(i)].value)      # вид симметрии
    rho0 = 1                               # плотность несжатого газа
    rho_ = float(sheet['E'+str(i)].value)  # плотность сжатого газа
    n0 = int(sheet['F'+str(i)].value)      # число точек на 0-слое
    n1 = int(sheet['G'+str(i)].value)      # число точек на ступеньке плотности в момент сжатия
    n2 = int(sheet['H'+str(i)].value)      # число точек на плюс-харктеристике области 3
    m = int(sheet['I'+str(i)].value)       # масса сжимаемого газа
#    c_ = pow(rho_, (gamma-1) / 2)
#    rw = calc_rw(nu, m, r1, rho0);
#    r_ = calc_r_(nu, m, rho_, rw);
#    t012 = tf - (rw - r_) / c_
#    dt0 = (tf - t012) / n0
    mainFunction(gamma, nu, rho_, n0, n1, n2, m)
    print ('старт поршня (пространственная координата)= ', p.loc[len(p)-1].x)
    print ('старт поршня (время)= ', p.loc[len(p)-1].t)
    print ('погрешность масс= ', dm(nu, rho_, p))
    print ('финал поршня (время)= ', p.loc[0].t)
    print ('финал поршня (пространственная координата)= ', p.loc[0].x)
    print ('работа поршня = ', Ap(nu, p))
    print('=== Длительность выполнения программыы, секунд= ', int(time.time() - start_time))
    print ('- & ', p.loc[0].t, ' & ', p.loc[0].x, ' & ', Ap(nu, p), ' & - & - & - \\')
    sheet = wb[str(i)]
    sheet.append(['t', 'r', 'R', 'L'])
    for j in range(0,len(p)):
        sheet.append([j, p.loc[j].t, p.loc[j].x, p.loc[j].R, p.loc[j].L])

# записываем файл
wb.save('examples.xlsx')

#mainFunction()
#nu = 0
#mainFunction()

#sys.exit(0)
#exit
    
