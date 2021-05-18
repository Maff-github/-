import numpy as np
import sympy as sp
from scipy.optimize import fsolve
from math import sin, cos, asin, acos, pi
from numpy import exp

#Вводим экспоненту, отражающую изменение фазы элемента антенны с круговой поляризацией
def e (beta):
    return exp(-1j*beta)
#по условию beta=pi/2
beta1=pi/2
exp1=e(beta1)
print('exp=',exp1)

#k- любое натуральное число, которое используется в неопределенностях
k=1
d_small=0.000001

th1=float(input ('Введите Тета(от 0 до 180 градусов)='))
th3=float(input ('Введите Тетаmn='))
fi1=float(input ('Введите Фи(от 0 до 360 градусов)='))
fi2=float(input ('Введите Фиmn='))
print()

th2=acos(sin(th3*pi/180)*sin(th1*pi/180)*
              cos((fi1-fi2)*pi/180)+cos(th3*pi/180)*
              cos(th1*pi/180))

print('Тета=',th1)
print('Тета"=',th2)
print('Тетаmn"=',th3)
print('Фи=',fi1)
print('Фиmn"=',fi2)
print()

#p>=1
p=1
#f1(theta,fi)= f2(theta,fi)
#Вводим нормированную действительную амплитуду 1 (ДН элемента 1)
def func_f1(theta1,fi):
    if th2>=(-pi/2) and th2<=pi/2:
        print('-pi/2<= Тета" <=pi/2')
        return cos(theta1)**p
    else:
        return 0
f1=func_f1(th2,fi1)
print('Нормированная ДН элемента f1(theta,fi)=',f1)
print()

#Вводим нормированную действительную амплитуду 2 (ДН элемента 2)
f2=func_f1(th2,fi1)
print('Нормированная ДН элемента f2(theta,fi)=',f2)
print()

#Вводим функцию Amn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - Amn считается без приращений.
def func_Amn(theta,theta1,thetamn,fi,fimn):
    try:
        Amn=(cos(theta*pi/180)*sin(thetamn*pi/180)
             *cos((fi-fimn)*pi/180)-sin(theta*pi/180)
             *cos(thetamn*pi/180))/sin(theta1*pi/180)
        print('Избежали попадания в неопределенность, Amn считается без приращений')
        return Amn         
    except ZeroDivisionError:
        print('Непределенность 0/0. Приращаем Тета, Тета",(Фи-Фиmn), Тетаmn')
        print('Ищем значения слева и справа, после чего берем среднее между ними')
        Amn_left=(cos((theta-d_small)*pi/180)*sin((thetamn-d_small)*pi/180)
        *cos(((fi-fimn)-d_small)*pi/180)-sin((theta-d_small)*pi/180)
        *cos((thetamn-d_small)*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left=',Amn_left)
        Amn_right=(cos((theta+d_small)*pi/180)*sin((thetamn+d_small)*pi/180)
        *cos(((fi-fimn)+d_small)*pi/180)-sin((theta+d_small)*pi/180)
        *cos((thetamn+d_small)*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right=',Amn_right)
        Amn=(Amn_left+Amn_right)/2
        return Amn
a=func_Amn(th1,th2,th3,fi1,fi2)
print('Alfa=',a)
print()
            
#Вводим функцию Bmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - Bmn считается без приращений.    
def func_Bmn (theta,theta1,thetamn,fi,fimn):
    try:
        Bmn=sin(thetamn*pi/180)*sin((fi-fimn)*pi/180)/sin(theta1*pi/180)
        print('Избежали попадания в неопределенность, Bmn считается без приращений')
        return Bmn      
    except ZeroDivisionError:        
        print('Непределенность 0/0. Приращаем Тета, Тета",(Фи-Фиmn), Тетаmn')
        print('Ищем значения слева и справа, после чего берем среднее между ними')
        Bmn_left=sin((thetamn-d_small)*pi/180)*sin(((fi-fimn)-d_small)
            *pi/180)/sin((theta1-d_small)*pi/180)
        print('Beta_left1=',Bmn_left)
        Bmn_right=sin((thetamn+d_small)*pi/180)*sin(((fi-fimn)+d_small)
            *pi/180)/sin((theta1+d_small)*pi/180)
        print('Beta_right1=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2
        return Bmn     
b=func_Bmn(th1,th2,th3,fi1,fi2)
print('Beta=',b)
print()

#Вводим функцию ctg проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - ctg считается без приращений.  
def func_ctgfi1 (theta,thetamn,fi,fimn):
    try:
        ctgfi1=(cos(thetamn*pi/180)*sin(theta*pi/180)
                *cos((fi-fimn)*pi/180)-sin(thetamn*pi/180)
                *cos(theta*pi/180))/(sin(theta*pi/180)*sin((fi-fimn)*pi/180)) 
        print('Избежали попадания в неопределенность, ctg считается без приращений')
        return ctgfi1
    except ZeroDivisionError:
        print('Непределенность 0/0. Приращаем Тета,(Фи-Фиmn), Тетаmn')
        print('Ищем значения слева и справа, после чего берем среднее между ними')
        ctgfi1_left=(cos((thetamn-d_small)*pi/180)*sin((theta-d_small)*pi/180)
                     *cos(((fi-fimn)-d_small)*pi/180)-sin((thetamn-d_small)
                     *pi/180)*cos((theta-d_small)*pi/180))/(sin((theta-d_small)*pi/180)
                     *sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left=',ctgfi1_left)
        ctgfi1_right=(cos((thetamn+d_small)*pi/180)*sin((theta+d_small)*pi/180)
                     *cos(((fi-fimn)+d_small)*pi/180)-sin((thetamn+d_small)
                     *pi/180)*cos((theta+d_small)*pi/180))/(sin((theta+d_small)*pi/180)
                     *sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right6=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2        
        return ctgfi1                
ctgf=func_ctgfi1(th1,th3,fi1,fi2)
print('CtgFi"=',ctgf)
print()

#Вводим функцию tgYmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - tgYmn считается без приращений.  
def func_tgYmn (theta,thetamn,fi,fimn):
    try:
        tgYmn=((cos(thetamn*pi/180)+cos(theta*pi/180))*sin((fi-fimn)
                *pi/180))/(sin(thetamn*pi/180)*sin(theta*pi/180)+(cos(theta
                *pi/180)*cos(thetamn*pi/180)+1)*cos((fi-fimn)*pi/180))
        print('Избежали попадания в неопределенность, tgYmn  считается без приращений')
        return tgYmn
    except ZeroDivisionError:
        print('Непределенность 0/0. Приращаем Тета,(Фи-Фиmn), Тетаmn')
        print('Ищем значения слева и справа, после чего берем среднее между ними')
        tgYmn_left=((cos((thetamn-d_small)*pi/180)+cos((theta-d_small)*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin((thetamn-d_small)*pi/180)
                    *sin((theta-d_small)*pi/180)+(cos((theta-d_small)*pi/180)
                    *cos((thetamn-d_small)*pi/180)+1)*cos(((fi-fimn)-d_small)*pi/180))
        print('tgYmn_left=',tgYmn_left)
        tgYmn_right=((cos((thetamn+d_small)*pi/180)+cos((theta+d_small)*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin((thetamn+d_small)*pi/180)
                    *sin((theta+d_small)*pi/180)+(cos((theta+d_small)*pi/180)
                    *cos((thetamn+d_small)*pi/180)+1)*cos(((fi-fimn)+d_small)*pi/180))
        print('tgYmn_right=',tgYmn_right)
        tgYmn=(tgYmn_left+tgYmn_right)/2
        return tgYmn
tgY=func_tgYmn (th1,th3,fi1,fi2)
print('tgYmn=',tgY)
print()           
 
#Подставив полученные значения получаем векторные диаграммы направленности элемента
Etmn=-(a*f1+b*f2*exp1)
print('Etmn=',Etmn)
Etmn_abs=abs(Etmn)
print('Etmn_abs=',Etmn_abs)
print()
Efmn=b*f1-a*f2*exp1
print('Efmn=',Efmn)
Efmn_abs=abs(Efmn)
print('Efmn_abs=',Efmn_abs)
print()    


#Вводим длину волны и расчитываем далее волновое число и коэффициент "а"
lam=float(input ('Введите длину волны в метрах='))
print('Длина волны, м =',lam)
k=2*pi/lam
print('Волновое число=',k)

a=18.8/k
print('Коэффициент "а"=',a)
print()

rmnr=a*(sin(th1)*sin(th3)*cos(fi1-fi2)+cos(th3)*cos(th1))
print('rmnr=',rmnr)

psimn=k*a*(sin(th1)*sin(th3)*cos(fi1-fi2)+cos(th3)*cos(th1))
print('psimn=',psimn)

bmn=abs(a-1j*b)
print('bmn=',bmn)


















