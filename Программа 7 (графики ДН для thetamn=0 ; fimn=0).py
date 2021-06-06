import matplotlib.pyplot as plt
import numpy as np
from math import sin, cos, asin, acos, pi, atan
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
d_small=0.0000001

#Вводим длину волны и расчитываем далее волновое число, исходя из рабочей частоты:
#Скорость света в м/с:
c=299792458
#Рабочая частота в Гц:
f0=1.5445e9
#Рабочая длина волны в м:
lam=c/f0
print('Рабочая длина волны, м =',lam)
#Волновое число:
k=2*pi/lam
print('Волновое число=',k)
#Радиус сферы в м:
a=2
print('Радиус сферы, м=',a)
print()

#Амплитудное распределение.
#Рассматривается равномерное амплитудное распределение:
Amn=1

#p>=1
p=10
#f1(theta,fi)= f2(theta,fi)
#Вводим нормированную действительную амплитуду 1 (ДН элемента 1)
def func_f1(theta1,fi):
    if th2>=(-pi/2) and th2<=pi/2:
        #print('-pi/2<= Тета" <=pi/2')
        return cos(theta1)**p
    else:
        return 0

#Вводим функцию Amn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - Amn считается без приращений.
def func_Amn(theta,theta1,thetamn,fi,fimn):
    try:
        Amn=(cos(theta*pi/180)*sin(thetamn*pi/180)
             *cos((fi-fimn)*pi/180)-sin(theta*pi/180)
             *cos(thetamn*pi/180))/sin(theta1*pi/180)
        #print('Избежали попадания в неопределенность, Amn считается без приращений')
        return Amn         
    except ZeroDivisionError:
        print('Непределенность 0/0. Приращаем Тета, Тета",(Фи-Фиmn), Тетаmn')
        #print('Ищем значения слева и справа, после чего берем среднее между ними')
        Amn_left=(cos((theta-d_small)*pi/180)*sin((thetamn-d_small)*pi/180)
        *cos(((fi-fimn)-d_small)*pi/180)-sin((theta-d_small)*pi/180)
        *cos((thetamn-d_small)*pi/180))/sin((theta1-d_small)*pi/180)
        #print('Alfa_left=',Amn_left)
        Amn_right=(cos((theta+d_small)*pi/180)*sin((thetamn+d_small)*pi/180)
        *cos(((fi-fimn)+d_small)*pi/180)-sin((theta+d_small)*pi/180)
        *cos((thetamn+d_small)*pi/180))/sin((theta1+d_small)*pi/180)
        #print('Alfa_right=',Amn_right)
        Amn=(Amn_left+Amn_right)/2
        return Amn
            
#Вводим функцию Bmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - Bmn считается без приращений.    
def func_Bmn (theta,theta1,thetamn,fi,fimn):
    try:
        Bmn=sin(thetamn*pi/180)*sin((fi-fimn)*pi/180)/sin(theta1*pi/180)
        #print('Избежали попадания в неопределенность, Bmn считается без приращений')
        return Bmn      
    except ZeroDivisionError:        
        #print('Непределенность 0/0. Приращаем Тета, Тета",(Фи-Фиmn), Тетаmn')
        #print('Ищем значения слева и справа, после чего берем среднее между ними')
        Bmn_left=sin((thetamn-d_small)*pi/180)*sin(((fi-fimn)-d_small)
            *pi/180)/sin((theta1-d_small)*pi/180)
        #print('Beta_left1=',Bmn_left)
        Bmn_right=sin((thetamn+d_small)*pi/180)*sin(((fi-fimn)+d_small)
            *pi/180)/sin((theta1+d_small)*pi/180)
        #print('Beta_right1=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2
        return Bmn     

#Вводим функцию ctg проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - ctg считается без приращений.  
def func_ctgfi1 (theta,thetamn,fi,fimn):
    try:
        ctgfi1=(cos(thetamn*pi/180)*sin(theta*pi/180)
                *cos((fi-fimn)*pi/180)-sin(thetamn*pi/180)
                *cos(theta*pi/180))/(sin(theta*pi/180)*sin((fi-fimn)*pi/180)) 
        #print('Избежали попадания в неопределенность, ctg считается без приращений')
        return ctgfi1
    except ZeroDivisionError:
        #print('Непределенность 0/0. Приращаем Тета,(Фи-Фиmn), Тетаmn')
        #print('Ищем значения слева и справа, после чего берем среднее между ними')
        ctgfi1_left=(cos((thetamn-d_small)*pi/180)*sin((theta-d_small)*pi/180)
                     *cos(((fi-fimn)-d_small)*pi/180)-sin((thetamn-d_small)
                     *pi/180)*cos((theta-d_small)*pi/180))/(sin((theta-d_small)*pi/180)
                     *sin(((fi-fimn)-d_small)*pi/180))
        #print('Ctgfi"_left=',ctgfi1_left)
        ctgfi1_right=(cos((thetamn+d_small)*pi/180)*sin((theta+d_small)*pi/180)
                     *cos(((fi-fimn)+d_small)*pi/180)-sin((thetamn+d_small)
                     *pi/180)*cos((theta+d_small)*pi/180))/(sin((theta+d_small)*pi/180)
                     *sin(((fi-fimn)+d_small)*pi/180))
        #print('Ctgfi"_right6=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2        
        return ctgfi1                

#Вводим функцию tgYmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем все аргументы.
#Если избежали попадания в неопределенность - tgYmn считается без приращений.  
def func_tgYmn (theta,thetamn,fi,fimn):
    try:
        tgYmn=((cos(thetamn*pi/180)+cos(theta*pi/180))*sin((fi-fimn)
                *pi/180))/(sin(thetamn*pi/180)*sin(theta*pi/180)+(cos(theta
                *pi/180)*cos(thetamn*pi/180)+1)*cos((fi-fimn)*pi/180))
        #print('Избежали попадания в неопределенность, tgYmn  считается без приращений')
        return tgYmn
    except ZeroDivisionError:
        #print('Непределенность 0/0. Приращаем Тета,(Фи-Фиmn), Тетаmn')
        #print('Ищем значения слева и справа, после чего берем среднее между ними')
        tgYmn_left=((cos((thetamn-d_small)*pi/180)+cos((theta-d_small)*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin((thetamn-d_small)*pi/180)
                    *sin((theta-d_small)*pi/180)+(cos((theta-d_small)*pi/180)
                    *cos((thetamn-d_small)*pi/180)+1)*cos(((fi-fimn)-d_small)*pi/180))
        #print('tgYmn_left=',tgYmn_left)
        tgYmn_right=((cos((thetamn+d_small)*pi/180)+cos((theta+d_small)*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin((thetamn+d_small)*pi/180)
                    *sin((theta+d_small)*pi/180)+(cos((theta+d_small)*pi/180)
                    *cos((thetamn+d_small)*pi/180)+1)*cos(((fi-fimn)+d_small)*pi/180))
        #print('tgYmn_right=',tgYmn_right)
        tgYmn=(tgYmn_left+tgYmn_right)/2
        return tgYmn          

th3=0
fi2=0
print('Тетаmn"=',th3)
print('Фиmn"=',fi2)
#Создаем пустые списки, которые заполним в дальнейшем
#Пустой список Фи:
list_fi1=[]
#Пустые список Тета для сечений Фи=0,45,90:
list_th1_0=[]
list_th1_45=[]
list_th1_90=[]
#Пустые списки для сечения Фи=0 градусов:    
list_abs_Etheta1_0=[]    
list_abs_Efi1_0=[]
list_abs_Etheta12_0=[]
list_abs_Efi12_0=[]
#Пустые списки для сечения Фи=45 градусов:
list_abs_Etheta1_45=[]    
list_abs_Efi1_45=[]
list_abs_Etheta12_45=[]
list_abs_Efi12_45=[]
#Пустые списки для сечения Фи=90 градусов:
list_abs_Etheta1_90=[]    
list_abs_Efi1_90=[]
list_abs_Etheta12_90=[]
list_abs_Efi12_90=[]

#Считается последовательно 3 сечения Фи: 0,45,90 градусов
for fi1 in range(0,135,45):
    print()
    print('Фи=',fi1)
    list_fi1.append(fi1)
    print()           
    #Для каждого из сечений задается угол Тета от -90 до 90 градусов с шагом в 1 градус
    for th1 in range(-90,91,1):        
        #print('Тета=',th1)        
        th2=acos(sin(th3*pi/180)*sin(th1*pi/180)*
            cos((fi1-fi2)*pi/180)+cos(th3*pi/180)*
            cos(th1*pi/180))
        
        f1=func_f1(th2,fi1)
        #print('Нормированная ДН элемента f1(theta,fi)=',f1)

        f2=func_f1(th2,fi1)
        #print('Нормированная ДН элемента f2(theta,fi)=',f2)
        
        alfamn=func_Amn(th1,th2,th3,fi1,fi2)
        #print('Alfamn=',alfamn)

        betamn=func_Bmn(th1,th2,th3,fi1,fi2)
        #print('Betamn=',betamn)

        ctgfi=func_ctgfi1(th1,th3,fi1,fi2)
        #print('CtgFi"=',ctgfi)

        tgY=func_tgYmn (th1,th3,fi1,fi2)
        #print('tgYmn=',tgY)
        
        #Подставив полученные значения получаем векторные диаграммы направленности элемента
        Ethetamn=-(alfamn*f1+betamn*f2*exp1)
        #print('Ethetamn=',Ethetamn)
        Ethetamn_abs=abs(Ethetamn)
        #print('Ethetamn_abs=',Ethetamn_abs)

        Efimn=betamn*f1-alfamn*f2*exp1
        #print('Efimn=',Efimn)
        Efimn_abs=abs(Efimn)
        print('Efmn_abs=',Efimn_abs)
        print('betamn=',betamn)
        print('alfamn=',alfamn)        
        print()

        #Скалярное произведение радиус-вектора r, длиной 1, и вектора rmn:
        rmnr=a*(sin(th1*pi/180)*sin(th3*pi/180)*cos((fi1-fi2)*pi/180)+cos(th3*pi/180)
                *cos(th1*pi/180))
        #print('rmnr=',rmnr)

        #Фазовое распределение на элементах, которое нужно для формирования
        #остронаправленной ДН сферической АР:
        r0=(sin(0*pi/180)*sin(th3*pi/180)*cos((0-fi2)*pi/180)+cos(th3*pi/180)
            *cos(0*pi/180))        
        psimn=k*a*r0
        #print('psimn=',psimn)

        Bmn=abs(alfamn-1j*betamn)
        #print('Bmn=',Bmn)
        #print()

        #Из формулы (1)
        Etheta1=Amn*exp(k*rmnr-psimn)*Ethetamn
        #print('Etheta1=',Etheta1)
        abs_Etheta1=abs(Etheta1)
        #print('abs_Etheta1=',abs_Etheta1)
        
        Efi1=Amn*exp(k*rmnr-psimn)*Efimn
        #print('Efi1=',Efi1)
        abs_Efi1=abs(Efi1)
        #print('abs_Efi1=',abs_Efi1)    
        
        #Из выведенной формулы (12)        
        try:            
            Etheta12=-Amn*f1*Bmn*exp(k*rmnr-psimn)*exp(-1j
                        *atan(1/ctgfi))*exp(-1j*atan(betamn/alfamn))                        
        except ZeroDivisionError:            
            Etheta12=-Amn*f1*Bmn*exp(k*rmnr-psimn)*exp(-1j
                        *atan((1/(ctgfi-d_small)+1/(ctgfi+d_small))/2))*exp(-1j
                        *atan((betamn/(alfamn-d_small)+betamn/(alfamn+d_small))/2))

        #print('Etheta12=',Etheta12)
        abs_Etheta12=abs(Etheta12)
        #print('abs_Etheta12=',abs_Etheta12)        

        try:
            Efi12=1j*Amn*f1*Bmn*exp(k*rmnr-psimn)*exp(-1j
                    *atan(1/ctgfi))*exp(-1j*atan(betamn/alfamn))
        except ZeroDivisionError:
            Efi12=1j*Amn*f1*Bmn*exp(k*rmnr-psimn)*exp(-1j
                        *atan((1/(ctgfi-d_small)+1/(ctgfi+d_small))/2))*exp(-1j
                        *atan((betamn/(alfamn-d_small)+betamn/(alfamn+d_small))/2))
        #print('Efi12=',Efi12)
        abs_Efi12=abs(Efi12)
        #print('abs_Efi12=',abs_Efi12)

        #Полученные значения функций распределяются по спискам,
        #в зависимости от значений Фи
        #В каждом списке 181 элемент
        if fi1==0:
            list_th1_0.append(th1)
            list_abs_Etheta1_0.append(abs_Etheta1)
            list_abs_Efi1_0.append(abs_Efi1)
            list_abs_Etheta12_0.append(abs_Etheta12)
            list_abs_Efi12_0.append(abs_Efi12)
            
        elif fi1==45:
            list_th1_45.append(th1)
            list_abs_Etheta1_45.append(abs_Etheta1)
            list_abs_Efi1_45.append(abs_Efi1)
            list_abs_Etheta12_45.append(abs_Etheta12)
            list_abs_Efi12_45.append(abs_Efi12)            
        else:
            list_th1_90.append(th1)
            list_abs_Etheta1_90.append(abs_Etheta1)
            list_abs_Efi1_90.append(abs_Efi1)
            list_abs_Etheta12_90.append(abs_Etheta12)
            list_abs_Efi12_90.append(abs_Efi12)

        #Построение диаграмм направленности для каждой проекции каждого сечения:
        def graphics(*args):            
            plt.figure(n)            
            plt.plot(x,y)
            plt.xlabel('Тета [град]')
            plt.ylabel(y_title_Etheta1)
            plt.title(top_title_Etheta1)
            plt.grid()          

            plt.figure(n+1)            
            plt.plot(x,z)            
            plt.xlabel('Тета [град]')
            plt.ylabel(y_title_Etheta12)
            plt.title(top_title_Etheta12)
            plt.grid()

            plt.figure(n+2)            
            plt.plot(x,u)            
            plt.xlabel('Тета [град]')
            plt.ylabel(y_title_Efi1)            
            plt.title(top_title_Efi1)
            plt.grid()

            plt.figure(n+3)            
            plt.plot(x,s)            
            plt.xlabel('Тета [град]')
            plt.ylabel(y_title_Efi12)            
            plt.title(top_title_Efi12)
            plt.grid()
            
            return args
        #В зависимости от сечения строятся определенные графики
        #Сечение Фи = 0 градусов:
        if fi1==0:
            n=1
            x=list_th1_0
            y=list_abs_Etheta1_0
            z=list_abs_Etheta12_0
            u=list_abs_Efi1_0
            s=list_abs_Efi12_0

            y_title_Etheta1='abs_Etheta1_0'
            y_title_Etheta12='abs_Etheta12_0'
            y_title_Efi1='abs_Efi1_0'
            y_title_Efi12='abs_Efi12_0'
            
            top_title_Etheta1='Проекция по Тета формулы (1) сечения Фи=0 град'
            top_title_Etheta12='Проекция по Тета формулы (12) сечения Фи=0 град'
            top_title_Efi1='Проекция по Фи формулы (1) сечения Фи=0 град'
            top_title_Efi12='Проекция по Фи формулы (12) сечения Фи=0 град'
        #Сечение Фи = 45 градусов:
        elif fi1==45:
            n=5
            x=list_th1_45
            y=list_abs_Etheta1_45
            z=list_abs_Etheta12_45
            u=list_abs_Efi1_45
            s=list_abs_Efi12_45

            y_title_Etheta1='abs_Etheta1_45'
            y_title_Etheta12='abs_Etheta12_45'
            y_title_Efi1='abs_Efi1_45'
            y_title_Efi12='abs_Efi12_45'
            
            top_title_Etheta1='Проекция по Тета формулы (1) сечения Фи=45 град'
            top_title_Etheta12='Проекция по Тета формулы (12) сечения Фи=45 град'
            top_title_Efi1='Проекция по Фи формулы (1) сечения Фи=45 град'
            top_title_Efi12='Проекция по Фи формулы (12) сечения Фи=45 град'
        #Сечение Фи = 90 градусов:
        else:
            n=9
            x=list_th1_90
            y=list_abs_Etheta1_90
            z=list_abs_Etheta12_90
            u=list_abs_Efi1_90
            s=list_abs_Efi12_90

            y_title_Etheta1='abs_Etheta1_90'
            y_title_Etheta12='abs_Etheta12_90'
            y_title_Efi1='abs_Efi1_90'
            y_title_Efi12='abs_Efi12_90'
            
            top_title_Etheta1='Проекция по Тета формулы (1) сечения Фи=90 град'
            top_title_Etheta12='Проекция по Тета формулы (12) сечения Фи=90 град'
            top_title_Efi1='Проекция по Фи формулы (1) сечения Фи=90 град'
            top_title_Efi12='Проекция по Фи формулы (12) сечения Фи=90 град'            

        graphics(n,x,y,z,u,s,y_title_Etheta1,top_title_Etheta1,
                 y_title_Etheta12,top_title_Etheta12,top_title_Efi1,
                 top_title_Efi12)

plt.show()
#После рассчета всех функций для данного сечения начинается расчет следующего сечения


    





    



