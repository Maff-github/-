import numpy as np
import sympy as sp
from scipy.optimize import fsolve
from numpy import sin, cos, arcsin, arccos, pi, exp

#Вводим экспоненту, отражающую изменение фазы элемента антенны с круговой поляризацией
def e (beta):
    return exp(-1j*beta)
#по условию beta=pi/2
beta1=pi/2
exp1=e(beta1)
print('exp=',exp1)


#k- любое натуральное число, которое используется в неопределенностях
k=1
d_small=0.0001

th1=float(input ('Введите Тета(от 0 до 180 градусов)='))
th3=float(input ('Введите Тетаmn='))
fi1=float(input ('Введите Фи(от 0 до 360 градусов)='))
fi2=float(input ('Введите Фиmn='))
print()
#Для проверки работы в синугулярностях убрать # для нужного значения Тета'
####Сингулярности для Альфа:####
###Знаменатель равен 0 при:###
#th2=0
#th2=pi*k
###Числитель равен 0 при:###
##Первый случай, когда th1 и th3 одновременно равны:##
# 1.1) 
#th1=pi*k/2
#th1=3*pi*k/2
#th3=pi*k/2
#th3=3*pi*k/2
# 1.2)
#th3=0
#th3=pi*k
#th1=0
#th1=pi*k
##Второй случай, когда th1 и fi1 одновременно равны:##
#fi1=pi*k/2+fi2
#fi1=3*pi*k/2+fi2
#th1=0
#th1=pi*k
##Третий случай, когда th3 и fi1 одновременно равны:##
#fi1=pi*k/2+fi2
#fi1=3*pi*k/2+fi2
#th3=pi*k/2
#th3=3*pi*k/2

####Сингулярности для Бета:####
#th2=0
#th2=pi*k
#th3=0
#th3=pi*k
#fi1=fi2
#fi1=pi*k+fi2

###Сингулярности для котангенса:##
#th1=0
#th1=pi*k
#th1=pi*k/2
#th1=3*pi*k/2
#th3=0
#th3=pi*k
#th3=pi*k/2
#th3=3*pi*k/2
#fi1=fi2
#fi1=pi*k+fi2
#fi1=pi*k/2+fi2
#fi1=3*pi*k/2+fi2
th2=arccos(sin(th3*pi/180)*sin(th1*pi/180)*
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
#Если попали в неопределенность - приращаем те аргументы, из-за которых она возникает.
#Если избежали попадания в неопределенность - Amn считается без приращений.
def func_Amn(theta,theta1,thetamn,fi,fimn):

#sin(Тета)=0. Деление на ноль.
#Берем среднее значения между правым и левым приращением Тета":
    if th2==0 or th2==pi*k:
        print('Обнуление sin(Тета"). Деление на 0. Приращаем Тета"')
        Amn_left=(cos(theta*pi/180)*sin(thetamn*pi/180)
        *cos((fi-fimn)*pi/180)-sin(theta*pi/180)
        *cos(thetamn*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left1=',Amn_left)
        Amn_right=(cos(theta*pi/180)*sin(thetamn*pi/180)
        *cos((fi-fimn)*pi/180)-sin(theta*pi/180)
        *cos(thetamn*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right1=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn

#sin(Тета)=0, sin(Тетаmn)=0, sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета, Тета" и Тетаmn:  
    elif ((th1==pi*k/2 or th1==3*pi*k/2)\
          and (th3==pi*k/2 or th3==3*pi*k/2)) \
          or ((th3==0 or th3==pi*k) and (th1==0 or th1==pi*k))\
          and ((fi1==pi/2+fi2 or fi1==3*pi/2+fi2)\
               and (th1==0 or th1==pi*k))\
          and ((fi1==pi/2+fi2 or fi1==3*pi/2+fi2)\
               and (th3==pi*k/2 or th1==3*pi*k/2)):      
        print('Непределенность 0/0. Приращаем Тета, Тета" и Тетаmn')
        Amn_left=(cos((theta-d_small)*pi/180)*sin((thetamn-d_small)
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)-sin((theta-d_small)
        *pi/180)*cos((thetamn-d_small)*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left2=',Amn_left)        
        Amn_right=(cos((theta+d_small)*pi/180)*sin((thetamn+d_small)
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)-sin((theta+d_small)
        *np.pi/180)*cos((thetamn+d_small)*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right2=',Amn_right)
        Amn=(Amn_left+Amn_right)/2        
        return Amn        

#sin(Тета)=0, cos(Тетаmn)=0, sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета, Тета" и Тетаmn:  
    elif ((th1==pi*k/2 or th1==3*pi*k/2)\
          and (th3==pi*k/2 or th3==3*pi*k/2)) \
          or ((th3==0 or th3==pi*k) and (th1==0 or th1==pi*k)):    
        print('Непределенность 0/0. Приращаем Тета, Тета" и Тетаmn')
        Amn_left=(cos((theta-d_small)*pi/180)*sin((thetamn-d_small)
        *pi/180)*cos((fi-fimn)*pi/180)-sin(theta*pi/180)
        *cos((thetamn-d_small)*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left2=',Amn_left)
        Amn_right=(cos((theta+d_small)*pi/180)*sin((thetamn+d_small)
        *pi/180)*cos((fi-fimn)*pi/180)-sin(theta*pi/180)
        *cos((thetamn+d_small)*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right2=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn

#cos(Фи-Фиmn)=0, sin(Тета)=0, sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета,(Фи-Фиmn) и Тета":
    elif (fi1==pi/2+fi2 or fi1==3*pi/2+fi2) and (th1==0 or th1==pi*k):         
        print('Непределенность 0/0. Приращаем Тета, Тета" и (Фи-Фиmn)')
        Amn_left=(cos((theta-d_small)*pi/180)*sin(thetamn*pi/180)
        *cos(((fi-fimn)-d_small)*pi/180)-sin((theta-d_small)*pi/180)
        *cos(thetamn*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left3=',Amn_left)
        Amn_right=(cos((theta+d_small)*pi/180)*sin(thetamn*pi/180)
        *cos(((fi-fimn)+d_small)*pi/180)-sin((theta+d_small)*pi/180)
        *cos(thetamn*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right3=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn

#cos(Фи-Фиmn)=0, cos(Тетаmn), sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тетаmn,(Фи-Фиmn) и Тета":  
    elif (fi1==pi/2+fi2 or fi1==3*pi/2+fi2)\
         and (th3==pi*k/2 or th1==3*pi*k/2):       
        print('Непределенность 0/0. Приращаем Тетаmn, Тета" и (Фи-Фиmn)')
        Amn_left=(cos(theta*np.pi/180)*sin((thetamn-d_small)*pi/180)
        *cos(((fi-fimn)-d_small)*pi/180)-sin(theta*pi/180)
        *cos((thetamn-d_small)*pi/180))/sin((theta1-d_small)*pi/180)
        print('Alfa_left4=',Amn_left)
        Amn_right=(cos(theta*pi/180)*sin((thetamn+d_small)*pi/180)
        *cos(((fi-fimn)+d_small)*pi/180)-sin(theta*pi/180)
        *cos((thetamn+d_small)*pi/180))/sin((theta1+d_small)*pi/180)
        print('Alfa_right4=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn                    

#Так как избежали попадания в неопределенность, то считаем Amn без приращений 
    else:      
        print('Избежали попадания в неопределенность, Amn считается без приращений')
        Amn=(cos(theta*pi/180)*sin(thetamn*pi/180)
         *cos((fi-fimn)*pi/180)-sin(theta*pi/180)
         *cos(thetamn*pi/180))/sin(theta1*pi/180)        
    return Amn        
a=func_Amn(th1,th2,th3,fi1,fi2)
print('Alfa=',a)
print()

#Вводим функцию Bmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем те аргументы, из-за которых она возникает.
#Если избежали попадания в неопределенность - Bmn считается без приращений.    
def func_Bmn (theta,theta1,thetamn,fi,fimn):    

#sin(Тета")=0. Деление на ноль.    
#Берем среднее значения между правым и левым приращением Тета":  
    if th2==0 or th2==np.pi*k:      
        print('Обнуление sin(Тета"). Деление на 0. Приращаем Тета"')
        Bmn_left=sin(thetamn*pi/180)*sin((fi-fimn)
        *pi/180)/sin((theta1-d_small)*pi/180)
        print('Beta_left1=',Bmn_left)
        Bmn_right=sin(thetamn*pi/180)*sin((fi-fimn)
        *pi/180)/sin((theta1+d_small)*pi/180)
        print('Beta_right1=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2
        return Bmn

#sin(Фи-Фиmn)=0, sin(Тетаmn)=0, sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями (Фи-Фиmn), Тета" и Тетаmn:
    elif (th3==0 and fi1==fi2) or (th3==0 and fi1==pi*k+fi2) \
         or (th3==pi*k and fi1== fi2) or (th3==pi*k and fi1==pi*k+fi2):         
        print('Непределенность 0/0. Приращаем Тетаmn, Тета" и (Фи-Фиmn)')
        Bmn_left=sin(thetamn-d_small*pi/180)*sin(((fi-d_small)-fimn)
        *pi/180)/sin((theta1-d_small)*pi/180)
        print('Beta_left2=',Bmn_left)
        Bmn_right=sin(thetamn+d_small*pi/180)*sin(((fi+d_small)-fimn)
        *pi/180)/sin((theta1+d_small)*pi/180)
        print('Beta_right2=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2                        
        return Bmn
    
#sin(Тетаmn)=0 и sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета" и Тетаmn:     
    elif th3==0 or th3==pi*k:        
        print('Непределенность 0/0. Приращаем Тетаmn и Тета"')
        Bmn_left=sin((thetamn-d_small)*pi/180)*sin((fi-fimn)
        *pi/180)/sin((theta1-d_small)*pi/180)
        print('Beta_left3=',Bmn_left)
        Bmn_right=sin((thetamn+d_small)*pi/180)*sin((fi-fimn)
        *pi/180)/sin((theta1+d_small)*pi/180)
        print('Beta_right3=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2                        
        return Bmn
    
#sin(Фи-Фиmn)=0 и sin(Тета")=0. Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями (Фи-Фиmn) и Тета":      
    elif fi1==fi2 or fi1==pi*k+fi2:       
        print('Непределенность 0/0. Приращаем (Фи-Фиmn) и Тета"')
        Bmn_left=sin(thetamn*pi/180)*sin(((fi-fimn)-d_small)
        *pi/180)/sin((theta1-d_small)*pi/180)
        print('Beta_left4=',Bmn_left)
        Bmn_right=np.sin(thetamn*np.pi/180)*np.sin(((fi-fimn)+d_small)
        *pi/180)/sin((theta1+d_small)*pi/180)
        print('Beta_right4=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2        
        return Bmn    

#Так как избежали попадания в неопределенность, то считаем Amn без приращений 
    else:         
        print('Избежали попадания в неопределенность, Bmn считается без приращений')
        Bmn=sin(thetamn*pi/180)*sin((fi-fimn)
            *pi/180)/sin(theta1*pi/180)        
        return Bmn
b=func_Bmn(th1,th2,th3,fi1,fi2)
print('Beta=',b)
print()


#Вводим функцию ctg проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем те аргументы, из-за которых она возникает.
#Если избежали попадания в неопределенность - ctg считается без приращений.  
def func_ctgfi1 (theta,thetamn,fi,fimn):    

#sin(Тета)=0 и sin(Фи-Фиmn)=0. Деление на 0.    
#Берем средние значения между правыми и левыми приращениями (Фи-Фиmn) и Тета:
    if (th1==0 and fi1==fi2) or (th1==0 and fi1==pi*k+fi2)\
       or (th1==pi*k and fi1==fi2) or (th1==pi*k and fi1==pi*k+fi2): 
        print('Обнуление sin(Тета) и sin(fi-fimn). Деление на 0. Приращаем Тета и (Фи-Фиmn)')
        ctgfi1_left=(cos(thetamn*pi/180)*sin((theta-d_small)
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)-sin(thetamn
        *pi/180)*cos((theta-d_small)*pi/180))/(sin((theta-d_small)
        *pi/180)*sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left1=',ctgfi1_left)
        ctgfi1_right=(cos(thetamn*pi/180)*sin((theta+d_small)
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)-sin(thetamn
        *pi/180)*cos((theta+d_small)*pi/180))/(sin((theta+d_small)
        *np.pi/180)*sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right1=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1

#sin(Тета)=0. Деление на 0.    
#Берем среднее значение между правым и левым приращением Тета:
    elif th1==0 or th1==pi*k:
        print('Обнуление sin(Тета). Деление на 0. Приращаем Тета')
        ctgfi1_left=(cos(thetamn*pi/180)*sin((theta-d_small)
        *pi/180)*cos((fi-fimn)*pi/180)-sin(thetamn
        *pi/180)*cos((theta-d_small)*pi/180))/(sin((theta-d_small)
        *pi/180)*sin((fi-fimn)*pi/180))
        print('Ctgfi"_left2=',ctgfi1_left)
        ctgfi1_right=(cos(thetamn*pi/180)*sin((theta+d_small)
        *pi/180)*cos((fi-fimn)*pi/180)-sin(thetamn
        *pi/180)*cos((theta+d_small)*pi/180))/(sin((theta+d_small)
        *pi/180)*sin((fi-fimn)*pi/180))
        print('Ctgfi"_right2=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1

#sin(Фи-Фиmn)=0. Деление на 0.    
#Берем среднее значение между правым и левым приращением Тета:    
    elif fi1==fi2 or fi1==pi*k+fi2:
        print('Обнуление sin(Фи-Фиmn). Деление на 0. Приращаем (Фи-Фиmn)')
        ctgfi1_left=(cos(thetamn*pi/180)*sin(theta
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)-sin(thetamn
        *pi/180)*cos(theta*pi/180))/(sin(theta
        *pi/180)*sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left3=',ctgfi1_left)
        ctgfi1_right=(cos(thetamn*pi/180)*sin(theta
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)-sin(thetamn
        *pi/180)*cos(theta*pi/180))/(sin(theta
        *pi/180)*sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right3=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1

#sin(Тета)=0, sin(Фи-Фиmn)=0, cos(Фи-Фиmn)=0, cos(Тетаmn)=0, cos(Тета)=0, sin(Тетаmn)=0.
#Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета, (Фи-Фиmn) и Тетаmn:     
    elif (th1==pi*k/2 and fi1==pi*k/2+fi2)\
       or (th1==pi*k/2 and fi1==3*pi*k+fi2)\
       or(th1==3*pi*k/2 and fi1==pi*k/2+fi2)\
       or(th1==3*pi*k/2 and fi1==3*pi*k/2+fi2)\
       or(th3==0 and fi1==pi*k/2+fi2)\
       or(th3==pi*k and fi1==pi*k/2+fi2)\
       or(th3==0 and fi1==3*pi*k/2+fi2)\
       or(th3==pi*k and fi1==3*pi*k/2+fi2):
        print('Непределенность 0/0. Приращаем Тетаmn, Тета, (Фи-Фиmn)')
        ctgfi1_left=(cos((thetamn-d_small)*pi/180)*sin((theta-d_small)
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)
        -sin((thetamn-d_small)*pi/180)*cos((theta-d_small)
        *pi/180))/(sin((theta-d_small)*pi/180)
        *sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left4=',ctgfi1_left)
        ctgfi1_right=(cos((thetamn+d_small)*pi/180)*sin((theta+d_small)
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)
        -sin((thetamn+d_small)*pi/180)*cos((theta+d_small)
        *pi/180))/(sin((theta+d_small)*pi/180)*sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right4=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1                            

#sin(Тета)=0, sin(Фи-Фиmn)=0, cos(Фи-Фиmn)=0, cos(Тета)=0.
#Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями Тета и (Фи-Фиmn): 
    elif (th1==pi*k/2 and fi1==pi*k/2+fi2)\
         or (th1==pi*k/2 and fi1==3*pi*k+fi2)\
         or(th1==3*pi*k/2 and fi1==pi*k/2+fi2)\
         or(th1==3*pi*k/2 and fi1==3*pi*k/2+fi2):
        print('Непределенность 0/0. Приращаем Тета и(Фи-Фиmn)')
        ctgfi1_left=(cos(thetamn*pi/180)*sin((theta-d_small)
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)
        -sin(thetamn*pi/180)*cos((theta-d_small)
        *pi/180))/(sin((theta-d_small)*pi/180)
        *sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left5=',ctgfi1_left)
        ctgfi1_right=(cos(thetamn*pi/180)*sin((theta+d_small)
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)
        -sin(thetamn*pi/180)*cos((theta+d_small)
        *pi/180))/(sin((theta+d_small)*pi/180)*sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right5=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1

#sin(Тета)=0, sin(Фи-Фиmn)=0, cos(Фи-Фиmn)=0, sin(Тетаmn)=0.
#Неопределенность типа 0/0.    
#Берем средние значения между правыми и левыми приращениями (Фи-Фиmn) и Тетаmn: 
    elif (th3==0 and fi1==pi*k/2+fi2)\
       or(th3==pi*k and fi1==pi*k/2+fi2)\
       or(th3==0 and fi1==3*pi*k/2+fi2)\
       or(th3==pi*k and fi1==3*pi*k/2+fi2):
        print('Непределенность 0/0. Приращаем Тетаmn и(Фи-Фиmn)')
        ctgfi1_left=(cos((thetamn-d_small)*pi/180)*sin(theta
        *pi/180)*cos(((fi-fimn)-d_small)*pi/180)
        -sin((thetamn-d_small)*pi/180)*cos(theta
        *pi/180))/(sin(theta*np.pi/180)
        *sin(((fi-fimn)-d_small)*pi/180))
        print('Ctgfi"_left6=',ctgfi1_left)
        ctgfi1_right=(cos((thetamn+d_small)*pi/180)*sin(theta
        *pi/180)*cos(((fi-fimn)+d_small)*pi/180)
        -sin((thetamn+d_small)*pi/180)*cos(theta
        *pi/180))/(sin(theta*pi/180)*sin(((fi-fimn)+d_small)*pi/180))
        print('Ctgfi"_right6=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1                                                              

#Так как избежали попадания в неопределенность, то считаем ctg без приращений 
    else:
        print('Избежали попадания в неопределенность, ctg считается без приращений')
        ctgfi1=(cos(thetamn*pi/180)*sin(theta*pi/180)
                *cos((fi-fimn)*pi/180)-sin(thetamn*pi/180)
                *cos(theta*pi/180))/(sin(theta*pi/180)*sin((fi-fimn)*pi/180))
        return ctgfi1
ctg=func_ctgfi1(th1,th3,fi1,fi2)
print('CtgFi"=',ctg)
print()









#Вводим функцию tgYmn проверяя ее на наличие попадания в неопределенности.
#Если попали в неопределенность - приращаем те аргументы, из-за которых она возникает.
#Если избежали попадания в неопределенность - tgYmn считается без приращений.  
def func_tgYmn (theta,thetamn,fi,fimn):

#sin(Тетаmn)=0, cos(Фи-Фиmn)=0, sin(Тета)=0, cos(Тета)*cos(Тетаmn)+1=0,
#(cos(Тетаmn)+cos(Тета))=0. Неопределенность 0/0.
#Берем средние значения между правыми и левыми приращениями Тетаmn, (Фи-Фиmn) и Тета:    
    if (th3==0 and fi1==pi*k/2+fi2 and th1==pi*k)\
       or (th3==0 and fi1==3*pi*k/2+fi2 and th1==pi*k)\
       or (th3==pi*k and fi1==pi*k/2+fi2 and th1==0) \
       or (th3==pi*k and fi1==3*pi*k/2+fi2 and th1==0):
        print('Неопределенность 0/0. Приращаем Тета, Тетаmn, (Фи-Фиmn)')
        tgYmn_left=((cos((thetamn-d_small)*pi/180)+cos((theta-d_small)*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin((thetamn-d_small)*pi/180)
                    *sin((theta-d_small)*pi/180)+(cos((theta-d_small)*pi/180)
                    *cos((thetamn-d_small)*pi/180)+1)*cos(((fi-fimn)-d_small)*pi/180))
        print('tgYmn_left1=',tgYmn_left)
        tgYmn_right=((cos((thetamn+d_small)*pi/180)+cos((theta+d_small)*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin((thetamn+d_small)*pi/180)
                    *sin((theta+d_small)*pi/180)+(cos((theta+d_small)*pi/180)
                    *cos((thetamn+d_small)*pi/180)+1)*cos(((fi-fimn)+d_small)*pi/180))
        print('tgYmn_right1=',tgYmn_right)
        thYmn=(tgYmn_left+tgYmn_right)/2
        return tgYmn

#sin(Тетаmn)=0, cos(Тета)*cos(Тетаmn)+1=0. Неопределенность 0/0.
#(cos(Тетаmn)+cos(Тета))=0
#Берем средние значения между правыми и левыми приращениями Тетаmn и Тета:     
    elif (th3==0 and th1==pi*k) or (th3==pi*k and th1==0):
        print('Неопределенность 0/0. Приращаем Тета и Тетаmn')
        tgYmn_left=((cos((thetamn-d_small)*pi/180)+cos((theta-d_small)*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin((thetamn-d_small)*pi/180)
                    *sin((theta-d_small)*pi/180)+(cos((theta-d_small)*pi/180)
                    *cos((thetamn-d_small)*pi/180)+1)*cos((fi-fimn)*pi/180))
        print('tgYmn_left2=',tgYmn_left)
        tgYmn_right=((cos((thetamn+d_small)*pi/180)+cos((theta+d_small)*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin((thetamn+d_small)*pi/180)
                    *sin((theta+d_small)*pi/180)+(cos((theta+d_small)*pi/180)
                    *cos((thetamn+d_small)*pi/180)+1)*cos((fi-fimn)*pi/180))
        print('tgYmn_right2=',tgYmn_right)
        tgYmn=(tgYmn_left+tgYmn_right)/2        
        return tgYmn

#sin(Тетаmn)=0, cos(Фи-Фиmn)=0. Деление на 0.
#Берем средние значения между правыми и левыми приращениями Тетаmn и (Фи-Фиmn):   
    elif (th3==0 and fi1==pi*k/2+fi2) or (th3==0 and fi1==3*pi*k/2+fi2)\
         or (th3==pi*k and fi1==pi*k/2+fi2) or (th3==pi*k and fi1==3*pi*k/2+fi2):
        print('Деление на 0. Приращаем (Фи-Фиmn) и Тетаmn')
        tgYmn_left=((cos((thetamn-d_small)*pi/180)+cos(theta*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin((thetamn-d_small)*pi/180)
                    *sin(theta*pi/180)+(cos(theta*pi/180)
                    *cos((thetamn-d_small)*pi/180)+1)*cos(((fi-fimn)-d_small)*pi/180))
        print('tgYmn_left3=',tgYmn_left)
        tgYmn_right=((cos((thetamn+d_small)*pi/180)+cos(theta*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin((thetamn+d_small)*pi/180)
                    *sin(theta*pi/180)+(cos(theta*pi/180)
                    *cos((thetamn+d_small)*pi/180)+1)*cos(((fi-fimn)+d_small)*pi/180))
        print('tgYmn_right3=',tgYmn_right)
        tgYmn=(tgYmn_left+tgYmn_right)/2
        return tgYmn

#sin(Тета)=0, cos(Фи-Фиmn)=0. Деление на 0.
#Берем средние значения между правыми и левыми приращениями Тета и (Фи-Фиmn):   
    elif (th1==0 and fi1==pi*k/2+fi2) or (th1==0 and fi1==3*pi*k/2+fi2)\
         or (th1==pi*k and fi1==pi*k/2+fi2) or (th1==pi*k and fi1==3*pi*k/2+fi2):
        print('Деление на 0. Приращаем (Фи-Фиmn) и Тета')
        tgYmn_left=((cos(thetamn*pi/180)+cos((theta-d_small)*pi/180))
                    *sin(((fi-fimn)-d_small)*pi/180))/(sin(thetamn*pi/180)
                    *sin((theta-d_small)*pi/180)+(cos((theta-d_small)*pi/180)
                    *cos(thetamn*pi/180)+1)*cos(((fi-fimn)-d_small)*pi/180))
        print('tgYmn_left4=',tgYmn_left4)
        tgYmn_right=((cos(thetamn*pi/180)+cos((theta+d_small)*pi/180))
                    *sin(((fi-fimn)+d_small)*pi/180))/(sin(thetamn*pi/180)
                    *sin((theta+d_small)*pi/180)+(cos((theta+d_small)*pi/180)
                    *cos(thetamn*pi/180)+1)*cos(((fi-fimn)+d_small)*pi/180))
        print('thYmn_right4=',thYmn_right)
        tgYmn=(tgYmn_left+tgYmn_right)/2
        return tgYmn
    
#Так как избежали попадания в неопределенность, то считаем ctg без приращений 
    else:
        print('Избежали попадания в неопределенность, tgYmn считается без приращений')
        tgYmn=((cos(thetamn*pi/180)+cos(theta*pi/180))
                    *sin((fi-fimn)*pi/180))/(sin(thetamn*pi/180)
                    *sin(theta*pi/180)+(cos(theta*pi/180)
                    *cos(thetamn*pi/180)+1)*cos((fi-fimn)*pi/180))
        return tgYmn
    

print('tgYmn=',func_tgYmn (th1,th3,fi1,fi2))
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



                
                

    

