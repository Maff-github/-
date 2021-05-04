import numpy as np
import sympy as sp
from scipy.optimize import fsolve

def f1(x1):    
    return x1
f1=f1(float(input('Введите нормированную амплитуду 1(ДН элемента)=')))
print('f1=',f1)

def f2(x2):            
    return x2
f2=f2(float(input('Введите нормированную амплитуду 2(ДН элемента)=')))
print('f2=',f2)    

def e (beta):
    return np.exp(-1j*beta)
beta1=float(input('Введите коэффициент, определяющий поляризацию поля='))
exp=e(beta1)
print('exp=',exp)

#k- любое натуральное число
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
#th2=np.pi*k
###Числитель равен 0 при:###
##Первый случай, когда th1 и th3 одновременно равны:##
# 1.1) 
#th1=np.pi*k/2
#th1=3*np.pi*k/2
#th3=np.pi*k/2
#th3=3*np.pi*k/2
# 1.2)
#th3=0
#th3=np.pi*k
#th1=0
#th1=np.pi*k
##Второй случай, когда th1 и fi1 одновременно равны:##
#fi1=np.pi*k/2+fi2
#fi1=3*np.pi*k/2+fi2
#th1=0
#th1=np.pi*k
##Третий случай, когда th3 и fi1 одновременно равны:##
#fi1=np.pi*k/2+fi2
#fi1=3*np.pi*k/2+fi2
#th3=np.pi*k/2
#th3=3*np.pi*k/2

####Сингулярности для Бета:####
#th2=0
#th2=np.pi*k
#th3=0
#th3=np.pi*k
#fi1=fi2
#fi1=np.pi*k+fi2

###Сингулярности для котангенса:##
#th1=0
#th1=np.pi*k
#th1=np.pi*k/2
#th1=3*np.pi*k/2
#th3=0
#th3=np.pi*k
#th3=np.pi*k/2
#th3=3*np.pi*k/2
#fi1=fi2
#fi1=np.pi*k+fi2
#fi1=np.pi*k/2+fi2
#fi1=3*np.pi*k/2+fi2
th2=np.arccos(np.sin(th3*np.pi/180)*np.sin(th1*np.pi/180)*
              np.cos((fi1-fi2)*np.pi/180)+np.cos(th3*np.pi/180)*
              np.cos(th1*np.pi/180))
print('Тета"=',th2)
print()

def func_Amn(theta,theta1,thetamn,fi,fimn):        
    if th2==0 or th2==np.pi*k:
        print('Попали в сингулярность, ищем значения слева и справа')
        Amn_left=(np.cos(theta*np.pi/180)*np.sin(thetamn*np.pi/180)
        *np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos(thetamn*np.pi/180))/np.sin((theta1-d_small)*np.pi/180)
        print('Alfa_left1=',Amn_left)
        Amn_right=(np.cos(theta*np.pi/180)*np.sin(thetamn*np.pi/180)
        *np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos(thetamn*np.pi/180))/np.sin((theta1+d_small)*np.pi/180)
        print('Alfa_right1=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn
    elif ((th1==np.pi*k/2 or th1==3*np.pi*k/2)\
          and (th3==np.pi*k/2 or th3==3*np.pi*k/2)) \
          or ((th3==0 or th3==np.pi*k) and (th1==0 or th1==np.pi*k))\
          and ((fi1==np.pi/2+fi2 or fi1==3*np.pi/2+fi2)\
               and (th1==0 or th1==np.pi*k))\
          and ((fi1==np.pi/2+fi2 or fi1==3*np.pi/2+fi2)\
               and (th3==np.pi*k/2 or th1==3*np.pi*k/2)):
        print('Попали в сингулярность, ищем значения слева и справа')
        Amn_left=(np.cos((theta-d_small)*np.pi/180)*np.sin((thetamn-d_small)
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)-np.sin((theta-d_small)
        *np.pi/180)*np.cos((thetamn-d_small)*np.pi/180))/np.sin((theta1-d_small)
                                                      *np.pi/180)
        print('Alfa_left2=',Amn_left)        
        Amn_right=(np.cos((theta+d_small)*np.pi/180)*np.sin((thetamn+d_small)
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)-np.sin((theta+d_small)
        *np.pi/180)*np.cos((thetamn+d_small)*np.pi/180))/np.sin((theta1+d_small)
                                                      *np.pi/180)
        print('Alfa_right2=',Amn_right)
        Amn=(Amn_left+Amn_right)/2        
        return Amn        
    elif ((th1==np.pi*k/2 or th1==3*np.pi*k/2)\
          and (th3==np.pi*k/2 or th3==3*np.pi*k/2)) \
          or ((th3==0 or th3==np.pi*k) and (th1==0 or th1==np.pi*k)):
        print('Попали в сингулярность, ищем значения слева и справа')
        Amn_left=(np.cos((theta-d_small)*np.pi/180)*np.sin((thetamn-d_small)
        *np.pi/180)*np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos((thetamn-d_small)*np.pi/180))/np.sin((theta1-d_small)
                                                      *np.pi/180)
        print('Alfa_left2=',Amn_left)
        Amn_right=(np.cos((theta+d_small)*np.pi/180)*np.sin((thetamn+d_small)
        *np.pi/180)*np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos((thetamn+d_small)*np.pi/180))/np.sin((theta1+d_small)
                                                      *np.pi/180)
        print('Alfa_right2=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn
    elif (fi1==np.pi/2+fi2 or fi1==3*np.pi/2+fi2) and (th1==0 or th1==np.pi*k):
        print('Попали в сингулярность, ищем значения слева и справа')
        Amn_left=(np.cos((theta-d_small)*np.pi/180)*np.sin(thetamn*np.pi/180)
        *np.cos(((fi-fimn)-d_small)*np.pi/180)-np.sin((theta-d_small)*np.pi/180)
        *np.cos(thetamn*np.pi/180))/np.sin((theta1-d_small)*np.pi/180)
        print('Alfa_left3=',Amn_left)
        Amn_right=(np.cos((theta+d_small)*np.pi/180)*np.sin(thetamn*np.pi/180)
        *np.cos(((fi-fimn)+d_small)*np.pi/180)-np.sin((theta+d_small)*np.pi/180)
        *np.cos(thetamn*np.pi/180))/np.sin((theta1+d_small)*np.pi/180)
        print('Alfa_right3=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn
    elif (fi1==np.pi/2+fi2 or fi1==3*np.pi/2+fi2)\
         and (th3==np.pi*k/2 or th1==3*np.pi*k/2):
        print('Попали в сингулярность, ищем значения слева и справа')
        Amn_left=(np.cos(theta*np.pi/180)*np.sin((thetamn-d_small)*np.pi/180)
        *np.cos(((fi-fimn)-d_small)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos((thetamn-d_small)*np.pi/180))/np.sin((theta1-d_small)
                                                               *np.pi/180)
        print('Alfa_left4=',Amn_left)
        Amn_right=(np.cos(theta*np.pi/180)*np.sin((thetamn+d_small)*np.pi/180)
        *np.cos(((fi-fimn)+d_small)*np.pi/180)-np.sin(theta*np.pi/180)
        *np.cos((thetamn+d_small)*np.pi/180))/np.sin((theta1+d_small)
                                                                *np.pi/180)
        print('Alfa_right4=',Amn_right)
        Amn=(Amn_left+Amn_right)/2                        
        return Amn                    
    else:
        print('Избежали попадания в сингулярность')
        Amn=(np.cos(theta*np.pi/180)*np.sin(thetamn*np.pi/180)
         *np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
         *np.cos(thetamn*np.pi/180))/np.sin(theta1*np.pi/180)        
    return Amn        
a=func_Amn(th1,th2,th3,fi1,fi2)
print('Alfa=',a)
print()
    
def func_Bmn (theta,theta1,thetamn,fi,fimn):
    if th2==0 or th2==np.pi*k:
        print('Попали в сингулярность, ищем значения слева и справа')
        Bmn_left=np.sin(thetamn*np.pi/180)*np.sin((fi-fimn)
        *np.pi/180)/np.sin((theta1-d_small)*np.pi/180)
        print('Beta_left1=',Bmn_left)
        Bmn_right=np.sin(thetamn*np.pi/180)*np.sin((fi-fimn)
        *np.pi/180)/np.sin((theta1+d_small)*np.pi/180)
        print('Beta_right1=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2
        return Bmn
    elif (th3==0 and fi1==fi2) or (th3==0 and fi1==np.pi*k+fi2) \
         or (th3==np.pi*k and fi1== fi2) or (th3==np.pi*k and fi1==np.pi*k+fi2):
        print('Попали в сингулярность, ищем значения слева и справа')
        Bmn_left=np.sin(thetamn-d_small*np.pi/180)*np.sin(((fi-d_small)-fimn)
        *np.pi/180)/np.sin((theta1-d_small)*np.pi/180)
        print('Beta_left2=',Bmn_left)
        Bmn_right=np.sin(thetamn+d_small*np.pi/180)*np.sin(((fi+d_small)-fimn)
        *np.pi/180)/np.sin((theta1+d_small)*np.pi/180)
        print('Beta_right2=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2                        
        return Bmn     
    elif th3==0 or th3==np.pi*k:
        print('Попали в сингулярность, ищем значения слева и справа')
        Bmn_left=np.sin((thetamn-d_small)*np.pi/180)*np.sin((fi-fimn)
        *np.pi/180)/np.sin((theta1-d_small)*np.pi/180)
        print('Beta_left3=',Bmn_left)
        Bmn_right=np.sin((thetamn+d_small)*np.pi/180)*np.sin((fi-fimn)
        *np.pi/180)/np.sin((theta1+d_small)*np.pi/180)
        print('Beta_right3=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2                        
        return Bmn
    elif fi1==fi2 or fi1==np.pi*k+fi2:
        print('Попали в сингулярность, ищем значения слева и справа')
        Bmn_left=np.sin(thetamn*np.pi/180)*np.sin(((fi-d_small)-fimn)
        *np.pi/180)/np.sin((theta1-d_small)*np.pi/180)
        print('Beta_left4=',Bmn_left)
        Bmn_right=np.sin(thetamn*np.pi/180)*np.sin(((fi+d_small)-fimn)
        *np.pi/180)/np.sin((theta1+d_small)*np.pi/180)
        print('Beta_right4=',Bmn_right)
        Bmn=(Bmn_left+Bmn_right)/2        
        return Bmn    
    else:
        print('Избежали попадания в сингулярность')
        Bmn=np.sin(thetamn*np.pi/180)*np.sin((fi-fimn)
            *np.pi/180)/np.sin(theta1*np.pi/180)        
        return Bmn
b=func_Bmn(th1,th2,th3,fi1,fi2)
print('Beta=',b)
print()

def func_ctgfi1 (theta,thetamn,fi,fimn):    
    if (th1==0 and fi1==fi2) or (th1==0 and fi1==np.pi*k+fi2)\
       or (th1==np.pi*k and fi1==fi2) or (th1==np.pi*k and fi1==np.pi*k+fi2): 
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos(thetamn*np.pi/180)*np.sin((theta-d_small)
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos((theta-d_small)*np.pi/180))/(np.sin((theta-d_small)
        *np.pi/180)*np.sin(((fi-fimn)-d_small)*np.pi/180))
        print('Ctgfi"_left1=',ctgfi1_left)
        ctgfi1_right=(np.cos(thetamn*np.pi/180)*np.sin((theta+d_small)
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos((theta+d_small)*np.pi/180))/(np.sin((theta+d_small)
        *np.pi/180)*np.sin(((fi-fimn)+d_small)*np.pi/180))
        print('Ctgfi"_right1=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1
    elif th1==0 or th1==np.pi*k:
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos(thetamn*np.pi/180)*np.sin((theta-d_small)
        *np.pi/180)*np.cos((fi-fimn)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos((theta-d_small)*np.pi/180))/(np.sin((theta-d_small)
        *np.pi/180)*np.sin((fi-fimn)*np.pi/180))
        print('Ctgfi"_left2=',ctgfi1_left)
        ctgfi1_right=(np.cos(thetamn*np.pi/180)*np.sin((theta+d_small)
        *np.pi/180)*np.cos((fi-fimn)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos((theta+d_small)*np.pi/180))/(np.sin((theta+d_small)
        *np.pi/180)*np.sin((fi-fimn)*np.pi/180))
        print('Ctgfi"_right2=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1
    elif fi1==fi2 or fi1==np.pi*k+fi2:
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos(thetamn*np.pi/180)*np.sin(theta
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos(theta*np.pi/180))/(np.sin(theta
        *np.pi/180)*np.sin(((fi-fimn)-d_small)*np.pi/180))
        print('Ctgfi"_left3=',ctgfi1_left)
        ctgfi1_right=(np.cos(thetamn*np.pi/180)*np.sin(theta
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)-np.sin(thetamn
        *np.pi/180)*np.cos(theta*np.pi/180))/(np.sin(theta
        *np.pi/180)*np.sin(((fi-fimn)+d_small)*np.pi/180))
        print('Ctgfi"_right3=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1
    elif (th1==np.pi*k/2 and fi1==np.pi*k/2+fi2)\
       or (th1==np.pi*k/2 and fi1==3*np.pi*k+fi2)\
       or(th1==3*np.pi*k/2 and fi1==np.pi*k/2+fi2)\
       or(th1==3*np.pi*k/2 and fi1==3*np.pi*k/2+fi2)\
       or(th3==0 and fi1==np.pi*k/2+fi2)\
       or(th3==np.pi*k and fi1==np.pi*k/2+fi2)\
       or(th3==0 and fi1==3*np.pi*k/2+fi2)\
       or(th3==np.pi*k and fi1==3*np.pi*k/2+fi2):
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos((thetamn-d_small)*np.pi/180)*np.sin((theta-d_small)
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)
        -np.sin((thetamn-d_small)*np.pi/180)*np.cos((theta-d_small)
        *np.pi/180))/(np.sin((theta-d_small)*np.pi/180)
        *np.sin(((fi-fimn)-d_small)*np.pi/180))
        print('Ctgfi"_left4=',ctgfi1_left)
        ctgfi1_right=(np.cos((thetamn+d_small)*np.pi/180)*np.sin((theta+d_small)
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)
        -np.sin((thetamn+d_small)*np.pi/180)*np.cos((theta+d_small)
        *np.pi/180))/(np.sin((theta+d_small)*np.pi/180)*np.sin(((fi-fimn)+d_small)*np.pi/180))
        print('Ctgfi"_right4=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1                            
    elif (th1==np.pi*k/2 and fi1==np.pi*k/2+fi2)\
         or (th1==np.pi*k/2 and fi1==3*np.pi*k+fi2)\
         or(th1==3*np.pi*k/2 and fi1==np.pi*k/2+i)\
         or(th1==3*np.pi*k/2 and fi1==3*np.pi*k/2+fi2):
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos(thetamn*np.pi/180)*np.sin((theta-d_small)
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)
        -np.sin(thetamn*np.pi/180)*np.cos((theta-d_small)
        *np.pi/180))/(np.sin((theta-d_small)*np.pi/180)
        *np.sin(((fi-fimn)-d_small)*np.pi/180))
        print('Ctgfi"_left5=',ctgfi1_left)
        ctgfi1_right=(np.cos(thetamn*np.pi/180)*np.sin((theta+d_small)
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)
        -np.sin(thetamn*np.pi/180)*np.cos((theta+d_small)
        *np.pi/180))/(np.sin((theta+d_small)*np.pi/180)*np.sin(((fi-fimn)+d_small)*np.pi/180))
        print('Ctgfi"_right5=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1
    elif (th3==0 and fi1==np.pi*k/2+fi2)\
       or(th3==np.pi*k and fi1==np.pi*k/2+fi2)\
       or(th3==0 and fi1==3*np.pi*k/2+fi2)\
       or(th3==np.pi*k and fi1==3*np.pi*k/2+fi2):
        print('Попали в сингулярность, ищем значения слева и справа')
        ctgfi1_left=(np.cos((thetamn-d_small)*np.pi/180)*np.sin(theta
        *np.pi/180)*np.cos(((fi-fimn)-d_small)*np.pi/180)
        -np.sin((thetamn-d_small)*np.pi/180)*np.cos(theta
        *np.pi/180))/(np.sin(theta*np.pi/180)
        *np.sin(((fi-fimn)-d_small)*np.pi/180))
        print('Ctgfi"_left6=',ctgfi1_left)
        ctgfi1_right=(np.cos((thetamn+d_small)*np.pi/180)*np.sin(theta
        *np.pi/180)*np.cos(((fi-fimn)+d_small)*np.pi/180)
        -np.sin((thetamn+d_small)*np.pi/180)*np.cos(theta
        *np.pi/180))/(np.sin(theta*np.pi/180)*np.sin(((fi-fimn)+d_small)*np.pi/180))
        print('Ctgfi"_right6=',ctgfi1_right)
        ctgfi1=(ctgfi1_left+ctgfi1_right)/2
        return ctgfi1                                                              
    else:
        print('Избежали попадания в сингулярность')
        ctgfi1=(np.cos(thetamn*np.pi/180)*np.sin(theta*np.pi/180)
                *np.cos((fi-fimn)*np.pi/180)-np.sin(thetamn*np.pi/180)
                *np.cos(theta*np.pi/180))/(np.sin(theta*np.pi/180)
                                           *np.sin((fi-fimn)*np.pi/180))
        return ctgfi1

ctg=func_ctgfi1(th1,th3,fi1,fi2)
print('CtgFi"=',ctg)
print()
Etmn=-(a*f1+b*f2*exp)
print('Etmn=',Etmn)
Etmn_abs=abs(Etmn)
print('Etmn_abs=',Etmn_abs)
print()
Efmn=b*f1-a*f2*exp
print('Efmn=',Efmn)
Efmn_abs=abs(Efmn)
print('Efmn_abs=',Efmn_abs)
print()



                
                

    

