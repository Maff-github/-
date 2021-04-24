import numpy as np
import sympy as sp
def f1(x1):
    return x1
x1=float(input('Введите нормированную амплитуду 1(ДН элемента)='))
print('f1=',f1(x1))
def f2(x2):
    return x1
x2=float(input('Введите нормированную амплитуду 2(ДН элемента)='))
print('f2=',f1(x2))    
def e (beta):
    return np.exp(-1j*beta)
beta1=float(input('Введите коэффициент, определяющий поляризацию поля='))
exp=e(beta1)
print('exp=',exp)
def Amn(theta,theta1,thetamn,fi,fimn):
    return (np.cos(theta*np.pi/180)*np.sin(thetamn*np.pi/180)
                *np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
                *np.cos(thetamn*np.pi/180))/np.sin(theta1*np.pi/180)
th1=float(input ('Введите Тета='))
th2=float(input('Введите Тета"='))
th3=float(input ('Введите Тетаmn='))
fi1=float(input ('Введите Фи='))
fi2=float(input ('Введите Фиmn='))
if th2==0:
    print('Предел справа для Amn =',sp.limit(Amn(th1,th2,th3,fi1,fi2),th2,0,'+'))
    print('Предел слева для Amn=',sp.limit(Amn(th1,th2,th3,fi1,fi2),th2,0,'-'))
    print('Предел для Amn=',sp.limit(Amn(th1,th2,th3,fi1,fi2),th2,0))
    print('Излучатель расположен в полюсе антенны')
else:
    Anm1=Amn(th1,th2,th3,fi1,fi2)
    print('Amn=',Anm1)      
def Bmn (theta,theta1,thetamn,fi,fimn):
    return np.sin(thetamn*np.pi/180)*np.sin((fi-fimn)
            *np.pi/180)/np.sin(theta1*np.pi/180)
if th2==0:
    print('Предел справа для Bmn =',sp.limit(Bmn(th1,th2,th3,fi1,fi2),th2,0,'+'))
    print('Предел слева для Bmn=',sp.limit(Bmn(th1,th2,th3,fi1,fi2),th2,0,'-'))
    print('Предел для Bmn=',sp.limit(Bmn(th1,th2,th3,fi1,fi2),th2,0))
    print('Излучатель расположен в полюсе антенны')
else:
    Bmn1=Bmn(th1,th2,th3,fi1,fi2)
    print('Bmn=',Bmn1)
    Etmn=abs(-(Amn(th1,th2,th3,fi1,fi2)
           *f1(x1)+Bmn(th1,th2,th3,fi1,fi2)*f2(x2)*e(beta1)))
    print('Etmn=',Etmn)
    Efmn=abs(Bmn(th1,th2,th3,fi1,fi2)*f1(x1)
          -Amn(th1,th2,th3,fi1,fi2)*f2(x2)*e(beta1))
    print('Efmn=',Efmn)








 
    

    
    
    
                                                
    

    
               
                    

    

    
        
        
        
        
        
            
        

       
 









    
    
 











