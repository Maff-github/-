import numpy as np
import sympy as sp

theta = float(input ('Введите Тета='))
theta1= float(input('Введите Тета"='))
thetamn = float(input ('Введите Тетаmn='))
fi = float(input ('Введите Фи='))
fimn = float(input ('Введите Фиmn='))
f1= float(input('Введите нормированную ДН элемента 1='))
f2= float(input('Введите нормированную ДН элемента 2='))
beta=float(input('Введите коэффициент, определяющий поляризацию поля='))

Amn=(np.cos(theta*np.pi/180)*np.sin(thetamn*np.pi/180)
         *np.cos((fi-fimn)*np.pi/180)-np.sin(theta*np.pi/180)
         *np.cos(thetamn*np.pi/180))/np.sin(theta1*np.pi/180)
Bmn= np.sin(thetamn*np.pi/180)*np.sin((fi-fimn)
            *np.pi/180)/np.sin(theta1*np.pi/180)

if theta1==0:
    print('Предел справа для Amn =',sp.limit(Amn,theta1,0,'+'))
    print('Предел слева для Amn=',sp.limit(Amn,theta1,0,'-'))
    print('Предел для Amn=',sp.limit(Amn,theta1,0))
    
    print('Предел справа для Bmn =',sp.limit(Bmn,theta1,0,'+'))
    print('Предел слева для Bmn=',sp.limit(Bmn,theta1,0,'-'))
    print('Предел для Bmn=',sp.limit(Bmn,theta1,0))
    
    print('Излучатель расположен в полюсе антенны')
else:
    print('Amn=', Amn)
    print('Bmn=', Bmn)
    Etmn=-(Amn*f1+Bmn*f2*np.exp(-1j*beta))
    Efmn=Bmn*f1-Amn*f2*np.exp(-1j*beta)
    print('Etmn=', abs(Etmn))
    print('Efmn=', abs(Efmn))

   



















