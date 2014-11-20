#! /usr/bin/env python


class QuantumInt:
    """ the quantum amplitue propagator """
    def __init__(self):
        
        
        return
        
    # f' = f(t)
    # f(t) = 0.5 * t
    # f(t+\delta t) = f(t) + f'(t+0.5*\delta t) \delta t
    def test(self):
        n = 1000
        dt = 0.001
        t = 0.0
        fold = 0.0
        fgrad = 0.5 * (dt/2)# * (dt/2)
        fnew = fold + fgrad * dt
        for i in xrange(n):
            t = t + i*dt
            fold = fnew
            fgrad = 0.5 * (t+dt/2)#*(t+dt/2)
            fnew = fold + fgrad * dt
            print t, fnew
        return

        
#include "stdio.h" 

#include "conio.h"  

float func(float x,float y) 

{  

return(2*x*y);  

}  

float runge_kutta(float x0,float xn,float y0,int n) 

{  

float x,y,y1,y2,h,xh; 

float d1,d2,d3,d4; 

int i; 

x=x0;  

y=y0;  

h=(xn-x0)/n; 

for(i=1;i<=n;i++)  
        {  

xh=x+h/2; 

d1=func(x,y);  

d2=func(xh,y+h*d1/2.0); 

d3=func(xh,y+h*d2/2.0); 

d4=func(xh,y+h*d3);  

y=y+h*(d1+2*d2+2*d3+d4)/6.0; 

x=x0+i*h; 

} 

return(y); 

}  

void main() 

{  

float x0,xn,y0,e; 

int n;  

printf("\ninput n:\n"); 

scanf("%d",&n);  

printf("input x0,xn:\n"); 

scanf("%f%f",&x0,&xn); 

printf("input y0:\n"); 

scanf("%f",&y0);  

e=runge_kutta(x0,xn,y0,n); 

printf("y(%f)=%6.6f",y0,e); 

}  

# http://zh.wikipedia.org/wiki/%E9%BE%99%E6%A0%BC%EF%BC%8D%E5%BA%93%E5%A1%94%E6%B3%95        
        
# main program
if __name__ == "__main__":
    q = QuantumInt()
    q.test()
    
    