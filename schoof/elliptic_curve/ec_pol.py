from ring.polynomial import *
from field.int_q import *
class ECPol:
    def __init__(self,a,b,p,h):
       self.a=a
       self.b=b
       self.F=F_qField(p)
       self.polRing=PolynomialRing(self.F)
       self.h=h
       self.polQuotRing=PolynomialQuotRing(self.polRing,h)
       self.f=self.polQuotRing(self.polRing([b,a,self.F(0),self.F(1)]))
    def add(self,el_1,el_2):
        if ( el_1.x == el_2.x ):
            n2=self.polQuotRing(self.polRing([self.F(2)]))
            n3=self.polQuotRing(self.polRing([self.F(3)]))
            q_a=self.polQuotRing(self.polRing([self.a]))
            r=(n3*(el_1.x^2)+q_a)/(n2*el_1.y*self.f)
        else:
            r=(el_1.y-el_2.y)/(el_1.x-el_2.x)
        x=r*r*self.f-el_1.x-el_2.x
        y=r*(el_1.x-x)-el_1.y
        return self(x,y)
    def equal(self,el_1,el_2):
        return el_1.x == el_2.x and el_1.y==el_2.y
    def __call__(self,x,y):
        return ECPolPoint(self,x,y)
class ECPolPoint:
    def __init__(self,ECPol,x,y):
        self.x=x
        self.y=y
        self.ECPol=ECPol
    def __add__(self,other):
        return self.ECPol.add(self,other)
    def __eq__(self,other):
        return self.ECPol.equal(self,other)
