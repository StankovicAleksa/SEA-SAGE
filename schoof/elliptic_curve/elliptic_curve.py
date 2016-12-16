r"""In this file we define fundamental classes that 
deal with Elliptic curve arithmetic
"""
from modular_arithmetic import *

class EllipticCurve:
    r"""
    Elliptic Curve given by equation: y^2=x^3+ax+b mod q
    with a,b satisfying 4a^3+27b^2!=0 mod q
    """

    def __init__(self,a,b,q):
        self.a = a
        self.b = b
        self.q = q
        if ( 4*a*a*a+27*b*b % q == 0 ):
            raise Exception("Initialized Elliptic curve is singular")
    def create_point(self,x,y):
        return EllipticCurvePoint(x,y,self)
    def create_inf_point(self):
        return EllipticCurvePoint(self,0,0,inf=True)
    class EllipticCurvePoint:
        r"""
        Represents a point on a curve. Operations of addition and multiplication
        are implemented

        """
        def __init__(self,x,y,curve,inf=False):
            # init curve info
            self.curve=curve
            self.a=Fq_Element(self.curve.a,curve.q)
            self.b=Fq_Element(self.curve.b,curve.q)
            self.q=self.curve.q
            # init point info
            if ( inf == True):
                self.inf=True
            else:
                self.inf=False
                if ( type (x) is not Fq_Element ):
                    self.x=Fq_Element(x,curve.q)
                else:
                    self.x=x
                if ( type (y) is not Fq_Element ):
                    self.y=Fq_Element(y,curve.q)
                else:
                    self.y=y
                if ( (y*y-x*x*x-self.a*x-self.b) % q != 0 ):
                    raise Exception("Initialized point is not on an given elliptic curve")
        def __add__(self,other):
            # check for neutral points
            if ( self.inf == True ):
                return other
            elif(other.inf==True):
                return self
            else:
                if ( self.x != other.x):
                     s=(self.y-other.y)/(self.x-other.x)
                     x=s*s-self.x-other.x
                     y=self.y+s*(x-self.x)
                     y=-y
                     return EllipticCurvePoint(x,y,self.curve)
                else:
                    if ( self.y == -other.y):
                       return EllipticCurvePoint(0,0,self.curve,inf=True)  # Return inf point(neutral)
                    else:
                        s=(Fq_Element(3,self.q)*self.x*self.x+self.a)/(Fq_Element(2,self.q)*self.y)
                        x=s*s-self.x-self.x
                        y=self.y+s(x-self.x)
                        return EllipticCurvePoint(x,-y,self.curve)

