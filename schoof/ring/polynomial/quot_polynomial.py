r"""Class that handles arithmetic on polynomials F[X]/h(X), with F 
being a given field and h a polynomial over F.
"""
import copy
from ring import Ring,RingElement
class PolynomialQuotRing(Ring):
    def __init__(self,polRing,q_pol):
        self.polRing=polRing
        self.q_pol=q_pol
    def add(self,el_1,el_2):
        pol=(el_1.pol+el_2.pol)%self.q_pol
        return PolynomialQuot(self,pol)
    def mul(self,el_1,el_2):
        pol=(el_1.pol*el_2.pol)%self.q_pol
        return PolynomialQuot(self,pol)
    def neg(self,el):
        pol=(-el.pol)
        return PolynomialQuot(self,pol)
    def div(self,el_1,el_2):
        return self.mul(el_1,self.invert(el_2))
    # Extended Euclidean algorithm
    def invert(self,el):
        b=el.pol
        n=self.q_pol
        x0, x1, y0, y1 = self.polRing.one(), self.polRing.zero(), self.polRing.zero(), self.polRing.zero() 
        while len(n)>0:
            q, b, n = b / n, n, b % n
            x0, x1 = x1, x0 - q * x1
            y0, y1 = y1, y0 - q * y1
        if ( len(b)!= 1 ):
            raise PolInvError(b)
        else:
            x0=x0*self.polRing([b[-1].__invert__()])
            return PolynomialQuot(self,x0)
    def equal(self,el_1,el_2):
        return el_1.pol == el_2.pol
    def one(self):
        return PolynomialQuot(self,self.polRing.one())
    def zero(self):
        return PolynomialQuot(self,self.polRing.zero())
    def __call__(self,arg):
        pol=arg
        if (type(arg) == type([])  ):
            pol=polRing(arg)
        return PolynomialQuot(self,pol)
class PolynomialQuot(RingElement):
    def __init__(self,polQuotRing,pol):
        r"""Initialize polynomial over field with given
        no_variables(in case no_var=1, we have P=F_q[X],
        in case no_var=2, we have P=F_q[X,Y], etc.).
        """
        RingElement.__init__(self,polQuotRing)
        self.pol=pol%polQuotRing.q_pol
        self.polQuotRing=polQuotRing
    def __getitem__(self,i):
        return self.pol[i]
    def __len__(self):
        return len(self.pol)
    def __mod__(self,other):
        return self.polQuotRing.mod(self,other)
    def __str__(self):
        return str(self.pol)
    def __invert__(self):
        return self.polQuotRing.invert(self)
    def __truediv__(self,other):
        return self.polQuotRing.div(self,other)
class PolInvError(Exception):
    def __init__(self,gcd):
        self.gcd=gcd
    def __str__(self):
        return "Polynomial has no inverse"
