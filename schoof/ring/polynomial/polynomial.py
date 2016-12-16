r"""Class that handles arithmetic on polynomials F[X], where F is given field
"""
import copy
from ring import Ring,RingElement
class PolynomialRing(Ring):
    def __init__(self,field):
        self.field=field
    def add(self,el_1,el_2):
        if ( len(el_1.coeff) < len(el_2.coeff)):
            c=el_1
            el_1=el_2
            el_2=c
        coeff=copy.copy(el_1.coeff)
        i=0
        for ar_el in el_2.coeff:
            coeff[i]+=ar_el
            i=i+1
        l=len(coeff)-1

        while l >= 0 and coeff[l] ==  self.field.zero() :
            del coeff[l]
            l=l-1
            
        return Polynomial(self,coeff)
    def mul(self,el_1,el_2):
        ar=[]
        # deg 
        deg=len(el_1)+len(el_2)-2
        for i in range(0,deg+1): # 0<= i <= deg
            # create c_ix^i
            val=self.field.zero()
            for j in range(0,i+1):
                if j < len(el_1) and i-j < len(el_2) :
                    val+=el_1[j]*el_2[i-j]
            ar.append(val)
        return Polynomial(self,ar)
    def neg(self,el):
        ar=copy.copy(el.coeff)
        for i in range(0,len(ar)):
            ar[i]=-ar[i]
        return Polynomial(self,ar)
    #TODO
    def mod(self,el_1,el_2):
        ar=copy.copy(el_1.coeff)
        inv=el_2[-1].__invert__()
        l1=len(el_1)
        l2=len(el_2)
        i=l1-1
        while ( i - l2+1 >= 0   ) :
            for j in range(0,l2):
                ar[i-l2+1+j]=ar[i-l2+1+j]-ar[i]*inv*el_2[j]
            i=i-1 
        while len(ar)>0 and  ar[-1] == self.field.zero() :
            del ar[-1]
        return Polynomial(self,ar)
    def div(self,el_1,el_2):
        if ( len(el_2) == 0 ):
            raise Exception("Dividing by 0")
        ar=copy.copy(el_1.coeff)
        ar_1=[]
        inv=el_2[-1].__invert__()
        l1=len(el_1)
        l2=len(el_2)
        i=l1-1
        while ( i - l2+1 >= 0   ) :
            ar_1.append(ar[i]*inv)
            for j in range(0,l2):
                ar[i-l2+1+j]=ar[i-l2+1+j]-ar[i]*inv*el_2[j]
            i=i-1 
        ar_1.reverse()
        return Polynomial(self,ar_1)
    def equal(self,el_1,el_2):
        if ( len(el_1) != len (el_2) ):
            return False
        for i in range(0,len(el_1)):
            if ( el_1[i] != el_2[i]):
                return False
        return True
    def __call__(self,coeff):
        return Polynomial(self,coeff)
    def one(self):
        return Polynomial(self,[self.field.one()])
    def zero(self):
        return Polynomial(self,[self.field.zero()])
class Polynomial(RingElement):
    def __init__(self,polRing,coeff):
        r"""Initialize polynomial over field with given
        no_variables(in case no_var=1, we have P=F_q[X],
        in case no_var=2, we have P=F_q[X,Y], etc.).
        """
        RingElement.__init__(self,polRing)
        self.coeff=coeff
        self.polRing=polRing
    def __getitem__(self,i):
        return self.coeff[i]
    def __len__(self):
        return len(self.coeff)
    def __mod__(self,other):
        return self.polRing.mod(self,other)
    def __str__(self):
        if ( len(self.coeff) == 0 ):
            return "[]"
        ret="["
        for i in self.coeff:
            ret=ret+str(i)+","
        return ret[:-1]+"]"
    def __truediv__(self,other):
        return self.polRing.div(self,other)
