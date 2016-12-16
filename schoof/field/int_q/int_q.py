r"""In this file required methods and classes for handling 
modular arithmetic are defined"""
from field import Field,FieldElement
class F_qField(Field):
    def __init__(self,q):
        self.q=q
    def add(self,el_1,el_2):
        val=el_1.val+el_2.val
        if ( val > self.q) :
            val-=self.q
        return self(val)
    def mul(self,el_1,el_2):
        val=el_1.val*el_2.val
        val%=self.q
        return self(val)
    def neg(self,el):
        return self(self.q-el.val)
    def eq(self,el_1,el_2):
        #print("EQ",el_1.val == el_2.val)
        #print("EQ",el_1.val,el_2.val)
        return el_1.val == el_2.val
    def __call__(self,val):
        return F_qElement(self,val)
    #(TODO) Faster way would be to use Extended Euclidean algorithm
    def invert(self,el):
        r"""Finds an inverse of a in F_q
        """
        return el^(self.q-2)
    def one(self):
        return self(1)
    def zero(self):
        return self(0)
class F_qElement(FieldElement):
    r"""Element of a finite field F_q, with q=p^d"""
    def __init__(self,field,val):
        self.field=field
        if ( val >= field.q or val < 0 ): # check this for speed
            val%=field.q
            if ( val< 0 ):
                val+=field.q
        self.val=val
    def __str__(self):
        return str(self.val)
    def __eq__(self,other):
        return self.field.eq(self,other)
# Greatest common divisor
#(TODO) Can be improved by eliminating recursion for 
# reducing stack memory usage
def gcd(a,b):
    if ( a < b ) :
        c=a
        a=b
        b=c
    r = a%b
    if ( r == 0 ):
        return b
    else:
        return gcd(b,r)

