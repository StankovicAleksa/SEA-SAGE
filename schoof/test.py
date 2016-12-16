import ring
from ring.integer import Integer
from ring.polynomial import *
from field.int_q import F_qField,F_qElement
# test Integer
a=Integer(2)
b=Integer(3)
c=a+b
#print(c.val)
c=a+b
#print(c.val)
F_q=F_qField(1009)
a=F_q(3)
inv=a.__invert__()
print(inv)
print(a*inv)

#b=F_q(1)
#c=a+b
#print(c.val)
#pol_f=PolynomialRing(F_q)
#pol1=pol_f([F_q(3),F_q.one()])
#pol2=pol_f([F_q(3),a])
#pol3=pol1%pol2
#print(pol3)
#q_pol=pol_f([F_q.one(),F_q.one()])
#PQ=PolynomialQuotRing(pol_f,pol1)
#pol=PQ(pol2)
#print(pol.__invert__())
# Looks good
