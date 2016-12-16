r"""Implementation of Schoof's algorithm"""
import random 
from ring.polynomial import *
from field.int_q import *
from elliptic_curve import *
import copy


def schoof(a,b,p):
    r"""Schoofs algorithm. Counts the number of rational points 
    on an elliptic curve given by Weierstrass representation:

        y^2=x^3+ax+b mod p   ( with 4a^3+27b^2 != 0 mod p )

    in polynomial time complexity. 
    """
    # Finite field over which Elliptic curve is defined
    F=F_qField(p)
    polRing=PolynomialRing(F)
    a=F(a)
    b=F(b)
    if ( F(4) *( a^3) + F(27) * (b^2) == F.zero() ):
        print("Given a and b do not satisfy 4a^3+27b^2 != 0 mod (p)")
        return -1
   
    # Initialization of division polynomials \psi_l
    psi=[]
    psi.append(polRing([])) ## introduce \psi_l=0 so we have nice indexing
    psi.append(polRing([F.one()])) # \psi_1 = 1
    psi.append(polRing([F(2)])) # \psi_2 = Y * ( 2 )
    psi.append(polRing([-(a^2),F(12)*b,F(6)*a,F(0),F(3)])) # \psi_3 = -a^2  + 12bX +6aX^2+X^4
    psi.append(polRing([F(4)])*polRing([ \
        -(a^3) -F(8)*(b^2),-F(4)*a*b,-F(5)*(a^2),F(20)*b,F(5)*a,F(0),F(1)
        ])) # \psi_4 =  4Y(-a^3-8b^2-4abX-5a^2X^2+20bX^3+5aX^4+X^6)
    # test of trace_mod_l
    return trace_mod_l(a,b,p,5,psi) 

    # Erathostenes Sieve
    # TODO Find bound!!! here we will use primes up to million
    bound = 1000000
    #bound = 100
    primes=[True]*bound
    for i in range(2,bound):
        if ( primes[i] ) :
            trace_mod_l(a,b,p,i,psi)
            for j in range(i,bound):
                if ( i*j >= bound ):
                    break
                else:
                    primes[i*j]=False     
    return 0
def trace_mod_l(a,b,p,l,psi):
    r"""Calculate trace of Frobenius endomorphism mod l. 
    
    trace_mod_l(a,b,F,l) 

    "a" and "b" are coefficients defining elliptic curve in Weierstrass form: y^2=x^3+ax+b,
    j is a field over which a curve is defined, and l prime against which we calculate trace.
    """
    F=F_qField(p)
    polRing=PolynomialRing(F)
    q=p%l
    # polynomial describing Elliptic curve
    f=polRing([b,a,F(0),F(1)])
    # Find division polynomial \psi_l ( without Y which we can deduce from l)
    for i in range(len(psi),l+1):
        if i % 2 == 1 :  #\psi_{2m+1} = \psi_{m+2}\psi_m^3-\psi_{m+1}^3\psi_{m-1}
            m=i//2 # int division
            if m % 2 ==0 :
                pol=f*f*psi[m+2]*(psi[m]^3)-(psi[m+1]^3)*psi[m-1]
            else:
                pol=psi[m+2]*(psi[m]^3)-f*f*(psi[m+1]^3)*psi[m-1]
            psi.append(pol)
        else :  #2Y\psi_{2m} = \psi_{m}(\psi_{m+2}\psi_{m-1}^2-\psi_{m-2}\psi_{m+1}^2)
            m=i//2 # int division
            pol=psi[m]*(psi[m+2]*(psi[m-1]^2)-psi[m-2]*(psi[m+1]^2))
            pol=pol*polRing([F(2).__invert__()])
            psi.append(pol)
    h=psi[l]
    print(psi[5])
    return
    while ( True ) :
        try:
            #print(h)
            # we do arithmethics mod polynomial h
            polQuotRing=PolynomialQuotRing(polRing,h)
            # find \pi(X,Y)=(X^p,Y^p) and \pi^2=(X^2p,Y^2p)
            pi_1=[polQuotRing(polRing([F(0)]*p+[F(1)])),polQuotRing(f)^((p-1)//2)]
            #print(pi_1[0])
            #print("OVDE")
            pi_2=[pi_1[0]^p,(pi_1[1]^p)*pi_1[1]]
            # putting the points in EC context
            ecPol=ECPol(a,b,p,h)
            epi_1=ecPol(pi_1[0],pi_1[1])
            epi_2=ecPol(pi_2[0],pi_2[1])
            q_xy=ecPol(polQuotRing(polRing([F(0),F(1)])),polQuotRing.one())
            for i in range(1,q):
                q_xy=q_xy+q_xy
            lhs=q_xy+epi_2
            #print(lhs.x)
            # check zero!
            rhs=epi_1
            #print(rhs.x)
            for i in range(1,l):
                if ( rhs == lhs ):
                    return i
                else :
                    rhs = rhs+epi_1
        except PolInvError as err:
            h=err.gcd 
#print(schoof(31,-12,97))
print(schoof(1,2,257))
#print(schoof(184,896,1009))
