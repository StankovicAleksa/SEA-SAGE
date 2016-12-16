import numpy
import pickle
# Generates all Atkin polynomials Phi_l with derivatives w.r.t. x,j up to second order
# for l prime and l < 230

# Use builtin database from sage 
# For this we are required to install kohel_database with command
# $: sage -i database_kohel
AtkinDB = AtkinModularPolynomialDatabase()
Z_Xj.<var_X,var_j>=PolynomialRing(ZZ,'var_X,var_j')#PolynomialRing(F_p,'Xj') # construct polynomial ring

l=3
while ( l < 230 ):
    # open file for writing
    fl=open("atkin_pols/phi/pol_{}".format(l),'wb')
    fl_x=open("atkin_pols/phi_x/pol_{}".format(l),'wb')
    fl_xx=open("atkin_pols/phi_xx/pol_{}".format(l),'wb')
    fl_xj=open("atkin_pols/phi_xj/pol_{}".format(l),'wb')
    fl_j=open("atkin_pols/phi_j/pol_{}".format(l),'wb')
    fl_jj=open("atkin_pols/phi_jj/pol_{}".format(l),'wb')

    # convert to dictionary
    Phi=Z_Xj(AtkinDB[l])
    Phi_x=Phi.derivative(var_X)
    Phi_xx=Phi_x.derivative(var_X)
    Phi_xj=Phi_x.derivative(var_j)
    Phi_j=Phi.derivative(var_j)
    Phi_jj=Phi_j.derivative(var_j)
    pickle.dump(Phi.dict(),fl,pickle.HIGHEST_PROTOCOL)
    pickle.dump(Phi_x.dict(),fl_x,pickle.HIGHEST_PROTOCOL)
    pickle.dump(Phi_xx.dict(),fl_xx,pickle.HIGHEST_PROTOCOL)
    pickle.dump(Phi_xj.dict(),fl_xj,pickle.HIGHEST_PROTOCOL)
    pickle.dump(Phi_j.dict(),fl_j,pickle.HIGHEST_PROTOCOL)
    pickle.dump(Phi_jj.dict(),fl_jj,pickle.HIGHEST_PROTOCOL)
    l=next_prime(l)

