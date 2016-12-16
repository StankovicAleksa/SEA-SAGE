# Count points on elliptic curve if the order is prime, if not 
# break and return unsafe warning
import numpy
def count_points(a,b,p):
    if not is_prime(p):
        raise ValueError("Given p={} is not a prime!".format(p))
    # Construct the field
    F_p=GF(p)
    a=F_p(a)
    b=F_p(b)
    disc= 4*a**3 + 27*b**2  # Discriminant of the elliptic curve
    if ( disc== 0 ):
        raise ValueError("Elliptic curve is singular")
    
    # Initialize Sagemath elliptic curve
    El_c=EllipticCurve(GF(p),[0,0,0,a,b])

    # calculate invariant j
    j=1728*(4*(a**3))/disc
    F_pX=PolynomialRing(F_p,'var_X')
    F_pj=PolynomialRing(F_p,'var_j')

    fc=F_pX("var_X^3+{}*var_X+{}".format(a,b))
    #global pol
    #pol=F_pX("3*var_X^4+6*var_X^2+7*var_X+16")
    #var_X=var('var_X')
    #var_j=var('var_j')
    F_pXj.<var_X,var_j>=PolynomialRing(F_p,'var_X,var_j')#PolynomialRing(F_p,'Xj') # construct polynomial ring
    # analyse trace of Frobenius trace against different primes l
    
    # set of results for Elkies primes and isogeny cycles
    E=[]
    # set of results for Atkin primes
    A=[]
    # use isogeny cycles for l <= 13 
    # sea for l>13 
    # Database of Atkin polynomials

    control_bound=pow(p,0.5)
    gathered_info=1

    # solve for l = 2 
    pq=fast_power(F_pX("1"),F_pX("var_X"),p,fc)-F_pX("var_X")
    if ( pq.gcd(fc).is_one() ):
        E.append([2,1])
    else:
        E.append([2,0])
        return -1
    l=2
    AtkinDB = AtkinModularPolynomialDatabase()
    res=1
    while ( True ):
        # make sure we work with l \in F_p
        # do the iteration for next prime
        l=next_prime(l)
        # No more Atkin polynomials available
        if ( l > 200 or gathered_info > control_bound) :
            break

        gathered_info*=l
        # Extract Atkin polynomial from the database 
        Phi=F_pXj(AtkinDB[l])
        X_Phi=F_pX(Phi(var_j=j))
        # calculate x^p mod X_Phi using fast exp
        x=F_pX("var_X")
        x_p=fast_power(F_pX("1"),x,p,X_Phi)
        x1=x_p-x
        
        
        # Find gcd(X_phi,x1) so that we can use Atkin classification
        G=x1.gcd(X_Phi)
        if ( G.is_constant() ): 
            #############################################
            # Procedure when l is Atkin prime
            #############################################
            print("Atkin prime {}".format(l))
            T_l=[]
            # We will calculate r as r=(l+1)/s, where s is nubmer of irreducible polynomials in G

            # We proceed by following approach with Berlekamp matrices
            x_ip=F_pX("1")
            rows=[]

            import copy
            for i in range(0,l+2):
                if ( i != 0 ):
                    x_ip=(x_ip*x_p)%X_Phi
                coef=x_ip.coefficients()
                ## extend length of vector to l in order to match matrix size
                while len(coef) != l+2 :
                    coef.append(F_p(0))
                coef[i]=coef[i]-1
                rows.append(coef)
            M=matrix(rows,ring=F_p)
            # debugging 
            #s=len(X_Phi.factor())
            s=l+2-M.rank()
            #print("{} {}".format(len(X_Phi.factor()),l+2-M.rank()))
            r=(l+1)/s
            # possible fix
            #if ( r % 2 == 0 ) :
                #r*=2
            # 
            l2=l**2
            F_l2=GF(l2)  # field F_{l^2}
            m_gen=F_l2.multiplicative_generator()
            # primitive r-th root of unity of F_l
            gamma=m_gen**((l2-1)/r)
            gamma_i=1
            ### We need only half of the values of r thanks to the symmetry
            for i in range(0,r//2):
                gamma_i=gamma_i*gamma
                tmp=gamma_i+2+(gamma_i**(-1)) # TODO solve and speed up the algorithm!!! 
                rhs=p*tmp.polynomial()[0]
                if ( rhs.is_square()) :
                    p_1=rhs.sqrt()
                    T_l.append(Integer(p_1))
                    if ( p_1 != 0 ):
                        T_l.append(Integer(-p_1))
            #T_l.sort()
            A.append([l,T_l])
            # Escape strategy!!! TODO
        else:                   
            #############################################
            # Procedure when l is Elkies prime
            #############################################
            # Just Elkies case is left 
            # search Elkies polynomial
            # We calculate l-th division polynomial
            # TODO probably can be faster !!!
            print("Elkies prime {}".format(l)) 
            m = 18 *b / a
            j_p = m * j
            k=j_p/(1728-j) 
            # uhm
            f_E=G.any_root()
            ## We now calculate derivatives of polynomial Phi with respect to both variables 
            j_Phi=F_pj(Phi(var_X=(f_E))) # Phi evaluated at f_E


            j_Phi_j=j_Phi.derivative()
            Phi_j=j_Phi_j(var_j=j)
            j_Phi_jj=j_Phi_j.derivative()
            Phi_jj=j_Phi_jj(var_j=j)
            p_Phi=Phi.derivative(var_X)
            pj_Phi=F_pj(p_Phi(var_X=f_E))
            pj_Phi_j=pj_Phi.derivative()
            Phi_fj=pj_Phi_j(var_j=j)
            pf_Phi=F_pX(p_Phi(var_j=j))
            Phi_f=pf_Phi(var_X=f_E)
            Phi_ff=pf_Phi.derivative()(var_X=f_E)
            ## Calculate f'
            f_p=-j_p*Phi_j/Phi_f
            # checks!!
            if ( l > 200 ):
                print (j_Phi)
                return

            continue;
            roots=j_Phi.roots() ## TODO probably can be faster!
            for root in roots: # try roots
                tj=root[0]
                # Calculate \tilde(\Phi_j),\tilde(\Phi_jf),...
                tPhi_j=j_Phi_j(var_j=tj)
                tPhi_jj=j_Phi_jj(var_j=tj)

                tPhi_fj=pj_Phi_j(var_j=tj)
                f_tPhi=F_pX(Phi(var_j=tj))
                f_tPhi_f=f_tPhi.derivative()
                tPhi_f=f_tPhi_f(var_X=f_E)
                tPhi_ff=f_tPhi_f.derivative()(var_X=f_E)
                
                # Computation of polynomials is finished
                tj_p=j_p/l * Phi_j/Phi_f * tPhi_f / tPhi_j
                tm=tj_p/tj
                tk=tj_p/(1728-tj)

                # Probably incorrect
                #ta=tm*tk/48
                #tb=(tm**2)*tk/864
                ta=(F_p(l)**4)*tm*tk/48
                tb=(F_p(l)**6)*(tm**2)*tk/864


                tr_f=-((f_p**2)*tPhi_ff+2*l*f_p*tj_p*tPhi_fj+(l**2)*(tj_p**2)*tPhi_jj)/(f_p*tPhi_f)
                r=tr_f-(j_p*Phi_jj+2*f_p*Phi_fj+(f_p**2)*Phi_ff/j_p)/(Phi_j)
                p_1=l*(r/2+(k-l*tk)/4+(l*tm-m)/3)
                d=(F_p(l)-1)/2
                t=numpy.empty([max(d+1,4)],dtype=type(j))
                t[0]=d
                t[1]=p_1/2
                t[2]=((1-10*d)*a-ta)/30
                t[3]=((1-28*d)*b-42*t[1]*a-tb)/70
                c=numpy.empty([max(d+1,3)],dtype=type(j))
                c[0]=F_p(0)
                c[1]=6*t[2]+2*a*t[0]
                c[2]=10*t[3]+6*a*t[1]+4*b*t[0]
                for i in range(2,d):
                    s=numpy.dot(c[0:i+1][::-1],c[0:i+1])
                    # formula not checked !!! 
                    c[i+1]=(3*s-(2*i-1)*(i-1)*a*c[i-1]-(2*i-2)*(i-2)*b*c[i-2])/((i-1)*(2*i+5))
                for i in range(3,d):
                    t[i+1]=(c[i]-(4*i-2)*a*t[i-1]-(4*i-4)*b*t[i-2])/(4*i+2)
                s=numpy.empty([d+1],dtype=type(j))
                s[0]=F_p(1)
                for n in range(1,d+1):
                    s[n]=F_p(0)
                    c=F_p(-1)
                    for i in range ( 1,n+1):
                        s[n]+=c*t[i]*s[n-i]
                        c=-c
                    s[n]=s[n]*F_p(-1)/n
                c=1
                for n in range(0,d+1): 
                    if ( c == 0 ) :
                        s[n]=-s[n]
                    c=(c+1)%2

                # Elkies polynomial is phi
                phi=F_pX(s[::-1].tolist())
               
                # construct division polynomials mod phi
                PQR=F_pX.quotient(phi)
                #PQR=F_pX
                x_pqr=PQR('var_X')
                div=numpy.empty([max(l+1,5)],dtype=type(x)) # div[i]: i-th division polynomial mod phi 
                div2=numpy.empty([max(l+1,5)],dtype=type(x)) # div2[i]: div[i]^2  mod phi
                div02=numpy.empty([max(l+1,5)],dtype=type(x)) # div02[i]: div[i]*div[i+2] mod phi
                # initialize starting div polynomials
                div[0]=PQR('0')
                div[1]=PQR('1')
                div[2]=PQR('1')
                div[3]=PQR('3*var_X^4+{}*var_X^2+{}*var_X+{}'.format(6*a,12*b,-a*a))
                div[4]=PQR("2*(var_X^6+{}*var_X^4+{}*var_X^3+{}*var_X^2+{}*var_X+{})".format(5*a,20*b,-5*a*a,-4*a*b,-(a**3)-8*b*b))

                f=PQR("var_X^3+{}*var_X+{}".format(a,b))
                # init div2
                for i in range(0,5):
                    div2[i]=div[i]*div[i]

                # init div02
                for i in range(0,3):
                    div02[i]=div[i]*div[i+2]
                #
                for i in range(5,l+1):
                    im=i//2
                    if ( i%2 ==1 ):
                        if ( im%2 == 0 ):
                            pol=16*f*f*div02[im]*div2[im]-div02[im-1]*div2[im+1]
                        else:
                            pol=div02[im]*div2[im]-16*f*f*div02[im-1]*div2[im+1]
                    else:
                        pol=(div02[im]*div2[im-1]-div02[im-2]*div2[im+1])
                    div[i]=pol
                    div2[i]=pol*pol
                    div02[i-2]=div[i-2]*div[i]
                if (div[l]==0):
                    x_p=x_pqr**(p)
                    y_p=f**((p-1)/2)
                    # check t_l==0
                    lam=0
                    # point at infinity
                    if ( x_p-x_pqr == 0 ):
                        if ( y_p == 1 ):
                            lam=1
                        else:
                            lam=l-1
                    else:
                        # do all checks
                        for i in range(2,(l-1)//2+1):
                            # TODO optimize
                            if ( i% 2 == 0 ):
                                pol_xi=x_p*div2[i]*4*f-x_pqr*div2[i]*4*f+div02[i-1]
                            else:
                                pol_xi=x_p*div2[i]-x_pqr*div2[i]+div02[i-1]*4*f
                            if ( pol_xi == 0 ):
                                # check the formulas
                                if ( i%2 ==0 ):
                                    pol_yi=4*f*y_p  *f*8*div2[i]*div[i]    -2*(   div[i+2]*div2[i-1]-div[i-2]*div2[i+1])
                                else:
                                    pol_yi=4*f*y_p      *div2[i]*div[i]    -4*f*(div[i+2]*div2[i-1]-div[i-2]*div2[i+1])
                                if ( pol_yi==  0 ):
                                    lam=i
                                else:
                                    lam=l-i
                                break
                    t_l=Integer(lam+GF(l)(p)/lam)
                    # Early-Abort strategy
                    if ( (p+1-t_l) % l == 0 or (p+1+t_l) %l ==0 ):
                        return -1
                    E.append([l,t_l])
                    break
                    # Indeed good
                    # Proceed to Schoof with phi
    # When to escape it??
    ################################################################
    #
    # Match-sort algorithm
    #
    ################################################################
    
    ################################################
    # Gather information from Elkies steps
    ################################################
    # check also what if there are no Elkies steps
    mE=1 
    # chinese remainder theorem
    l_set=[]
    l_mods=[]
    for i in E:
        l_set.append(i[0])
        l_mods.append(i[1])
        mE*=i[0]
    tE=crt(l_mods,l_set)
    


    ################################################
    # Gather information from Atkin steps
    ################################################
    A1=[]
    A2=[]
    sz1=1
    sz2=1
    m1=1
    m2=1
    # split A into sets A1 and A2
    for A_l in A:
        if ( sz1*A_l[0] > 100000 or sz2*A_l[0] > 100000 ):
            break
        if ( sz1 <= sz2) :
            A1.append(A_l)
            sz1*=len(A_l[1])
            m1*=A_l[0]
        else :
            A2.append(A_l)
            sz2*=len(A_l[1])
            m2*=A_l[0]
    # Special case when there is only one or none Atkin primes l.
    # In this case BS-GS is not applicable.
    if ( len(A) == 1 ):
        sz2=0
    elif ( len(A) == 0 ):
        sz1=0
        sz2=0

    # construct sets R_1 and R_2
    r1_set=numpy.empty([sz1],dtype=Integer)
    r2_set=numpy.empty([2*sz2],dtype=Integer)
    t_set1=numpy.empty([max(sz1,2*sz2)],dtype=Integer) # auxilliary set for constructing r1_set and r2_set
    t_set2=numpy.empty([max(sz1,2*sz2)],dtype=Integer) # auxilliary set for constructing r1_set and r2_set
    m1mE=m1*mE
    m2mE=m2*mE

    
    ################################################
    # Construct R_1
    ################################################

    #initialization
    src_sz=1
    dest_sz=0
    t_set1[0]=1
    t_m1=1
    switch = True
    for el in A1: 
        dest_sz=0
        if ( switch ):
            src=t_set1
            dest=t_set2
        else:
            src=t_set2
            dest=t_set1

        for el_a in el[1]:
            for i in range(0,src_sz):
                dest[dest_sz]=crt(src[i],el_a,t_m1,el[0])
                dest_sz+=1

        switch=not switch
        src_sz=dest_sz
        t_m1*=el[0]
    if ( switch ):
        src=t_set1
    else:
        src=t_set2
    m1Ring=IntegerModRing(m1)
    m2mE_inv1=1/m1Ring(m2mE)
    bound1=(m1-1)/2
    for i in range (0,sz1):
        # make sure |r1_set[i]|<=(m1-1)/2
        r1_set[i]=Integer((src[i]-m1Ring(tE))*m2mE_inv1)
        if ( r1_set[i] < 0 ):
            r1_set[i]+=m1
        if ( r1_set[i]> bound1):
            r1_set[i]-=m1
    r1_set.sort()
    ################################################
    # Construct R_2
    ################################################

    #initialization
    src_sz=1
    dest_sz=0
    t_set1[0]=1
    t_m2=1
    switch = True
    for el in A2: 
        dest_sz=0
        if ( switch ):
            src=t_set1
            dest=t_set2
        else:
            src=t_set2
            dest=t_set1

        for el_a in el[1]:
            for i in range(0,src_sz):
                dest[dest_sz]=crt(src[i],el_a,t_m2,el[0])
                dest_sz+=1

        switch=not switch
        src_sz=dest_sz
        t_m2*=el[0]
    if ( switch ):
        src=t_set1
    else:
        src=t_set2
    
    m2Ring=IntegerModRing(m2)
    m1mE_inv2=1/m2Ring(m1mE)

    for i in range (0,sz2):
        r2_set[i]=Integer((src[i]-m2Ring(tE))*m1mE_inv2)
        r2_set[i]%=m2
        if (r2_set[i]<0 ):
            r2_set[i]+=m2
    for i in range(0,sz2):
        r2_set[i+sz2]=r2_set[i]-m2
    r2_set.sort()
    # Match-sort algorithm
    # search for solution

    # run for multiple points
    enough_iter=10
    n_iter=0
    Q1=numpy.empty([sz1],dtype=type(El_c.random_point()))
    #Q1=numpy.empty([sz1],dtype=Integer)
    while( n_iter < enough_iter):
        # Find random point P
        P=El_c.random_point()
        Q=(p+1-tE)*P
        m1mE_P=m1mE*P
        m2mE_P=m2mE*P
        # CHECK ATKIN STEPS!
        # Atkin steps were not required!
        if (Q == 0 ):
            print("ATKIN NOT NECESSARY ")
            return p+1-tE
        for i in range(0,len(r1_set)):
            #Q1[i]=r1_set[i]*m2
            Q1[i]=(Q-r1_set[i]*m2mE_P)
            Q1[i].r1=r1_set[i]
        Q1.sort()
        for i in range(0,len(r2_set)):
            index=binary_search(Q1,r2_set[i]*m1mE_P)
            if ( index!=-1):
                r1=Q1[index].r1
                r2=r2_set[i]
                return p+1-(tE+mE*(r1*m2+r2*m1))
        n_iter+=1
        
# calculate a^b mod m
def fast_power(one,a,b,m):
    if ( b== 0 ) :
        return one 
    r=fast_power(one,a,b//2,m)
    r=(r*r)%m
    if ( b%2 == 1 ):
        r=r*a
        r%=m
    return r
def binary_search(array,elem):
    low=0
    high=len(array)-1
    while ( low < high ):
       mid=(low+high)//2
       if array[mid] < elem:
           low=mid+1
       elif array[mid] == elem:
           return mid
       else:
           high=mid-1
    return -1
        
       
import time
a=231231231252423323213154362654220789709865
b=23107312023123123313213123124523452345
p=next_prime(2231238224987808690822570987327809797089711231)
## Elliptic curve is of the form y^2=x^3+ax+b
def do_test(a,b,p):
    start=time.time()
    print("my count: {}".format(count_points(a,b,p)))
    print(time.time()-start)
    start=time.time()
    E=EllipticCurve(GF(p),[0,0,0,a,b])
    print("candidate: {}".format(E.cardinality()))
    print(time.time()-start)
do_test(a,b,p)
