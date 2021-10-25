reset()
load('constants.sage')
load('simultaneous_modular_operations.sage')
load('matrix_of_relations.sage')
load('change_of_basis.sage')
load('modular_composition_basecase.sage')

# Some small batches of tests:
# to active test, set to True; to disable, set to False
SimModOpsTests=True # simultaneous modular operations
MatRelTests=True # matrix of relations
ChBasTests=True # change of basis
ModCompTests=True # modular composition base case
# (complete test, with all set to True, should not take much more than a minute
# except on old or slow architectures)

if any([SimModOpsTests,MatRelTests,ChBasTests,ModCompTests]):
    print("Launching tests. If nothing prints below, this means that all tests passed.")
    print("Note: failure because of an unlucky random choice is very unlikely, but can happen.")
    K = GF(65537)
    # using double letters for variables
    # to make sure we spot any potential undefined (single letter) x / y / X / Y in functions
    Kxx.<xx> = K[]
    Kyy.<yy> = K[]
    Kxxyy.<XX,YY> = K[]

if SimModOpsTests:
    # Brent Kung
    for n in [10,50,200]:
        f = Kxx.random_element(degree=n)
        a = Kxx.random_element(degree=n-1)
        for d in [n//2,n-1,n+n//2]:
            g = Kyy.random_element(degree=d-1)
            b = ModularComposition_BrentKung(f,a,g,d)
            if b != (g(a) % f):
                print("ModularComposition_BrentKung Fail with",n,d)

    # Small minpoly (actually not small)
    g = Kyy.random_element(degree=n-1)
    b = ModularComposition_SmallMinimalPolynomial(f,a,g,n,[K.random_element() for i in range(n)])
    if b != (g(a) % f):
        print("ModularComposition_SmallMinimalPolynomial(notsmall) Fail with",n,d)

    # Small minpoly (actually small)
    f2 = xx**n
    a = xx**(n/2)
    b = ModularComposition_SmallMinimalPolynomial(f2,a,g,2,[K.random_element() for i in range(n)])
    if b != (g(a) % f2):
        print("ModularComposition_SmallMinimalPolynomial(notsmall) Fail with",n,d)

    # SimultaneousBivariateModularComposition
    f = Kxx.random_element(degree=n)
    a = Kxx.random_element(degree=n-1)
    m = 5
    r = ceil(n/m)
    s = 6
    gg = []
    for i in range(s):
        gg.append(Kxxyy.random_element(degree=m+r,terms=Infinity))
        gg[i] = gg[i].truncate(XX,m).truncate(YY,r)
    bb = SimultaneousBivariateModularComposition(f,a,gg,m,r)
    if any([bb[i] != gg[i](xx,a) % f for i in range(s)]):
        print("SimultaneousBivariateModularComposition Fail with",n,m,r,s)

    # BivariateModularComposition
    for n in [10,50,200]:
        f = Kxx.random_element(degree=n)
        a = Kxx.random_element(degree=n-1)
        for m in [5,10,30]:
            for d in [5,10,30]:
                g = Kxxyy.random_element(degree=m+d,terms=Infinity)
                g = g.truncate(XX,m).truncate(YY,d)
                b = BivariateModularComposition(f,a,g,m,d)
                if b != (g(xx,a) % f):
                    print("BivariateModularComposition Fail with",n,m,d)

    # SimultaneousTruncatedModularMultiplication
    for n in [10,50,200]:
        f = Kxx.random_element(degree=n)
        for r in [1,5,20]:
            pp = [Kxx.random_element(degree=n-1) for i in range(r)]
            for s in [1,3,15]:
                qq = [Kxx.random_element(degree=n-1) for j in range(s)]
                for m in [1,5,30]:
                    R = SimultaneousTruncatedModularMultiplication(f,pp,qq,m)
                    if any([R[i,j] != ((pp[i]*qq[j]) % f).truncate(m) for i in range(r) for j in range(s)]):
                        print("SimultaneousTruncatedModularMultiplication Fail with",n,r,s,m)

    # TruncatedPowers
    for n in [10,50,200]:
        f = Kxx.random_element(degree=n)
        a = Kxx.random_element(degree=n-1)
        b = Kxx.random_element(degree=n-1)
        for m in [1,5,15,50]:
            for d in [1,15,50]:
                rr = TruncatedPowers(f,a,b,m,d)
                bakf = b
                for k in range(d):
                    if rr[k] != bakf.truncate(m):
                        print("TruncatedPowers Fail with",n,m,d)
                    bakf = (bakf*a) % f
        for m in [1,5,15,50]:
            #if m-1 < n: # requirement of algo: deg(b) < n = deg(f)
            b = xx**(m-1) % f
            for d in [1,15,50]:
                rr = TruncatedPowers(f,a,b,2*m-1,d)
                bakf = b
                for k in range(d):
                    if rr[k] != bakf.truncate(2*m-1):
                        print("Special TruncatedPowers Fail with",n,2*m-1,d)
                    bakf = (bakf*a) % f

    # BlockTruncatedPowers
    for n in [10,50,200]:
        f = Kxx.random_element(degree=n)
        while f.constant_coefficient() == 0:
            f = f+K.random_element()
        a = Kxx.random_element(degree=n-1)
        aa = [Kxx(1)]
        for k in range(1,50):
            aa.append((a*aa[k-1]) % f)
        for m in [1,5,15,50]:
            for d in [1,15,50]: # if more than 50 update aa above
                A = BlockTruncatedPowers(f,a,m,d)
                if any([A[i,k] != ((xx**i*aa[k]) % f).truncate(m) for i in range(m) for k in range(d)]):
                        print("BlockTruncatedPowers Fail with",n,m,d)

if MatRelTests:
    # CandidateBasis and BivariateModularCompositionWithRelationMatrix
    for n in [5,10,20]:
        f = Kxx.random_element(degree=n)
        while f(0) == 0:
            f = Kxx.random_element(degree=n)
        a = Kxx.random_element(degree=n-1)
        while f.gcd(a) != 1:
            a = Kxx.random_element(degree=n-1)
        for m in [n//5,n//2,n]:
            for d in [ceil(n/m),n//2,n]:
                if m*d >= n: # otherwise failure is basically guaranteed
                    fail = False
                    R,flag = CandidateBasis(f,a,m,d)
                    if not R.is_weak_popov(row_wise=False):
                        fail = True
                        print("CandidateBasis FAIL: not weak Popov; n,m,d=",n,m,d)
                    for j in range(m):
                        if sum([R[i,j](a)*xx**i % f for i in range(m)]) != 0:
                            fail = True
                            if flag == CERT:
                                print("CandidateBasis FAIL: column",j,"not relation; n,m,d=",n,m,d)
                            else:
                                print("CandidateBasis: flag NOCERT and column",j,"not relation; n,m,d=",n,m,d)
                    if not fail:
                        for dg in [4,9,30]:
                            g = Kxxyy.random_element(m+dg,terms=Infinity)
                            g = g.truncate(XX,m).truncate(YY,dg)
                            b = BivariateModularCompositionWithRelationMatrix(f,a,g,R)
                            if b != g(xx,a) % f:
                                print("BivariateModularCompositionWithRelationMatrix(univ g) FAIL; m,n,d,dg",m,n,d,dg)
                        for dg in [4,9,30,100]:
                            g = Kyy.random_element(degree=dg-1)
                            b = BivariateModularCompositionWithRelationMatrix(f,a,g,R)
                            if b != g(a) % f:
                                print("BivariateModularCompositionWithRelationMatrix(univ g) FAIL; m,n,d,dg",m,n,d,dg)

    # MatrixOfRelations
    # --> has been tested through degenerate examples such as those in file examples.sage

if ChBasTests:
    # ChangeOfBasis
    for n in [5,10,20]:
        f = Kxx.random_element(degree=n)
        while f(0) == 0:
            f = Kxx.random_element(degree=n)
        gamma = Kxx.random_element(degree=n-1)
        while f.gcd(gamma) != 1: # correctness clear when gcd != 1
            gamma = Kxx.random_element(degree=n-1)
        a = Kxx.random_element(degree=n-1)
        for m in [n//5,n//2,n]:
            for d in [ceil(n/m),n//2,n]:
                if m*d >= n: # otherwise failure is basically guaranteed
                    result = ChangeOfBasis(f,gamma,a,m,d)
                    if result == FAIL:
                        print("ChangeOfBasis possible FAIL: returned Fail on random input; n,m,d=",n,m,d)
                    else:
                        R,mu,alpha = result
                        if not R.is_popov(row_wise=False):
                            print("ChangeOfBasis FAIL: not Popov form; n,m,d=",n,m,d)
                        if sum(R.column_degrees()) != n:
                            print("ChangeOfBasis FAIL: R does not have degdet = n; n,m,d=",n,m,d)
                        for j in range(m):
                            if sum([R[i,j](gamma)*xx**i % f for i in range(m)]) != 0:
                                print("ChangeOfBasis FAIL: column",j,"not relation; n,m,d=",n,m,d)
                        if mu.degree() != n:
                            print("ChangeOfBasis FAIL: target minpoly has degree != n; n,m,d=",n,m,d)
                        if mu(gamma) % f != 0:
                            print("ChangeOfBasis FAIL: target minpoly is not annihilating; n,m,d=",n,m,d)
                        if alpha.degree() >= n:
                            print("ChangeOfBasis FAIL: target invmodcomp has degree >= n; n,m,d=",n,m,d)
                        if alpha(gamma) % f != a:
                            print("ChangeOfBasis FAIL: target invmodcomp does not realize inverse composition; n,m,d=",n,m,d)

if ModCompTests:
    for n in [5,10,20,50,200]:
        f = Kxx.random_element(degree=n)
        a = Kxx.random_element(degree=n-1)
        for deg_g in [n//2,n-1,n,n+1,2*n,5*n]:
            g = Kyy.random_element(degree=deg_g)
            rr = [K.random_element() for k in range(n+ceil(n**ETA))]
            b = ModularCompositionBaseCase(f,a,g,rr)
            if b == FAIL:
                print("ModularCompositionBaseCase returned FAIL on random input; n,m,deg_g=",n,deg_g)
            if g(a) % f != b:
                print("ModularCompositionBaseCase FAIL: not composition; n,m,deg_g=",n,deg_g)
