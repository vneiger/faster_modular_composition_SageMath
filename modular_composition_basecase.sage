##############################
#  Algorithm from Section 8  #
##############################

# Note: there are two versions below;
# - ModularCompositionBaseCase is the version from Section 8 of the paper
# - ModularCompositionBaseCaseVerbose is exactly the same computations,
#   plus additional printing statements for illustration purposes
#   (this is used to make the examples verbose in `examples.sage`)

# Note: the constants FAIL and ETA are defined in constants.sage,
# which should therefore be loaded if using the functions below


def ModularCompositionBaseCase(f,a,g,rr):
    """Univariate modular composition g(a) rem f

    :f: in K[x] of degree n
    :a: in K[x] of degree < n
    :g: in K[y]
    :rr: list of n+ceil(n^ETA) in K
    :returns: the modular composition g(a) % f, or Fail

    """
    Kx = f.parent()
    x = Kx.gen()
    Ky = g.parent()
    y = Ky.gen()
    # Steps 1-4
    if n==1:
        return g(a)
    g = g(y-rr[0])
    a = a + rr[0]
    if f.gcd(a) != 1:
        return FAIL
    f = f(x+rr[1])
    a = a(x+rr[1])
    if f.constant_coefficient() == 0:
        return FAIL
    m = ceil(n**ETA)
    d = ceil(n/m)
    # Steps 5-6
    gamma = sum([rr[k+2]*x**k for k in range(n)])
    cob = ChangeOfBasis(f,gamma,a,m,d)
    if cob == FAIL:
        return FAIL
    R_gamma,mu,alpha = cob
    if mu.constant_coefficient() == 0:
        return FAIL
    # Step 8
    R_alpha = MatrixOfRelations(mu,alpha,m,d,rr[n+2:])
    if R_alpha == FAIL:
        return FAIL
    # Step 9
    beta = BivariateModularCompositionWithRelationMatrix(mu,alpha,g,R_alpha)
    # Step 11
    b = BivariateModularCompositionWithRelationMatrix(f,gamma,beta,R_gamma)
    return b(x-rr[1])


    

def ModularCompositionBaseCaseVerbose(f,a,g,rr,verbose=0):
    """Univariate modular composition g(a) rem f

    :f: in K[x] of degree n
    :a: in K[x] of degree < n
    :g: in K[y]
    :rr: list of n+ceil(n^ETA) in K
    :verbose: level of verbosity, 0 (none) or 1 (medium) or 2 (high)
    :returns: the modular composition g(a) % f, or Fail

    """
    if verbose > 0:
        print("===============================================")
        print("====  Entering ModularCompositionBaseCase  ====")
        print("===============================================")
        print("Parameters n , m , d =",f.degree(),",",ceil(n**ETA),",",ceil(n/ceil(n**ETA)))
    if verbose == 2:
        print("Input polynomials f,a,g:\n\tf =",f,"\n\ta =",a,"\n\tg =",g)
    if verbose > 0:
        print("Initialization and polynomial shifts...")
    Kx = f.parent()
    x = Kx.gen()
    Ky = g.parent()
    y = Ky.gen()
    # Steps 1-4
    if n==1:
        return g(a)
    g = g(y-rr[0])
    a = a + rr[0]
    if f.gcd(a) != 1:
        return FAIL
    f = f(x+rr[1])
    a = a(x+rr[1])
    if f.constant_coefficient() == 0:
        return FAIL
    m = ceil(n**ETA)
    d = ceil(n/m)
    # Steps 5-6
    if verbose > 0:
        print("========")
        print("Change of basis...")
    gamma = sum([rr[k+2]*x**k for k in range(n)])
    cob = ChangeOfBasis(f,gamma,a,m,d)
    if cob == FAIL:
        return FAIL
    R_gamma,mu,alpha = cob
    if mu.constant_coefficient() == 0:
        return FAIL
    if verbose == 1:
        print("--> found matrix of relations R_gamma, its degree matrix is")
        print(R_gamma.degree_matrix())
        print("--> found minimal polynomial mu_gamma")
        print("--> found inverse composition alpha")
    if verbose == 2:
        print("--> found matrix of relations R_gamma:")
        print(R_gamma)
        print("--> found minimal polynomial mu_gamma:")
        print("\t",mu)
        print("--> found inverse composition alpha:")
        print("\t",alpha)
    # Step 8
    if verbose > 0:
        print("========")
        print("Matrix of relations for alpha,mu_gamma...")
    R_alpha = MatrixOfRelations(mu,alpha,m,d,rr[n+2:])
    if R_alpha == FAIL:
        return FAIL
    if verbose == 1:
        print("--> found matrix of relations R_alpha, its degree matrix is")
        print(R_alpha.degree_matrix())
    if verbose == 2:
        print("--> found matrix of relations R_alpha:")
        print(R_alpha)
    # Step 9
    if verbose > 0:
        print("========")
        print("Composition beta = g(alpha) rem mu_gamma...")
    beta = BivariateModularCompositionWithRelationMatrix(mu,alpha,g,R_alpha)
    if verbose == 1:
        print("--> done")
    if verbose == 2:
        print("--> beta =", beta)
    # Step 11
    if verbose > 0:
        print("========")
        print("Inverse change of basis b = beta(gamma) rem f...")
    b = BivariateModularCompositionWithRelationMatrix(f,gamma,beta,R_gamma)
    if verbose == 1:
        print("--> done, now returning.")
    if verbose == 2:
        print("--> b =", b)
        print("========")
        print("==> now returning result, b(x-rr[1]) =",b(x-rr[1]))
        print("========")
    return b(x-rr[1])
