######################################
#  Algorithms from Sections 4 and 5  #
######################################

# Note: the constants FAIL, CERT, and NOCERT are defined in constants.sage,
# which should therefore be loaded if using the functions below

def BivariateModularCompositionWithRelationMatrix(f,a,g,R):
    """Bivariate composition when a relation matrix is known

    :f: in K[x] of degree n
    :a: in K[x] of degree <n
    :g: in K[x,y] of degree <(m,.)
    :R: in K[y]^{mxm}, matrix of relations of a modulo f
    :returns: g(x,a) rem f

    """
    K = f.base_ring()
    Ky.<y> = K[]
    m = R.nrows()
    d = R.degree()
    if g.parent().ngens() == 1: # g univariate in y
        Kxy.<X,Y> = K[]
        dg = g.degree()
        vg = Matrix(Ky, m, 1, [g] + [0]*(m-1))
    else:
        X,Y = g.parent().gens()
        dg = g.degree(X)
        vg = Matrix(Ky, m, 1, [g.coefficient(X**k)(0,y) for k in range(m)])
    K = Matrix.block([[R,-vg]]).minimal_kernel_basis(row_wise=False,shifts=[d]*m+[dg])
    vgg = R * (K[:m,0] % K[m,0])
    vgg = [vgg[k,0] // K[m,0] for k in range(m)]
    gg = sum([vgg[k](Y)*X**k for k in range(m)])
    return BivariateModularComposition(f,a,gg,m,d)

def CandidateBasis(f,a,m,d):
    """Candidate matrix of relations

    :f: in K[x] of degree n with f(0) != 0
    :a: in K[x] of degree < n with gcd(a,f) = 1
    :m: positive integer, at most n
    :d: positive integer
    :returns: (R,CERT) or (R,NOCERT), see paper for more details

    """
    Ky.<y> = f.base_ring()[]
    A = BlockTruncatedPowers(f,a.inverse_mod(f),m,2*(d+1))
    F = Matrix.block(Ky,[[Matrix.zero(Ky,m,m),-Matrix.identity(Ky,m)]])
    for k in range(2*d):
        for i in range(m):
            for j in range(m):
                F[j,i] += -A[i,k+1][j] * y**k
    P = F.minimal_approximant_basis(2*d,row_wise=False)
    R = P[:m,:m]
    if sum([R[i,i].degree() for i in range(m)]) == n \
            and min(P[:,m:].column_degrees()) >= R.degree():
        return (R, CERT)
    else:
        return (R, NOCERT)

def MatrixOfRelations(f,a,m,d,rr):
    """Computes a certified matrix of relations via Sylvester matrix

    :f: in K[x] of degree n with f(0) != 0
    :a: in K[x] of degree < n with gcd(a,f) = 1
    :m: positive integer, at most n
    :d: positive integer
    :rr: list of m-2 elements of K
    :returns: either Fail or a matrix of relations

    """
    Kxy.<X,Y> = f.base_ring()[]
    R,Flag = CandidateBasis(f,a,m,d)
    if Flag == CERT:
        return R
    if m == 1:
        if ModularComposition_BrentKung(f,a,R[0,0]) != 0:
            return FAIL
        else:
            return R
    # arriving here, Flag == NOCERT and m > 1
    dR = R.degree()
    r = sum([R[i,0](Y) * X**i for i in range(m)])
    comb = R[:,1] + sum([rr[i]*R[:,i+2] for i in range(m-2)])
    s = sum([comb[i,0](Y) * X**i for i in range(m)])
    if BivariateModularComposition(f,a,r,m,dR+1) != 0 \
            or BivariateModularComposition(f,a,s,m,dR+1) != 0:
        return FAIL
    if m == 2:
        return R
    if r.gcd(s) != 1:
        return FAIL
    #dr = r.degree(X)
    #ds = s.degree(X)
    #S = Matrix(R.base_ring(), dr+ds, dr+ds)
    #for j in range(ds):
    #    for i in range(dr):
    #        S[j+i,j] = R[i,0]
    #for j in range(dr):
    #    for i in range(ds):
    #        S[j+i,ds+j] = comb[i,0]
    S = r.sylvester_matrix(s,X).antitranspose()
    return S

#####################################
#  Verifying a matrix of relations  #
#####################################
# The following algorithm is NOT in the paper:
# we do not know how to perform this verification in sufficiently good complexity.
# It is only used here for testing and for experimental purposes, see for
# example the file `examples.sage`.

def VerifyMatrixOfRelations(f,a,R,verbose=0):
    """Verifying whether the matrix R is a matrix of relations

    :f: in K[x] of degree n with f(0) != 0
    :a: in K[x] of degree < n with gcd(a,f) = 1
    :R: matrix in K[y]^(m x m)
    :returns: boolean, True iff R is a matrix of relations

    """
    m = R.nrows()
    d = R.degree()
    if R.ncols() != m:
        if verbose>0:
            print("VerifyMatrixOfRelations Fail: input matrix is not square")
        return False
    if R.determinant() == 0:
        if verbose>0:
            print("VerifyMatrixOfRelations Fail: input matrix is not nonsingular")
        return False
    # we are going to evaluate the whole matrix R at y=a modulo f
    # for this, we first transform R into a polynomial with matrix coefficients
    matPol.<Y> = PolynomialRing(MatrixSpace(f.base_ring(), m, m))
    # then we do Horner evaluation
    # first take coefficient of degree d = deg(R), seen as matrix over K[x]
    E = Matrix(f.parent(), m, m, R.coefficient_matrix(d))
    for k in range(d-1,-1,-1): # k from d-1 down to 0 (included)
        E = ((a * E) + R.coefficient_matrix(k)) % f
    # check that each column represents a polynomial 0 mod f
    for j in range(m):
        if sum([E[i,j]*x**i for i in range(m)]) % f != 0:
            if verbose>0:
                print("VerifyMatrixOfRelations Fail: column",j,"not relation")
            return False
    return True

