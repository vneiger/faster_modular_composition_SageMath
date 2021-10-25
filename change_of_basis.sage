###############################
#  Algorithms from Sections 6 #
###############################

# Note: the constant FAIL is defined in constants.sage,
# which should therefore be loaded if using the functions below

def ChangeOfBasis(f,gamma,a,m,d):
    """Change of basis

    :f: in K[x] of degree n, with f(0) != 0
    :gamma: in K[x] of degree < n
    :a: in K[x] of degree < n
    :m: positive integer at most n
    :d: positive integer
    :returns: either Fail or (Popov basis of relations, inverse modular composition, minpoly)

    """
    K = f.base_ring()
    Ky.<y> = K[]
    # Step 1
    if f.gcd(gamma) != 1:
        return FAIL
    # Step 2
    gammainv = gamma.inverse_mod(f)
    rr = TruncatedPowers(f, gammainv, (gammainv*a) % f, m, 2*d)
    s = sum([Matrix(K,m,1,rr[k].padded_list(m)) * y**k for k in range(2*d)])
    #s = Matrix(Ky,m,1)
    #for k in range(2*d):
    #    for i in range(m):
    #        s[i,0] += rr[k][i] * y**k
    # Initialize matrix F of step 4
    F = Matrix.block(Ky,[[Matrix.zero(Ky,m,m),-Matrix.identity(Ky,m),s]])
    # Step 3
    Gamma = BlockTruncatedPowers(f, gammainv, m, 2*(d+1))
    for k in range(2*d):
        for i in range(m):
            for j in range(m):
                F[j,i] += -Gamma[i,k+1][j] * y**k
    # Step 4
    t = [0]*(2*m) + [2*d]
    PP = F.minimal_approximant_basis(2*d,shifts=t,row_wise=False,normal_form=True)
    P = PP[:2*m,:2*m]
    vaa = PP[:m,2*m]
    # Step 5
    R = P[:m,:m]
    if sum([R[i,i].degree() for i in range(m)]) < n:
        return FAIL
    if min(P[:,m:].column_degrees()) < R.degree():
        return FAIL
    # Step 6
    # using antitranposes because
    # - Sage's HNF is row-wise (row transformation) while we need column-wise
    # - Sage's HNF is upper triangular (which becomes lower triangular
    # if we just use normal transposes)
    T = R.antitranspose().hermite_form().antitranspose()
    mu = T[0,0]
    if mu.degree() != n:
        return FAIL
    # Step 7
    alpha = vaa[0,0] - sum([T[0,j]*vaa[j,0] for j in range(1,m)])
    alpha = alpha % mu
    # Return
    return (R,mu,alpha)
