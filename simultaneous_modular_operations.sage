###############################
#  Algorithms from Section 3  #
###############################

# Note: the constant FAIL is defined in constants.sage,
# which should therefore be loaded if using the functions below

# The following algorithm `ModularPower` is added here because in SageMath the
# direct writing "a**k % f" does the remainder only after computing the "full"
# power a**k (which can be substantially slower if k is not small)
# Input: f nonzero and a both in K[x], k nonnegative integer
# Output: a**k % f
def ModularPower(f,a,k):
    return f.parent().quotient(f)(a**k).lift()

def ModularComposition_BrentKung(f,a,g,d):
    """Brent and Kung's algorithm [1978]

    :f: in K[x], of degree n
    :a: in K[x], of degree < n
    :g: in K[y], of degree < d
    :d: positive integer
    :returns: g(a) rem f

    """
    K = f.base_ring()
    Kx = f.parent()
    n = f.degree()
    r = ceil(sqrt(d))
    s = ceil(d/r)
    aa = [Kx(1)]
    for i in range(1,r+1):
        aa.append((a*aa[i-1]) % f)  # aa[i] = a**i rem f
    A = Matrix(K, r, n)
    for i in range(r):
        for j in range(n):
            A[i,j] = aa[i][j] # coeff of degree j of aa[i]
    G = Matrix(K, s, r)
    for i in range(s):
        for j in range(r):
            G[i,j] = g[i*r+j]
    B = G*A
    b = Kx(0)
    for i in range(s): # Horner evaluation
        b = Kx(B[s-i-1,:].list()) + ((b*aa[r]) % f)
    return b

def PowerProjection(f,a,d,rr):
    """Power projection, naive implementation

    :f: in K[x], of degree n
    :a: in K[x], of degree < n
    :d: positive integer
    :rr: list of n elements of K
    :returns: ell(1),ell(a),...,ell(a**(d-1) mod f)
    with ell given by the vector rr
    """
    Kx = f.parent()
    n = f.degree()
    apow = Kx(1)
    pp = []
    for i in range(d):
        pp.append(sum([rr[i]*apow[i] for i in range(n)]))
        apow = (apow * a) % f
    return pp

def ModularComposition_SmallMinimalPolynomial(f,a,g,d,rr):
    """Modular composition when minpoly has small degree

    :f: in K[x], of degree n
    :a: in K[x], of degree < n
    :g: in K[y], of degree < n
    :d: positive integer in {1...n}
    :returns: g(a) rem f or Fail

    """
    Ky = g.parent()
    vv = PowerProjection(f,a,2*d,rr)
    # minpoly via rational reconstruction
    F = Matrix([[Ky(list(reversed(vv)))],[-1]])
    P = F.minimal_approximant_basis(2*d)
    mu = P[0,0]
    t = ModularComposition_BrentKung(f,a,mu,n)
    if t != 0:
        return FAIL
    b = ModularComposition_BrentKung(f,a,g % mu,mu.degree())
    return b

def SimultaneousBivariateModularComposition(f,a,gg,m,r):
    """SimultaneousBivariateModularComposition

    :f: in K[x], of degree n
    :a: in K[x], of degree < n
    :gg: list of elements in K[x,y]_{<(m,r)}
    :m: positive integer
    :r: positive integer
    :returns: (g[i](x,a) % f for i = 0 ... s-1)

    """
    Kx = f.parent()
    x = Kx.gen()
    (X,Y) = gg[0].parent().gens()
    s = len(gg)
    aa = [Kx(1)]
    for i in range(1,r):
        aa.append((a*aa[i-1]) % f)
    A = Matrix(Kx,r,ceil(n/m))
    for i in range(r):
        for j in range(ceil(n/m)):
            A[i,j] = aa[i].truncate(m)
            aa[i] = aa[i].shift(-m)
    G = Matrix(Kx,s,r)
    for i in range(s):
        for j in range(r):
            G[i,j] = (gg[i].coefficient(Y**j))(x,0)
    B = G*A
    bb = []
    for i in range(s):
        bb.append(sum([B[i,j].shift(j*m) for j in range(ceil(n/m))]) % f)
    return bb

def BivariateModularComposition(f,a,g,m,d):
    """Nusken and Ziegler's algorithm for bivariate modular composition

    :f: in K[x], of degree n
    :a: in K[x], of degree < n
    :g: in K[x,y], of degree < m,d
    :m: positive integer
    :d: positive integer
    :returns: g(x,a) rem f

    """
    Kx = f.parent()
    Kxy = g.parent()
    X,Y = Kxy.gens()
    r = ceil(sqrt(d))
    s = ceil(d/r)
    gg = []
    for i in range(s):
        gtrunc = g.truncate(Y,r)
        gg.append(gtrunc)
        g = (g - gtrunc) // Y^r
    bb = SimultaneousBivariateModularComposition(f,a,gg,m,r)
    apow = ModularPower(f,a,r)
    b = Kx(0)
    for i in range(s): # Horner evaluation
        b = bb[s-i-1] + ((b*apow) % f)
    return b

    
def SimultaneousTruncatedModularMultiplication(f,pp,qq,m):
    """Simultaneous Truncated Modular Multiplication

    :f: in K[x] of degree n
    :pp: r elements in K[x] of degree < n
    :qq: s elements in K[x] of degree < n
    :m: positive integer
    :returns: r x s matrix with entry i,j equal to (pp[i]*qq[j] % f).truncate(m)

    """
    Kx = f.parent()
    r = len(pp)
    s = len(qq)
    (ell,t) = (n-m-1).quo_rem(m)
    if m >= n:
        (ell,t) = (0,n-m-1)
    pprev = [p.reverse(n-1) for p in pp]
    inv_rev_f = f.reverse(n).inverse_series_trunc(n-1)
    qqexp = [(q.reverse(n-1) * inv_rev_f).truncate(n-1) for q in qq]
    P1 = Matrix(Kx,r,ell+1)
    for i in range(r):
        pi = pprev[i].shift(-t)
        for j in range(ell+1):
            P1[i,j] = pi.truncate(m)
            pi = pi.shift(-m)
    P2 = P1[:,:ell]  # submatrix of all rows and first ell columns
    Q1 = Matrix(Kx,ell+1,s)
    for j in range(s):
        qj = qqexp[j]
        for i in range(ell+1):
            Q1[ell-i,j] = qj.truncate(m)
            qj = qj.shift(-m)
    Q2 = Q1[1:,:]  # submatrix of all columns and last ell rows
    C = Matrix(Kx,r,s)  # correcting term
    cpi = [prev.truncate(t) for prev in pprev]
    cqj = [qexp.shift(-ell*m-1).truncate(m-1+t) for qexp in qqexp]
    for i in range(r):
        for j in range(s):
            C[i,j] = (cpi[i]*cqj[j]).shift(-t+1).truncate(m)
    H = (P1*Q1).truncate(m) + (P2*Q2).shift(-m).truncate(m) + C
    R = Matrix(Kx,r,s)  # output matrix
    for i in range(r):
        for j in range(s):
            R[i,j] = (pp[i]*qq[j] - H[i,j].reverse(m-1)*f).truncate(m)
    return R

def TruncatedPowers(f,a,b,m,d):
    """Computes truncations of powers of the form b*a**k % f

    :f: in K[x] of degree n
    :a: in K[x] of degree <n
    :b: in K[x] of degree <n
    :m: positive integer
    :d: positive integer
    :returns: list of (b*a**k % f).truncate(m) for 0 <= k < d

    """
    Kx = f.parent()
    r = ceil(sqrt(d))
    s = ceil(d/r)
    aa = [Kx(1)]
    for i in range(1,r+1):
        aa.append((a*aa[i-1]) % f)
    aaa = [b]
    for j in range(1,s):
        aaa.append((aa[r]*aaa[j-1]) % f)
    C = SimultaneousTruncatedModularMultiplication(f,aa,aaa,m)
    rr = []
    for j in range(s-1):
        for i in range(r):
            rr.append(C[i,j])
    for i in range(d-(s-1)*r):
        rr.append(C[i,s-1])
    return rr

def BlockTruncatedPowers(f,a,m,d):
    """Computes truncations of powers of the form x**i a**k modulo f

    :f: in K[x] of degree n
    :a: in K[x] of degree <n
    :m: positive integer
    :d: positive integer
    :returns: matrix whose entry (i,j) is (x**i * a**k % f).truncate(m) for 0<=i<m, 0<=k<d

    """
    Kx = f.parent()
    x = Kx.gen()
    rr = TruncatedPowers(f,a,x**(m-1) % f,2*m-1,d)
    f0_inv = - 1 / f.constant_coefficient()
    A = Matrix(Kx, m, d)
    for k in range(d):
        A[m-1,k] = rr[k]
        for i in range(m-1,0,-1): # i from m-1 down to 1
            c = A[i,k].constant_coefficient() * f0_inv
            A[i-1,k] = (A[i,k] + c * f.truncate(m+i)).shift(-1)
    A = A.truncate(m)
    return A
