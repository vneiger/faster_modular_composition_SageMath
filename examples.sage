reset()
load('constants.sage')
load('simultaneous_modular_operations.sage')
load('matrix_of_relations.sage')
load('change_of_basis.sage')
load('modular_composition_basecase.sage')

################################
#  Examples and illustrations  #
################################


# There are four examples: RandomSmall, Random, SingularSmall, Singular.
#  --> "Small" examples have small n, and are very verbose (they print out
#  matrices of relations and polynomials completely)
#  --> "NonSmall" examples have slightly larger n, and are a bit less verbose
#  (they print only degrees in matrices of relations, and they print completely
#  some relevant polynomials)
#  --> "Random" examples have input polynomials taken at random from a finite field
#  of size 997, making them these examples "generic" (for these one could
#  instead only rely on the algorithm MatrixOfRelations, without going through a
#  random change of basis).
#  --> "Singular" examples have input polynomials chosen on purpose so that the
#  MatrixOfRelations algorithm would fail, and therefore for these it is crucial
#  to go through the change of basis.

#  To activate an example, set the corresponding flag to True; to disable, set to False.
ExRandomSmall=True
ExRandom=True
ExSingularSmall=True
ExSingular=True

# Here the base field can be changed. Note that setting a too small field
# may lead to frequent failures in the algorithms.
K = GF(997)
Kx.<x> = K[]
Ky.<y> = K[]

def TryMatrixOfRelations(f,a):
    """Trying matrix of relations

    :f: in K[x] of degree n with f(0) != 0
    :a: in K[x] of degree < n with gcd(a,f) = 1

    Function for trying to compute a matrix of relations directly and observing
    whether this could have been sufficient to perform modular composition,
    without the need to resort to a change of basis.

    """
    print("====================================")
    print("====  Trying MatrixOfRelations  ====")
    print("====================================")
    print("========")
    print("Computing matrix of relations directly, without change of basis")
    print("--> This is done for experimental purpose, to observe whether this does or "   \
            "does not provide a matrix that can be used for composition.")
    n = f.degree()
    m = ceil(n**ETA)
    d = ceil(n/m)
    print("========")
    print("Chosen parameters: m =",m,"and d =",d)
    if f.gcd(a) != 1:
        print("TryMatrixOfRelations: requirement gcd(f,a) = 1 not satisfied, aborting.")
        return
    if f(0) == 0:
        print("TryMatrixOfRelations: requirement f(0) != 0 not satisfied, aborting.")
        return
    print("========")
    print("Entering CandidateBasis...")
    R,Flag = CandidateBasis(f,a,m,d)
    print("Returned flag:",Flag)
    print("Degree matrix of the returned matrix:")
    print(R.degree_matrix())
    if Flag == CERT:
        print("==> for this input we could have relied on CandidateBasis, " \
                "without resorting to ChangeOfBasis, since it returns a " \
                "certified basis of relations.")
        # if Flag == CERT, the matrix is certified to be a basis of relations;
        # still, one can change DEBUG (from `constants.sage`) from "False" into
        # "True" to check that it is indeed the case; this can be useful for
        # tracking possible bugs
        if DEBUG:
            print("Verification of the result...")
            correct = VerifyMatrixOfRelations(f,a,R,verbose=1)
            if not correct:
                print("--> verification failed; there is a bug in MatrixOfRelations")
            else:
                print("--> verification succeeded; all is well")
    if Flag == NOCERT:
        print("Flag NoCert, verification of the result...")
        correct = VerifyMatrixOfRelations(f,a,R,verbose=1)
        if correct:
            print("--> this verification succeeded")
            print("==> this verification shows that for this input we could have relied on " \
                    "CandidateBasis, without resorting to ChangeOfBasis; yet it is not known " \
                    "how to perform this verification in good complexity.")
            print("We therefore call MatrixOfRelations, expecting to find " \
                    "a certified matrix of relations.")
        else:
            print("--> this verification failed")
            print("==> for this input we cannot rely on CandidateBasis " \
                "(without resorting to ChangeOfBasis) to perform modular composition; "  \
                "let us now see if MatrixOfRelations could have sufficed.")
        print("========")
        print("Entering MatrixOfRelations...")
        rr = [f.base_ring().random_element() for i in range(m-2)]
        R = MatrixOfRelations(f,a,m,d,rr)
        if R == FAIL:
            print("--> this call returned Fail.")
            print("==> for this input one does not get a matrix of relations from " \
                    "either CandidateBasis or MatrixOfRelations; " \
                    "on the other hand, going through ChangeOfBasis will still "  \
                    "allow us to perform modular composition.")
        else:
            print("--> returned a matrix (and not Fail), so this matrix is certified " \
                    "to be a matrix of relations.")
            print("==> for this input we could have relied on MatrixOfRelations, " \
                "without resorting to ChangeOfBasis.")
            # change the variable DEBUG (from `constants.sage`) from "False"
            # into "True" to still check that the matrix is indeed a matrix of
            # relations; this can be useful for tracking possible bugs
            if DEBUG:
                print("verifying result of MatrixOfRelations...")
                correct = VerifyMatrixOfRelations(f,a,R,verbose=1)
                if not correct:
                    print("--> verification failed; there is a bug in MatrixOfRelations")
                else:
                    print("--> verification succeeded; all is well")

if ExRandomSmall:
    n = 10
    f = Kx.random_element(degree=n)
    a = Kx.random_element(degree=n-1)
    print("\n\n")
    print("========")
    print("RandomSmall example: random polynomials")
    print(". over",K)
    print(". with f(x) =",f)
    print(". with a(x) =",a)
    TryMatrixOfRelations(f,a)
    g = Ky.random_element(degree=n)
    rr = [K.random_element() for i in range(n+ceil(n**ETA))]
    b = ModularCompositionBaseCaseVerbose(f,a,g,rr,verbose=2)
    print("Testing result...")
    if b != FAIL and b == g(a) % f :
        print("SUCCESS. The algorithm has not failed and has returned the correct polynomial,\n\tb = g(a) rem f =",b)

if ExRandom:
    n = 700
    f = Kx.random_element(degree=n)
    a = Kx.random_element(degree=n-1)
    print("\n\n")
    print("========")
    print("Random example: random polynomials")
    print(". over",K)
    print(". with f(x) of degree",f.degree())
    print(". with a(x) of degree",a.degree())
    TryMatrixOfRelations(f,a)
    g = Ky.random_element(degree=n)
    rr = [K.random_element() for i in range(n+ceil(n**ETA))]
    b = ModularCompositionBaseCaseVerbose(f,a,g,rr,verbose=1)
    print("========")
    print("Testing result...")
    if b != FAIL and b == g(a) % f :
        print("SUCCESS. The algorithm has not failed and has returned the correct polynomial.")

if ExSingularSmall:
    n = 12
    f = x**n # power series case
    a = 1 + x + x**(n//2)
    print("\n\n")
    print("========")
    print("SingularSmall example, power series case with")
    print("f(x) =",f)
    print("a(x) =",a)
    TryMatrixOfRelations(f(x-1),a(x-1))  # shift x-1 to ensure f(0) != 0
    g = Ky.random_element(degree=n)
    rr = [K.random_element() for i in range(n+ceil(n**ETA))]
    b = ModularCompositionBaseCaseVerbose(f,a,g,rr,verbose=2)
    print("Testing result...")
    if b != FAIL and b == g(a) % f :
        print("SUCCESS. The algorithm has not failed and has returned the correct polynomial,\n\tb = g(a) rem f =",b)

if ExSingular:
    n = 700
    f = x**n # power series case
    a = 1 + x + x**(n//2)
    print("\n\n")
    print("========")
    print("Singular example, power series case with")
    print("f(x) =",f)
    print("a(x) =",a)
    TryMatrixOfRelations(f(x-1),a(x-1))  # shift x-1 to ensure f(0) != 0
    g = Ky.random_element(degree=n)
    rr = [K.random_element() for i in range(n+ceil(n**ETA))]
    b = ModularCompositionBaseCaseVerbose(f,a,g,rr,verbose=1)
    ("========")
    print("Testing result...")
    if b != FAIL and b == g(a) % f :
        print("SUCCESS. The algorithm has not failed and has returned the correct polynomial.")
