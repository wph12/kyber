import numpy as np
from numpy.polynomial.polynomial import Polynomial

from polynomials import add_poly, mul_poly_simple, sign_extend

def add_vec(v0, v1, q):
    """
    adds 2 vectors of polynomials in polynomial ring R := Z_q[x]/f(x)
    params: 
        v0(int[k][n]):  vector representing a polynomial. easier to think of this as a list of polynomial[k]
        v1(int[k][n]):  vector representing a polynomial
        q(int): plain modulo for operations in Z_q
    returns:
        result(int[k][n]): vector representing a polynomial
    """
    assert(len(v0) == len(v1)) # sizes need to be the same

    result = []

    for i in range(len(v0)):
        result.append(add_poly(v0[i], v1[i], q))
    
    return result


def mul_vec_simple(v0, v1, f, q):
    """
    computes dot product of 2 vectors in polynomial ring R := Z_q[x]/f(x)
    params: 
        v0(int[k][n]):  vector representing a polynomial. easier to think of this as a list of polynomial[k]
        v1(int[k][n]):  vector representing a polynomial
        f (int[n+1]): reduction polynomial's coefficients
        q(int): plain modulo for operations in Z_q
    returns:
        result(int[n]): vector representing a polynomial
    """
    assert(len(v0) == len(v1)) # sizes need to be the same

    degree_f = len(f) - 1
    result = [0 for i in range(degree_f - 1)]

    # textbook vector inner product
    for i in range(len(v0)):
        result = add_poly(result, mul_poly_simple(v0[i], v1[i], f, q), q)
    
    return result


def mul_mat_vec_simple(M, a, f, q):
    """
    computes matrix vector multiplication Ma in polynomial ring R := Z_q[x]/f(x)
    
    params: 
        M(int[k][k][n]):  matrix of polynomials. easier to think of this as a matrix of polynomial[k][k]
        a(int[k][n]):  vector representing a polynomial. easier to think of this as a list of polynomial[k]
        f (int[n+1]): reduction polynomial's coefficients
        q(int): plain modulo for operations in Z_q
    returns:
        result(int[k][n]): vector representing a polynomial
    """
    result = []
    
    # textbook matrix-vector multiplication
    for i in range(len(M)):
        result.append(mul_vec_simple(M[i], a, f, q))
    
    return result


def transpose(m):
    result = [[None for i in range(len(m))] for j in range(len(m[0]))]

    for i in range(len(m)):
        for j in range(len(m[0])):
            result[j][i] = m[i][j]
    
    return result

#####################################################################################
# TESTING CODE HERE
np.random.seed(0xdeadbeef)

def test_mul_vec(N, k, f, q):
  degree_f = len(f) - 1

  for i in range(N):
    m = (np.random.random([k, k, degree_f]) * q).astype(int)
    a = (np.random.random([k, degree_f]) * q).astype(int)

    m_mul_a = mul_mat_vec_simple(m, a, f, q)

    m_poly = list(map(lambda x: list(map(Polynomial, x)), m))
    a_poly = list(map(Polynomial, a))
    prod = np.dot(m_poly, a_poly)
    m_mul_a_ref = list(map(lambda x: list(map(lambda y: int(y) % q, sign_extend((x % Polynomial(f)).coef, degree_f))), prod))

    assert(m_mul_a == m_mul_a_ref)

def main():
    test_mul_vec(100, 2, [1, 0, 0, 0, 1], 17)

if __name__ == "__main__":
    main()
