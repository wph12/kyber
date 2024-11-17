import numpy as np
from numpy.polynomial.polynomial import Polynomial

def add_poly(a, b, q):
    """
    adds two polynomials modulo q
    params:
        a (list[int]): list of integers representing polynomial a's coefficients with a[i] being the coefficient of x^i
        b (list[int]): list of integers representing polynomial b's coefficients
        q (int): plain modulo for operations in Z_q
    returns:
        result (list[int]): list of integers representing polynomial a + b (mod q)
    """
    result = [0] * max(len(a), len(b))
    for i in range(max(len(a), len(b))):
        if i < len(a):
            result[i] += a[i]
        if i < len(b):
            result[i] += b[i]
            result[i] %= q
    return result


def inv_poly(a, q):
    """
    makes all elements of polynomial negative
    params:
        a (list[int]): list of integers representing polynomial coefficients with a[i] being the coefficient of x^i
        q (int): plain modulo for polynomial operations
    returns:
        result (list[int]): list of integers representing -a 
    """
    return list(map(lambda x: -x % q, a))


def sub_poly(a, b, q):
    """
    gives a - b (mod q)
    params:
        a (list[int]): list of integers representing polynomial a's coefficients with a[i] being the coefficient of x^i
        b (list[int]): list of integers representing polynomial b's coefficients
        q (int): plain modulo for polynomial operation
    returns:
        result (list[int]): list of integers representing polynomial a - b (mod q)
    """
    return add_poly(a, inv_poly(b, q), q)


def mul_poly_simple(a, b, f, q):
    """
    computes multiplication of polynomials a and b over polynomial ring R defined by 
    R := Z_q[x]/f(x)   

    params:
        a (list[int]): list of integers representing polynomial a's coefficients with a[i] being the coefficient of x^ia 
        b (list[int]): list of integers representing polynomial b's coefficients
        f (list[int]): reduction polynomial's coefficients
        q (int): plain modulo for polynomial operation
    returns:
        result (list[int]): list of integers representing resulting polynomial a * b (mod q)
    """
    tmp = [0] * (len(a) * 2 - 1) # the product of two degree n polynomial cannot exceed 2n
    
    # schoolbook multiplication
    for i in range(len(a)):
        # perform a_i * b
        for j in range(len(b)):
            tmp[i + j] += a[i] * b[j]
    
    # take polynomial modulo f
    # since Kyber's poly modulus is always x^n + 1,
    # we can efficiently compute the remainder
    degree_f = len(f) - 1
    for i in range(degree_f, len(tmp)):
        tmp[i - degree_f] -= tmp[i]
        tmp[i] = 0

    # take coefficients modulo q
    tmp = list(map(lambda x: x % q, tmp))
    return tmp[:degree_f]


#####################################################################################
#TESTING 
def sign_extend(poly, degree):
    if len(poly) >= degree:
        return poly
    return [0] * (degree - len(poly))

def test_mul_poly(N, f, q):
    degree_f = len(f) - 1

    for i in range(N):
        a = (np.random.random(degree_f) * q).astype(int)
        b = (np.random.random(degree_f) * q).astype(int)
        
        a_mul_b = mul_poly_simple(a, b, f, q)
        
        # NumPy reference poly multiplication
        # note that we need to convert the coefficients to int and extend the list to match the fixed size of our impl
        a_mul_b_ref = list(map(lambda x: int(x) % q, ((Polynomial(a) * Polynomial(b)) % Polynomial(f)).coef))
        a_mul_b_ref = sign_extend(a_mul_b_ref, degree_f)

        assert(a_mul_b == a_mul_b_ref)

np.random.seed(0x71904)
test_mul_poly(100, [1, 0, 0, 0, 1], 17)
