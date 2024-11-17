from polynomials import *
from module import * 
import math

def naive_sample_binom(n ,k, eta):
    """
    samples array of k integers from central binomial distribution
    note: numpy is not cryptographically secure

    params:
        n (int): number of elements in a polynomial. given by deg(f)
        k (int): length of array to sample
        eta (int): maximum absolute polynomial coefficient
    returns:
        result (int[k][n]): polynomial vector with coefficients of range [-eta, eta]
    """
    depth = 2 * eta
    binom_coefficients = np.array([math.comb(depth, k) for k in range(depth + 1)])
    choice_arr = [0]
    for i in range(1, eta+1):
        choice_arr.insert(0, -i)
        choice_arr.append(i)    
    
    probs_arr = binom_coefficients / np.sum(binom_coefficients)
    result = []
    for i in range(k):
        result.append(np.random.choice(choice_arr, size=n, p= probs_arr).tolist())
    
    return result


def naive_keygen(q, f, n, k, eta_1=2, eta_2=2):
    """
    naive key generation function for kyber-PKE. operations are done in polynomial ring R := Z_q[x]/f(x)
    note: numpy is not cryptographically secure

    params:
        q (int): plain modulo for operations in Z_q
        f (int[n]): reduction polynomial's coefficients
        n (int): number of elements in a polynomial. given by deg(f)
        k (int): number of polynomials in a polynomial vector
        eta_1 (int): maximum absolute polynomial coefficient for secret key, error
        eta_2 (int): maximum absolute polynomial coefficient for noise, error 1, error 2
    returns: 
        A (int[k][k][n] aka polynomial[k][k]): k by k matrix of polynomials. Part of the public key
        t (int[k][n] aka polynomial[k]): polynomial vector. part of the public key.
        s (int[k][n]): polynomial vector representing the secret key
    """ 
    A = (np.random.random([k, k, n]) * q).astype(int)
    s = naive_sample_binom(n,k, eta_1)
    e = naive_sample_binom(n,k, eta_2)
    t = add_vec(mul_mat_vec_simple(A, s, f, q), e, q)
    return A, s, t

def encrypt(A, t, m_b, q, f, n, k, eta_1 = 2, eta_2 = 2):
            
    """
    encryption function for kyber-PKE. operations are done in polynomial ring R := Z_q[x]/f(x)

    params:
        A (int[k][k][n] aka polynomial[k][k]): k by k matrix of polynomials. Part of the public key
        t (int[k][n] aka polynomial[k]): polynomial vector. part of the public key.
        m_b (int[n]): message to encrypt, represented as a binary polynomial.
        q (int): plain modulo for operations in Z_q
        f (int[n]): reduction polynomial's coefficients
        n (int): number of elements in a polynomial. given by deg(f)
        k (int): number of polynomials in a polynomial vector
        eta_1 (int): maximum absolute polynomial coefficient for secret key, error
        eta_2 (int): maximum absolute polynomial coefficient for noise, error 1, error 2
       
    returns: 
        u (int[k][n] aka polynomial[k]): polynomial vector. part of the ciphertext
        v (int[n]): polynomial. part of the ciphertext
    """ 
    half_q = int(q / 2 + 0.5)
    m = list(map(lambda x: x * half_q, m_b))

    r = naive_sample_binom(n,k, eta_1)
    e_1 = naive_sample_binom(n,k, eta_2)
    e_2 = naive_sample_binom(n,1, eta_2)[0]

    u = add_vec(mul_mat_vec_simple(transpose(A), r, f, q), e_1, q)
    v = sub_poly(add_poly(mul_vec_simple(t, r, f, q), e_2, q), m, q)

    return u, v


def decrypt(s, u, v, f, q):
    """
    decryption function for kyber-PKE. operations are done in polynomial ring R := Z_q[x]/f(x)

    params:
        s (int[k][n]): polynomial vector representing the secret key
        u (int[k][n] aka polynomial[k]): polynomial vector. part of the ciphertext
        v (int[n]): polynomial. part of the ciphertext 
        f (int[n]): reduction polynomial's coefficients
        q (int): plain modulo for operations in Z_q
 
    returns: 
        m_b (int[n]): polynomial representing plaintext message
    """ 
    m_n = sub_poly(v, mul_vec_simple(s, u, f, q), q)

    half_q = int(q / 2 + 0.5)
    def round(val, center, bound):
        dist_center = np.abs(center - val)
        dist_bound = min(val, bound - val)
        return center if dist_center < dist_bound else 0

    m_n = list(map(lambda x: round(x, half_q, q), m_n))
    m_b = list(map(lambda x: x // half_q, m_n))
    
    return m_b


# TESTING FUNCTIONS HERE:
def test_enc_dec(N, k, f, q):
    degree_f = len(f) - 1
    failed = 0
    # A = (np.random.random([k, k, degree_f]) * q).astype(int)
    # s = (np.random.random([k, degree_f]) * 3).astype(int) - 1
    # e = (np.random.random([k, degree_f]) * 3).astype(int) - 1
    # t = add_vec(mul_mat_vec_simple(A, s, f, q), e, q)

    
    for i in range(N):
        A,s,t = naive_keygen(q, f, degree_f, k)
        # print(A)
        # print(s)
        # print(t)
        m_b = (np.random.random(degree_f) * 2).astype(int)

        r = (np.random.random([k, degree_f]) * 3).astype(int) - 1
        e_1 = (np.random.random([k, degree_f]) * 3).astype(int) - 1
        e_2 = (np.random.random([degree_f]) * 3).astype(int) - 1

        u, v = encrypt(A, t, m_b, q, f, degree_f, k)
        m_b2 = decrypt(s, u, v, f, q)

        if m_b.tolist() != m_b2:
            failed += 1
    
    print(f"[k={k}, f={f}, q={q}] Test result: {failed}/{N} failed decryption!")


test_enc_dec(100, 2, [1, 0, 0, 0, 1], 67)