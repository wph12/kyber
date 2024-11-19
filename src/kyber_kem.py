from kyber_pke import *
from utils import *
from hashlib import shake_128, sha3_256, shake_256, sha3_512
import secrets

class kyber_kem():
    def __init__(self, q, n, k, eta_1, eta_2, d):
        """
        Toy implementation of ML-KEM. 
        Compression and NTT have NOT been implemented.
        All operations are done in polynomial ring R := Z_q[x]/f(x).
        NOTE: KEM implementation only works if n is a multiple of 8 due to bytewise implementation. Sorry!

        params:
            q (int): plain modulo for operations in Z_q
            f (int[n+1]): reduction polynomial's coefficients
            n (int): number of elements in a polynomial. given by deg(f). Should be a multiple of 8
            k (int): number of polynomials in a polynomial vector
            d: number of bits for each integer coefficient. this is given by ceil(log_2 (q)). in the default ML-KEM with q=3329, bits_per_int = 12 bits
            eta_1 (int): maximum absolute polynomial coefficient for secret key, error
            eta_2 (int): maximum absolute polynomial coefficient for noise, error 1, error 2
        """
        if(n < 8 or n % 8 != 0):
            raise ValueError("Sorry, code only works if n is a multiple of 8 due to this project's byte-wise implementation")
        self.n = n


        func = [0] * (n+1)
        func[0] = 1
        func[n] = 1
        self.f = func

        self.q = q
        self.k = k
        self.eta_1 = eta_1
        self.eta_2 = eta_2
        self.d = d

    def _H(self,input_bytes):
        """
        Hash function H. uses sha3_256.

        params: 
            input_bytes(bytes) : input bytes to hash
            num_bytes(int): number of bytes in the output
        returns:
            hash digest of length num_bytes
        """
        return sha3_256(input_bytes).digest()

    def _J(self,input_bytes, num_bytes):
        """
        Hash function J. uses shake_256.

        params: 
            input_bytes(bytes) : input bytes to hash
            num_bytes(int): number of bytes in the output
        returns:
            hash digest of length num_bytes
        """
        num_bytes = self.n//8
        return shake_256(input_bytes).digest(num_bytes)

    def _G(self,input_bytes):
        """
        Hash function G. uses sha3_512.

        params: 
            input_bytes(bytes) : input bytes to hash
            num_bytes(int): number of bytes in the output
        returns:
            G_digest(bytes,bytes): tuples of hash digests of K, R
        """
        num_bytes_K = self.n//8
        G_digest = sha3_512(input_bytes).digest()
        return G_digest[:num_bytes_K], G_digest[num_bytes_K:num_bytes_K + 32]

    def get_random_bytes(self,num_bytes):
        """
        get random bytestream 

        params:
            num_bytes (int): number of bytes
        
        returns: 
            random bytestream containing num_bytes bytes
        """
        return secrets.token_bytes(num_bytes)
    
    def expand(self, rho):
        """
        expand rho into A by hashing rho with a counter

        params:
            rho (bytes): random bitstream containing n bits

        returns:
            A (int[k][k][n] or polynomial[k][k]): k by k matrix of polynomials. Part of the public key
        """
        A = [[[0]* self.n for _ in range(self.k)] for _ in range(self.k)]
        num_bytes = (self.d * self.n) // 8
        for i in range(self.k):
            for j in range(self.k):
                bytes_ = shake_128(rho + bytearray(i+1) + bytearray(j+1)).digest(num_bytes) #expand with SHAKE-128
                A[i][j] = bytes_to_poly(bytes_, self.d, self.n, self.q)
        
        return A


    def kem_keygen(self):
        """
        key generation function for kyber-KEM. 

        returns: 
            ek (bytes, bytes): (rho, t)
                rho (bytes): random bitstream containing n bits to derive matrix A. part of the public key
                t (bytes): bytes representing polynomial vector. part of the public key.
            dk (dict): dictionary of private key parts. contains s, ek, H(ek), z_
                s (int[n][k] or polynomial[k]): polynomial vector. secret key.
                ek (bytes): ek from above
                H(ek) (bytes): hash of ek using hash function H
                z_ (bytes): z_. used for decapsulation failures
        """
        
        z_ =  self.get_random_bytes(self.n//8)

        #PKE Encryption
        rho = self.get_random_bytes(self.n//8)
        A = self.expand(rho)
        s = sample_binom(self.n,self.k, self.eta_1)
        e = sample_binom(self.n,self.k, self.eta_2)
        t = add_vec(mul_mat_vec_simple(A, s, self.f, self.q), e, self.q)

        t_bytes = poly_vector_to_bytes(t,self.d,self.n,self.q)
        ek = (rho, t_bytes)
        dk = {
            "s" : s,
            "ek": ek,
            "H(ek)" : self._H(rho + t_bytes),
            "z_" : z_
        }

        return ek, dk

    def kem_encrypt(self, m, ek, R):
        """
        key encrypt. used in encapsulation/decapsulation for kyber-KEM

        params:
            m (bytes): message plaintext
            ek ((bytes, bytes)): PKE encryption key
            R (bytes): seed for deterministic psuedorandom sampling from CBD

        returns: 
            u_bytes (bytes): bytes representing polynomial vector. part of the ciphertext
            v_bytes (bytes): bytes representing polynomial. part of the ciphertext
        """
        m_poly = bytes_to_poly(m,1,self.n,self.q)
        rho, t_bytes = ek
        A = self.expand(rho)
        t = bytes_to_poly_vector(t_bytes, self.d , self.n ,self.q, self.k)
        
        r = naive_sample_binom(self.n, self.k, self.eta_1, R, 0) #we can use naive_sample_binom as psuedorandom generator. 
        e_1 = naive_sample_binom(self.n, self.k, self.eta_2, R, 1) #its fine here since R is securely random
        e_2 = naive_sample_binom(self.n, self.k, self.eta_2, R, 2)[0]

        half_q = int(self.q / 2 + 0.5)
        rounded_m = list(map(lambda x: x * half_q, m_poly))
        u = add_vec(mul_mat_vec_simple(transpose(A), r, self.f, self.q), e_1, self.q)
        v = sub_poly(add_poly(mul_vec_simple(t, r, self.f, self.q), e_2, self.q), rounded_m, self.q)

        u_bytes = poly_vector_to_bytes(u, self.d, self.n, self.q)
        v_bytes = poly_to_bytes(v, self.d, self.n, self.q)

        return u_bytes, v_bytes



    def kem_encapsulate(self, ek, bad_rng = False, fixed_bytes = None): 
        """
        key encapsulation function for kyber-KEM

        params:
            ek ((bytes, bytes)): PKE encryption key
            bad_rng (bool): whether to use bad RNG when generating message

        returns: 
            u_bytes (bytes): bytes representing polynomial vector. part of the ciphertext
            v_bytes (bytes): bytes representing polynomial. part of the ciphertext
            K (bytes): symmetric key to be used
        """
        rho, t_bytes = ek
        h = self._H(rho + t_bytes)

        if(bad_rng):
            m = fixed_bytes + self.get_random_bytes(self.n//8 - len(fixed_bytes))
        else:
            m = self.get_random_bytes(self.n//8)


        K, R = self._G(m + h)

        u_bytes, v_bytes = self.kem_encrypt(m, ek, R)

        return u_bytes, v_bytes, K
    

    
    def kem_decapsulate(self, u_bytes, v_bytes, dk):
        """
        key decapsulation function for kyber-KEM. 

        params: 
            u_bytes (bytes): bytes of polynomial u from ciphertext
            u_bytes (bytes): bytes of polynomial v from ciphertext

            dk (dict): dictionary of private key parts. contains s, ek, H(ek), z_
                s (int[n][k] or polynomial[k]): polynomial vector. secret key.
                ek (bytes): ek from above
                H(ek) (bytes): hash of ek using hash function H
                z_ (bytes): z_. used during decapsulation failures
        returns:
            K (bytes): secret key OR garbage random values
        """
        u = bytes_to_poly_vector(u_bytes, self.d, self.n, self.q, self.k)
        v = bytes_to_poly(v_bytes, self.d, self.n, self.q)
        s = dk['s']
        ek = dk['ek']
        z_ = dk['z_']
        m_poly = decrypt(s, u, v, self.f, self.q)
        m = poly_to_bytes(m_poly, 1, self.n, self.q)
        
        K_prime, R_prime = self._G(m + dk["H(ek)"])

        u_bytes_prime, v_bytes_prime = self.kem_encrypt(m, ek, R_prime)

        K_bar = self._J(z_, self.n//8)
        
        if (u_bytes == u_bytes_prime and v_bytes==v_bytes_prime):
            return K_prime
        else:
            print("TOP SECRET: DECAPSULATION FAILURE ENCOUNTERED (THIS IS ONLY KNOWN BY ALICE)")
            return K_bar

def main():
    KYBER_512 = kyber_kem(3329, 256, 2, 2, 2, 12)
    ek, dk = KYBER_512.kem_keygen()
    rho, t = ek
    print("PUBLIC KEY:")
    print_bytes(rho + t)

    u_bytes, v_bytes, K = KYBER_512.kem_encapsulate(ek)
    print("CIPHERTEXT:")
    print_bytes(u_bytes + v_bytes)
    print("SECRET KEY GENERATED: ")
    print_bytes(K)

    K_prime = KYBER_512.kem_decapsulate(u_bytes, v_bytes, dk)
    print("SECRET KEY DECAPSULATED: ")
    print_bytes(K_prime)

if __name__ == "__main__":
    main()