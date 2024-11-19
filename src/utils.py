import math
import binascii

def bytes_to_poly(input_bytes, bits_per_int, n, q):
    """
    converts bytestream into polynomial  
    params:
        input_bytes (bytes): bytestream
        bits_per_int: number of bits per integer coefficient. this is given by ceil(log_2 (q)). in the default ML-KEM with q=3329, bits_per_int = 12 bits
        n (int): number of elements in a polynomial, or degree of the function polynomial
        q (int): plain modulo for operations in Z_q

    return:
        polynomial (int[n])
    """

    if (len(input_bytes) != math.ceil((bits_per_int * n) / 8)):
        raise ValueError("byte length mismatch")

    poly = [0 for _ in range(n)]
    mask = (1<<bits_per_int)-1
    b_int = int.from_bytes(input_bytes, "little")
    for i in range(n):
        poly[i] = (b_int & mask) % q
        b_int >>= bits_per_int
    return poly
    


def poly_to_bytes(poly, bits_per_int, n, q):
    """
    converts polynomial into bytes

    params:
        polynomial (int[n]): polynomial to convert
        bits_per_int: number of bits per integer coefficient. this is given by ceil(log_2 (q)). in the default ML-KEM with q=3329, bits_per_int = 12 bits
        n (int): number of elements in a polynomial, or degree of the function polynomial
        q (int): plain modulo for operations in Z_q

    return:
        bytestream
    """
    if len(poly) != n:
        raise ValueError("Polynomial length does not match specified n.")
    
    big_int = 0
    for i in range(n):
        coeff = poly[n-i-1] % q
        if coeff >= (1 << bits_per_int):
            raise ValueError("Coefficient is too large for specified bits per integer.")
        
        big_int = (big_int << bits_per_int) | coeff
    
    num_bytes = math.ceil((bits_per_int * n) / 8)
    
    return big_int.to_bytes(num_bytes, 'little')

def poly_vector_to_bytes(poly_vector, bits_per_int, n, q):
    """
    converts polynomial vector into bytes

    params:
        poly_vector(int[k][n]): polynomial vector to convert
        bits_per_int: number of bits per integer coefficient. this is given by ceil(log_2 (q)). in the default ML-KEM with q=3329, bits_per_int = 12 bits
        n (int): number of elements in a polynomial, or degree of the function polynomial
        q (int): plain modulo for operations in Z_q

    return:
        bytestream
    """
    bytestream = b''
    for poly in poly_vector:
        bytestream += poly_to_bytes(poly, bits_per_int, n, q)

    return bytestream

def bytes_to_poly_vector(bytes_, bits_per_int, n, q, k):
    """
    converts polynomial vector into bytes

    params:
        bytes: bytestream to convert
        bits_per_int: number of bits per integer coefficient. this is given by ceil(log_2 (q)). in the default ML-KEM with q=3329, bits_per_int = 12 bits
        n (int): number of elements in a polynomial, or degree of the function polynomial
        q (int): plain modulo for operations in Z_q
        k (int): number of polynomials in a vector

    return:
        poly_vector(int[k][n]): polynomial vector
    """
    poly_vector = [[0] * n] * k
    num_bytes_per_poly = math.ceil((bits_per_int * n) / 8)
    
    for i in range(k):
        start_idx = i * num_bytes_per_poly
        end_idx = start_idx + num_bytes_per_poly
        poly_vector[i] = bytes_to_poly(bytes_[start_idx:end_idx], bits_per_int, n, q)
    
    return poly_vector
    

def print_bytes(input_bytes):
    print(binascii.hexlify(input_bytes))

# def string_to_bytes(string):
#     s=binascii.unhexlify(string)
#     return bytearray(s)

#TESTS
def main():
    bytez = poly_to_bytes([11,10,58], 12, 3,67)

    print_bytes(bytez)
    print(bytes_to_poly(bytez, 12, 3,67))

    vector_bytez = poly_vector_to_bytes([[11,10,58],[5,6,7]], 12, 3, 67)
    print_bytes(vector_bytez)
    print(bytes_to_poly_vector(vector_bytez, 12, 3, 67, 2))

if __name__ == "__main__":
    main()
