import random
import numpy as np

def bits2int(bits):
    """
    Convert an array of bits to integer
    """
    r = 0
    s = 1
    for b in bits:
        if b & 1:
            r += s
        s <<= 1
    return r

def parity(v):
    """
    Count number of '1' (modulo 2) in binary representation of 'v'
    """
    return bin(v).count('1') & 1

def normalize_state(state):
    """
    Normalize a quantum state
    
    The norm of the state is equal to unity, and the first nonzero element
    is real and positive.
    """
    norm = np.linalg.norm(state)
    for s in state:
        if s != 0:
            norm *= s/abs(s)
            break;
    return state/norm

def state_support(state):
    """
    Count number of nonzero elements of the state
    """
    n = 0
    for s in state:
        if s != 0:
            n += 1
    return n

def gf2_rank(rows):
    """
    Find rank of a matrix over GF(2)

    The rows of the matrix are given as nonnegative integers, thought
    of as bit-strings.

    This function modifies the input list. Use gf2_rank(rows.copy())
    instead of gf2_rank(rows) to avoid modifying rows.
    
    Source: https://stackoverflow.com/a/56858995
    """
    rank = 0
    while rows:
        pivot_row = rows.pop()
        if pivot_row:
            rank += 1
            lsb = pivot_row & -pivot_row
            for index, row in enumerate(rows):
                if row & lsb:
                    rows[index] = row ^ pivot_row
    return rank

def gf2_mul_mat_vec(m, v):
    """
    Multiply a matrix by a vector over GF(2)

    The matrix 'm' is stored as an array of nonnegative integers,
    representing rows. The vector 'v' is a single
    integer, representing a bit-string
    """
    return bits2int(map(parity, m & v))

def gf2_mul_vec_vec(v1, v2):
    """
    Calculate a dot product between two vectors over GF(2)

    Both vectors are represented by integers.
    """
    return parity(v1 & v2)

def qbinomial(n, k, q = 2):
    """
    Calculate q-binomial coefficient
    """
    c = 1
    for j in range(k):
        c *= q**n - q**j
    for j in range(k):
        c //= q**k - q**j
    return c

def number_of_states(n):
    """
    Calculate the number of n-qubit states having 2**k nonzero elements
    
    This function returns a list of values for k = 0, ..., n.
    The total number of n-qubit stabilizer states can be obtained by
    sum(number_of_states(n)).
    """
    return [2**n * 2**((k+1)*k//2) * qbinomial(n, k) for k in range(n+1)]

def random_stabilizer_state(n):
    """
    Generate random n-qubit stabilizer state
    """
    dtype = np.int32
    nmax = np.log2(np.iinfo(dtype).max + 1)
    if not (0 < n and n <= nmax):
        raise ValueError('Number of qubits must be in range(1, %d)!' % (nmax + 1))
    dimn = 2**n
    state = np.zeros(dimn, dtype=np.complex64)
    k = random.choices(range(0, n+1), weights=number_of_states(n))[0]
    dimk = 2**k
    if k == 0:
        state[random.randrange(dimn)] = 1
        return state
    rank = 0
    while rank < k:
        R = np.random.randint(dimk, size=n, dtype=dtype)
        rank = gf2_rank(list(R))
    t = np.random.randint(dimn, dtype=dtype)
    Q = np.random.randint(dimk, size=k, dtype=dtype)
    c = np.random.randint(dimk, dtype=dtype)
    for x in range(dimk):
        y = gf2_mul_mat_vec(R, x) ^ t # y = R @ x + t
        ib = gf2_mul_vec_vec(c, x) # ib = c @ x, 'imaginary' bit
        mb = gf2_mul_vec_vec(x, gf2_mul_mat_vec(Q, x)) # mb = x @ Q @ x, 'minus' bit
        state[y] = 1j**ib * (-1)**mb
    state = normalize_state(state)
    return state
