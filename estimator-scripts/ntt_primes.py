import math
from sympy import nextprime
from sympy import prevprime
from sympy import isprime
   
 

def closest_primes(a):
    bigger = nextprime(a)
    smaller = prevprime(a)
    return [bigger, smaller]

#NOTE: Actually shouldn't it be the number of roots of unity instead of N?
def next_NTT_modulus(N, start_Q):
    t = 2
    while t * N + 1 < start_Q:
        t = t * 2

    Q = t * N + 1
    """
    Note that Q - 1 = 0 mod N is
    Q - 1 = t * N, for some t
    Thus Q = t * N + 1, for some t.
    Now we need to find the largest t s.t. Q is prime and smaller than the limit.
    """
    while (isprime(Q) == False):
        t = t + 1
        Q = t * N + 1
    return Q

def prev_NTT_modulus(N, start_Q):
    t = 2
    while t * N + 1 < start_Q:
        t = t * 2

    Q = t * N + 1
    """
    Note that Q - 1 = 0 mod N is
    Q - 1 = t * N, for some t
    Thus Q = t * N + 1, for some t.
    Now we need to find the largest t s.t. Q is prime and smaller than the limit.
    """
    while (isprime(Q) == False):
        t = t - 1
        Q = t * N + 1
    return Q


def ntt_friendly_moduli(N, start_Q, end_Q):
    current = start_Q
    out = []
    while(current < end_Q):
        current = next_NTT_modulus(N, current)
        out.append(current)
    return out


def decomposition(x, base):
    k = math.ceil(math.log(x, base)) 
    powers = [0]*k
    for i in range(0,k):
        powers[i] = x % base
        x = int(x/base)
    return powers

def l_two_norm(inp):
    norm = 0
    for x in inp:
        norm = norm + x**2
    return math.sqrt(norm)

start_Q = 2**47
end_Q = 2**49
base = 2**12
# If you wanna have a NTT friendly prime for the negacyclic ring Z_Q[X]/(X^N + 1), then we have 2 * the_actual_N
N =  2**11
two_N = 2 * N

#primes = closest_primes(N)
#print("Closest Prime to N: " + str(primes))
#print("Distances from N: N+" + str((primes[0]-N)) + ", N-" + str(-(primes[1]-N))) 

x = prev_NTT_modulus(two_N, 2**54)
print("Smaller Q: " + str(x))
print("Distance from start_Q: " + str(start_Q - x))  
print("- Prime: " + str(x))
decomp = decomposition(x, base)
print("  Base-" + str(base) + " decomposition of " + str(x) + ": " + str(decomp)) 
print("  log(" + str(x) + ", 2)): " + str(math.log(x, 2)) + ", ceil(log(" + str(x) + ", " + str(base)+ ")): " + str(math.ceil(math.log(x, base)))) 
print("  L_2 norm: " + str(l_two_norm(decomp))) 

#Q = next_NTT_modulus(N, start_Q)
#print("Larger Q: " + str(Q))
#print("Distance from start_Q: " + str(Q - start_Q)) 
#print("Estimate Magnitude: " + str(math.log(Q, 2)))  

primes = ntt_friendly_moduli(two_N, start_Q, end_Q)
print("NTT Friendly primes within range (" + str(start_Q) + ", " + str(end_Q) + ")")
for x in primes:
    print("- Prime: " + str(x))
    decomp = decomposition(x, base)
    print("  Base-" + str(base) + " decomposition of " + str(x) + ": " + str(decomp)) 
    print("  log(" + str(x) + ", 2)): " + str(math.log(x, 2)) + ", ceil(log(" + str(x) + ", " + str(base)+ ")): " + str(math.ceil(math.log(x, base)))) 
    print("  L_2 norm: " + str(l_two_norm(decomp))) 