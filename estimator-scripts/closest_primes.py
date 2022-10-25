 
import math
from sympy import nextprime
from sympy import prevprime


def closest_primes(a):
    bigger = nextprime(a)
    smaller = prevprime(a)
    return [bigger, smaller]


Q = 2**36
print("Q: " + str(Q))
primes = closest_primes(Q)
print("Closest Prime to Q: " + str(primes))
print("Distances from Q: " + str(primes[0]) + " = Q+" + str((primes[0]-Q)) + ", " + str(primes[1]) + " = Q-" + str(-(primes[1]-Q))) 