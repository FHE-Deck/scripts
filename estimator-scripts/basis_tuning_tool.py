import math

 
def powers_of(A, basis):
    if A == 0:
        return 0
    return math.ceil(math.log(A, basis))


def print_powers_of(modulus, base):
    magnitude = powers_of(modulus, base)
    ell_basis = modulus
    for i in range(1, magnitude):
        ell_basis_new = math.ceil(math.log(modulus, base**(i)))
        if ell_basis_new < ell_basis:
            print("ceil(log(" + str(base) + "**("+ str(i) +"))) = " + str(ell_basis_new))
            ell_basis = ell_basis_new


Q = 2**33
print_powers_of(Q, 2)

