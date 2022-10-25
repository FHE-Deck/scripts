import math

from sympy import N

from ntt_primes import Q
   

# A Function to estimate the magnitude in the form of a floating point number:
def estimate_magnitude_float(A):
    if A == 0.0:
        return A
    return math.log(A, 2)

class ntru_param: 
    P = None
    Q = None
    N = None
    gadget_basis = None
    gadget_ell = None
    stddev_1 = None
    stddev_2 = None
    variance_sk = None
    expected_sk = None
    hamming_sk = None
    ring = None

    def __init__(self, P, Q, N, key_type, ring, hamming_sk, gadget_basis, stddev_1, stddev_2): 
        # Power of two Modulus
        self.P = P
        # Prime Modulus  
        self.Q = Q
        # Ring Degree
        self.N = N
        # Gadget Basis
        self.gadget_basis = gadget_basis
        self.gadget_ell =  math.ceil(math.log(P, gadget_basis)) 
        # Standard Deviation
        self.stddev_1 = stddev_1
        self.stddev_2 = stddev_2
        self.key_type = key_type
        # Variance of the secret key 
        if(key_type == "ternary"):
            self.variance_sk = 2/3.0
            self.expected_sk = 0
        elif(key_type == "binary"):
            self.variance_sk = 1/4.0
            self.expected_sk = 1/2.0
        # Hamming weight of the secret key
        self.hamming_sk = hamming_sk
        # For the ring only two options are considered:
        # - Cyclic X^N - 1
        # - Negacyclic: X^N + 1
        self.ring = ring
 


class lwe_param:
    Q = None
    n = None
    gadget_basis = None
    gadget_ell = None
    stddev = None
    stddev = None
    variance_sk = None
    expected_sk = None
    hamming_sk = None
    def __init__(self,  Q, n, key_type, hamming_sk, gadget_basis, stddev):
        # Modulus
        self.Q = Q
        # Dimension
        self.n = n
        # Gadget Basis
        self.gadget_basis = gadget_basis
        self.gadget_ell =  math.ceil(math.log(Q, gadget_basis)) 
        # Standard Deviation
        self.stddev = stddev
        self.key_type = key_type
        # Variance of the secret key 
        if(key_type == "ternary"):
            self.variance_sk = 2/3.0
            self.expected_sk = 0
        elif(key_type == "binary"):
            self.variance_sk = 1/4.0
            self.expected_sk = 1/2.0
        # Hamming weight of the secret key
        self.hamming_sk = hamming_sk


def estimate_paramaters(t_1, t_2, lwe_param, ntru_param): 
    # Power of two modulus
    P = ntru_param.P
    # Prime Modulus  
    Q = ntru_param.Q
    # Ring degree
    N = ntru_param.N
    # LWE dimension (note that the keyswitching is given in the Q modulus)
    n = lwe_param.n

    #NOTE LWE modulus (before blind rotating)
    # For negacyclic its 2*N, for Cyclic its N
    if ntru_param.ring == "cyclic":
        q = N
    elif ntru_param.ring == "negacyclic":
        q = 2*N
    else:
        raise Exception("Suported only cyclic and negacyclic Rings")

    #NOTE: NTRU secret key parameters:
    # For ternary secret keys: we have variance 2/3.0, and expected value 0.
    # We can also take larger variance for the secret keys 
    # (that isn't the bottle neck - the bottle nech is the hamming weight)
    variance_sk = ntru_param.variance_sk
    #expected_sk = ntru_param.expected_sk
    hamming_sk = ntru_param.hamming_sk

    #NOTE: LWE secret key parameters:
    # This is the set size for the blind rotation
    # For binary its u=1
    # For ternary its u=2
    if lwe_param.key_type == "binary":
        u = 1
    elif lwe_param.key_type == "ternary":
        u = 2
    else:
        raise Exception("Supported only Binary or Ternary LWE Keys")
    # Variance of the LWE secret key.
    # This is used in mod switching from q to 2N and 2N_F
    # For Ternary secret key the variance is 2/3 and expectation is 0
    # # For binary variance is 1/4 and expectation is 1/2
    variance_s = lwe_param.variance_sk
    expected_s = lwe_param.expected_sk
    # At default hamming weigth is n
    hamming_s = lwe_param.hamming_sk

    
    #RGSW Basis and ell parameter
    basis_ext_NTRU = ntru_param.gadget_basis
    ell_ext_NTRU = ntru_param.gadget_ell 
    print("- ell_ext_NTRU: " + str(ell_ext_NTRU)) 
    #key switching Basis and ell parameter for second blind rotation
    basis_keySwitchKey = lwe_param.gadget_basis
    ell_keySwitchKey = lwe_param.gadget_ell  
    print("- ell_keySwitchKey: " + str(ell_keySwitchKey)) 

    # Fresh Standard Deviation
    stddev_1 = ntru_param.stddev_1
    stddev_2 = ntru_param.stddev_2
    # This bound 
    variance_1  = stddev_1**2
    variance_2  = stddev_2**2
    #Variance on the blind rotation RGSW samples
    var_brKey_1 = variance_1
    var_brKey_2 = variance_2

    # Standard Deviation of the Key Switching key
    stddev_ksk = lwe_param.stddev
    #variance on the LWE KeySwitch key samples
    var_keySwitchKey =  stddev_ksk**2


    # Error variance of the basic blind rotation
    B_BR_1 = var_brKey_1 +  n * u * N  * ell_ext_NTRU * (basis_ext_NTRU**2)/4 * var_brKey_1
    B_BR_2 = var_brKey_2 +  n * u * N  * ell_ext_NTRU * (basis_ext_NTRU**2)/4 * var_brKey_2
    B_BR_prime_mod_1 = (Q/P)**2 * B_BR_1

    B_BR_prime_mod_2 = (Q/P)**2 * B_BR_2  

    #B_BR =  n * u * N  * ell_ext_NTRU * (basis_ext_NTRU**2-1)/12 * var_brKey
    # Below is a Blind rotation that is (most likely) faster but has larger error
    # B_BR = var_brKey +  5 * n * N  * ell_ext_NTRU * (basis_ext_NTRU**2-1)/12 * var_brKey 

    B_KS = N * ell_keySwitchKey * (basis_keySwitchKey**2)/4 * var_keySwitchKey
    
    # The error variance when decrypting right after blind rotation
    # Note: I take a version where: we have c = e_1* g/f + e_2
    # Where e_2 = 0, and g is just like f. So, the decryption error is: e_1 * g.
    # Meaning, that its B_BR * variance_sk * hamming_sk, because each coef of the outcome is: 
    # sum_{i=1}^hamming_sk var(sk_i) * Var(e_1,j)
    # If you want to have e_2, then we have 2times the noise
    B_out =  hamming_sk * variance_sk * (B_BR_prime_mod_1 + B_BR_prime_mod_2)

    # The error variance of a ciphertext after key switching from NTRU to LWE
    # The same as above with B_out
    B_ksk = B_out + B_KS

    # The error after modulus switching from Q to q (just before blind rotation)
    B_in = (q**2/Q**2) * B_ksk + hamming_s * (variance_s)  


    print("")
    print("=== Correctness:")

    print("- Standard Deviation After blind rotation: " + str(math.sqrt(B_out)))
    print("- Standard Deviation Before blind rotation: " + str(math.sqrt(B_in)))

    # Computing the probability of having correct ciphertext after blind rotation
    #NOTE: Need to handle the case where the variances are zero (just in case) in whihc case we would divide by zero
    if B_out != 0.0:
        pr_corr_B_out = math.erf(Q/(2.0 * t_1 * math.sqrt(2) * math.sqrt(B_out)))
    else:
        pr_corr_B_out = 1.0 
    print("- Pr correctness of after blind rotation (for t_1 plaintext space): " + str(pr_corr_B_out))
    pr_error = 1.0 -  pr_corr_B_out 
    print("    Pr error: " + str(pr_error) + " approx 2.0**(-" + str(estimate_magnitude_float(pr_error)) + "))")


    # Probability of having a correct LWE ciphertext before blind rotating it
    if B_in != 0.0:
        pr_corr_B_in = math.erf(q/(2.0 * t_2 * math.sqrt(2) * math.sqrt(B_in)))
    else:
        pr_corr_B_in = 1.0 
    print("- Pr correctness of LWE ciphertext before blind rotation (for t_2 plaintext space): " + str(pr_corr_B_in))
    pr_error = 1.0 -  pr_corr_B_in 
    print("    Pr error: " + str(pr_error) + " approx 2.0**(-" + str(estimate_magnitude_float(pr_error)) + "))")


    share_of_key_s = (hamming_s * variance_s  * 100) / B_in
    print("- Share of the LWE secret key s in the variance: " + str(share_of_key_s) + "%")
  
    """
    print("============= Fatigue point at: " + str(estimate_magnitude_int(int(0.004 * N**(2.484)))))

    # Attempt at estimation
    S = math.log(math.sqrt(2 * N * variance_sk), N)
    print("S: " + str(S))
    Q_param = math.log(Q, N)
    print("Q_param: " + str(Q_param))
    print("Q_param: " + str(math.log(Q, 2)/math.log(N, 2)))
    beta = ((8 * S)/(Q_param**2 + 1) - 0.48) * N/2
    print("beta: " + str(beta))
    beta = 405
    print("beta: " + str(beta))
    exp_cost = 0.292 * beta + 16.4 + math.log(8 * N, 2)
    T = 2**exp_cost
    print("============= Security Level Estimate: " + str(estimate_magnitude_int(int(T))))

    print("From LWE Estimator cost model: " + str(estimate_magnitude_int(int(1.62651036352611e42))))
    # it actually is the same

    """ 
    print("")
    print("--- Size ================")

    # How many bytes I need to store a single element in Q:
    bits_Q = math.ceil((math.log2(Q)))
    bytes_Q = math.ceil(bits_Q/8)  
    # How many bytes I need to store a single element in Q:
    bits_P = math.ceil((math.log2(P)))
    bytes_P = math.ceil(bits_P/8)  
    # How many bytes I need to store a single element in q:
    bits_q = math.ceil((math.log2(q)))
    bytes_q = math.ceil(bits_q/8)  

    # Size of a NTRU the ciphertext: 
    ct_size_bytes = N * bytes_Q
    ct_size_kilo_bytes = ct_size_bytes/1000.0
    print("- NTRU Ciphertext size in KB: " + str(ct_size_kilo_bytes))
    #Size of a LWE ciphertexts:
    ct_LWE_size_bytes = (n+1) * bytes_q
    ct_LWE_size_kilo_bytes = ct_LWE_size_bytes/1000.0
    print("- LWE Ciphertext size in KB: " + str(ct_LWE_size_kilo_bytes))
 
    ct_KSK_size_bytes = (n+1) * bytes_Q
    ct_KSK_size_kilo_bytes = ct_KSK_size_bytes/1000.0
    ks_key_size_kilo_bytes =  N * ell_keySwitchKey * ct_KSK_size_kilo_bytes
    ks_key_size_mega_bytes = ks_key_size_kilo_bytes/1000.0
    print("- Key Switching Key is MB: " + str(ks_key_size_mega_bytes))

    
    # Size of a NTRU (Power of two) ciphertext: 
    ct_size_bytes_P = N * bytes_P 
    ct_size_kilo_bytes_P = ct_size_bytes_P/1000.0
    # Key Switching key (N * ell_keySwitchKey * LWE sample): 
    # Size of the blind rotation keys n * u * (NTRU samples)  
    br_key_size_kilo_bytes = n * u * ell_ext_NTRU * ct_size_kilo_bytes_P
    br_key_size_mega_bytes = br_key_size_kilo_bytes/1000.0
    
    print("- Blind Rotation Key size in MB: " + str(br_key_size_mega_bytes))

    sum_of_all_key_mega_bytes = br_key_size_mega_bytes + ks_key_size_mega_bytes
    sum_of_all_key_giga_bytes = sum_of_all_key_mega_bytes/1000.0

    print("- Sum of all keys in MB: " + str(sum_of_all_key_mega_bytes))
    print("- Sum of all keys in GB: " + str(sum_of_all_key_giga_bytes)) 


 
 


print("============================== Set: Set with N=11 - Binary ====================================== ")
# Message space after Blind Rotation
k_1 = 6
t_1 = 2**k_1
#Message space before blind rotation
k_2 = 6
t_2 = 2**k_2
# Power of two NTRU modulus
P = 2**30
# Prime NTRU modulus
Q = 2**25-39
# Ring Size (Should be prime)
N = 2**11-9
ntru_key_type = "ternary"
ring = "negacyclic"
hamming_sk = N 
ntru_g_basis = 2**6
ntru_stddev_1 = math.sqrt(2/3)
ntru_stddev_2 = 1/4
ntru_par = ntru_param(P, Q, N, ntru_key_type, ring, N, ntru_g_basis, ntru_stddev_1, ntru_stddev_2)
lwe_key_type = "binary"
n = 2**9 + 125
hamming_s = n
lwe_g_basis = 2
lwe_stddev = 2**10 
lwe_par = lwe_param(Q, n, lwe_key_type, hamming_s, lwe_g_basis, lwe_stddev)
#estimate_paramaters(t_1, t_2, lwe_par, ntru_par)


 

print("============================== Set: Set with N=12 - Binary ====================================== ")
# Message space after Blind Rotation
k_1 = 9
t_1 = 2**k_1
#Message space before blind rotation
k_2 = 9
t_2 = 2**k_2
# Power of two NTRU modulus
P = 2**45
# Prime NTRU modulus
Q = 2**33-9
# Ring Size (Should be prime)
N = 2**12-3
ntru_key_type = "ternary"
ring = "negacyclic"
hamming_sk = N 
ntru_g_basis = 2**15
ntru_stddev_1 = math.sqrt(2/3)
ntru_stddev_2 = 1/4 
ntru_par = ntru_param(P, Q, N, ntru_key_type, ring, N, ntru_g_basis, ntru_stddev_1, ntru_stddev_2)
lwe_key_type = "binary"
n = 2**9 + 238
hamming_s = n
lwe_g_basis = 2
lwe_stddev = 2**14 
lwe_par = lwe_param(Q, n, lwe_key_type, hamming_s, lwe_g_basis, lwe_stddev)
#estimate_paramaters(t_1, t_2, lwe_par, ntru_par)




print("============================== Set: Set with N=13 - Binary ====================================== ")
# Message space after Blind Rotation
k_1 = 11
t_1 = 2**k_1
#Message space before blind rotation
k_2 = 11
t_2 = 2**k_2
# Power of two NTRU modulus
P = 2**45
# Prime NTRU modulus
Q = 2**34-41
# Ring Size (Should be prime)
N = 2**13-3
ntru_key_type = "ternary"
ring = "negacyclic"
hamming_sk = N 
#ntru_g_basis = 2**16
ntru_g_basis = 2**15
ntru_stddev_1 = math.sqrt(2/3)
ntru_stddev_2 = 1/4 
ntru_par = ntru_param(P, Q, N, ntru_key_type, ring, N, ntru_g_basis, ntru_stddev_1, ntru_stddev_2)
lwe_key_type = "binary"
n = 2**9 + 315
hamming_s = n
lwe_g_basis = 2
lwe_stddev = 2**14 
lwe_par = lwe_param(Q, n, lwe_key_type, hamming_s, lwe_g_basis, lwe_stddev)
#estimate_paramaters(t_1, t_2, lwe_par, ntru_par)




print("============================== Set: Set with N = 14 - Binary ====================================== ")
# Message space after Blind Rotation
k_1 = 11
t_1 = 2**k_1
#Message space before blind rotation
k_2 = 11
t_2 = 2**k_2
# Power of two NTRU modulus
P = 2**45
# Prime NTRU modulus
Q = 2**36-5
# Ring Size (Should be prime)
N = 2**14-3
ntru_key_type = "ternary"
ring = "negacyclic"
hamming_sk = N 
ntru_g_basis = 2**15    
ntru_stddev_1 = math.sqrt(2/3)
ntru_stddev_2 = 1/4 
ntru_par = ntru_param(P, Q, N, ntru_key_type, ring, N, ntru_g_basis, ntru_stddev_1, ntru_stddev_2)
lwe_key_type = "binary"
n = 2**9 + 390
hamming_s = n
lwe_g_basis = 2
lwe_stddev = 2**14 
lwe_par = lwe_param(Q, n, lwe_key_type, hamming_s, lwe_g_basis, lwe_stddev)
estimate_paramaters(t_1, t_2, lwe_par, ntru_par)


 