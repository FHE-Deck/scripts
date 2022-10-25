import math
   
# A Function to estimate the magnitude of the input:
def estimate_magnitude_int(A):
    return math.log(A, 2) 

# A Function to estimate the magnitude in the form of a floating point number:
def estimate_magnitude_float(A):
    if A == 0.0:
        return A
    return math.log(A, 2)



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
 
def contribution_of_the_secret_key(bound_out, hamming_weigth_s, variance_vec_s, expected_vec_s):
    #TODO Actually this is wrong!!! Need to change that...
    # Computing the contribution of the secret key
    sk_variance_in_mod_switch = 2/3.0 + 2/3.0 * hamming_weigth_s  * (variance_vec_s + expected_vec_s**2)
    sk_contribution = (sk_variance_in_mod_switch * 100)/ bound_out   
    print("--- Contribution of the mod switch secret key part: " + str(int(sk_contribution)) + "%")




def print_correctness(error_params, q, t): 
    bound = q/(2.0 * t) - float(error_params['mean'])
    pr_error_out = math.erfc(bound/(math.sqrt(2) * error_params['stddev'])) 
    print("[Mean: " + str(float(error_params['mean'])) + "], [stddev: " + str(error_params['stddev']) + "], [log(stddev): " + str(estimate_magnitude_float(error_params['stddev'])) +"], [Pr error: " + str(pr_error_out) + "], [log_2(pr_error_out): " + str(estimate_magnitude_float(pr_error_out)) + "]")
   

def bound_on_Gaussian(stddev, bound):
    pr_error_out = math.erfc(bound/(math.sqrt(2) * stddev))
    print("Pr error: " + str(pr_error_out) + ", log_2(pr_error_out) = " + str(estimate_magnitude_float(pr_error_out)))

def key_switching_error(error_params, lwe_dim, stddev_ksk, basis_ksk, ell_ksk):    
    #key switching Basis and ell parameter for second blind rotation  
    error_params['stddev'] = math.sqrt(error_params['stddev']**2 + lwe_dim * ell_ksk * basis_ksk**2 * stddev_ksk**(2))
    return error_params


def compute_c_epsilon_m(sec_par, m):
    pi = 3.14159
    delta = 2**(-sec_par-1)  
    return math.log(2 * m *  + (2 * m)/delta)/pi  


def compute_gaussian_LHL_stddev(sec_par, bound_on_e, Q, base, m):
    #pi = 3.14159
    #delta = 2**(-sec_par-1)  
    #in_the_sqrt = math.log(2 * m *  + (2 * m)/delta)/pi  
    in_the_sqrt = compute_c_epsilon_m(sec_par, m)
    norm = 0 
    if base**(math.ceil(math.log(Q, base))) == Q:
        norm = math.sqrt(base**2 + 1)  
    else:
        decomp = decomposition(Q, base)
        norm = l_two_norm(decomp)
    return norm * math.sqrt(1 + bound_on_e) * math.sqrt(in_the_sqrt)
 

def compute_classic_LHL_stddev(sec_par, m):
    gamma = 2**(-sec_par-4)
    epsilon = 2 * gamma 
    return 4 * ((1 - gamma) * epsilon**2)**(-1/m)


def compute_stddev_for_gaussian_sampling(sec_par, Q, base): 
    ell = math.ceil(math.log(Q, base))
    stddev = 0
    if base**ell == Q: 
        in_the_sqrt = compute_c_epsilon_m(sec_par, ell)
        stddev = math.sqrt(base**2 + 1) * in_the_sqrt
    else:
        #decomp = decomposition(Q, base)
        #norm = l_two_norm(decomp)
        #pi = 3.14159 
        #delta = 2**(-sec_par-1)  
        #in_the_sqrt = math.log(2 * ell *  + (2 * ell)/delta)/pi  
        in_the_sqrt = compute_c_epsilon_m(sec_par, ell)
        stddev = math.sqrt(2*base) * (2 * base + 1) * in_the_sqrt 
    return stddev


def blind_rotation_error(params):  
    pi = 3.14159
    # RLWE degree
    N = params['N']
    # LWE dimension
    n = params['n']   
    #RGSW Basis and ell parameter   
    var_br = params['stddev_br']**(2) 

    if params['version'] == 'det':
        # This is normal deterministic base decomposition  
        var_RGSW_gadget = params['basis_RGSW']**2
        # Blind rotation
        B_BR = 2 * n *  N  * params['ell_RGSW'] * var_RGSW_gadget  * var_br
        B_BR = math.sqrt(B_BR)
        return {'stddev':B_BR,'mean':0}
    elif params['version'] == 'simul':
        # Simulation: var_RGSW_gadget from the GLHL
        #delta = 2**(-params['cp_secpar']-1) 
        #bound_error_RGSW = params['bound_on_Gaussian'] * params['stddev_br']
        #in_the_sqrt = math.log(4 * params['n'] * params['N'] * params['ell_RGSW'] * (1 + 1/delta))/pi
        #stddev_RGSW_gadget = (params['basis_RGSW'] + 1) * math.sqrt(1 + bound_error_RGSW) * math.sqrt(in_the_sqrt)
        stddev_gaussian_lhl = compute_gaussian_LHL_stddev(params['cp_secpar'], params['bound_on_Gaussian'] * params['stddev_br'], params['Q'], params['basis_RGSW'], 2 * params['n'] * params['N'] * params['ell_RGSW'])
        stddev_gaussian_sampling = compute_stddev_for_gaussian_sampling(params['cp_secpar'], params['Q'], params['basis_RGSW'])
        if(stddev_gaussian_lhl >= stddev_gaussian_sampling):
            stddev_RGSW_gadget = stddev_gaussian_lhl 
            print("For Gadget Multiplication: using stddev form the Gaussian Leftover Hash Lemma. stddev: " + str(stddev_RGSW_gadget))
        else:
            stddev_RGSW_gadget = stddev_gaussian_sampling
            print("For Gadget Multiplication: using stddev form the Gaussian Sampling Algorithm. stddev: " + str(stddev_RGSW_gadget)) 

        var_RGSW_gadget = stddev_RGSW_gadget**2 *  params['basis_RGSW']**2 
        B_BR = 2 * n *  N  * params['ell_RGSW'] * var_RGSW_gadget  * var_br
 
        # Simulation: var_masking from core randomization lemma (we take the max from LHL and GLHL)
        #bound_error_masking = params['bound_on_Gaussian'] * params['stddev_masking']
        # Compute the LHL
        #gamma = 2**(-params['cp_secpar']-4)
        #epsilon = 2 * gamma 
        #stddev_lhl_masking = 4 * ((1 - gamma) * epsilon**2)**(-1/params['ell_masking'])
        stddev_lhl_masking = compute_classic_LHL_stddev(params['cp_secpar'], params['ell_masking'])
        # Compute GLHL
        #in_the_sqrt = math.log(2 * params['ell_masking'] * (1 + 1/gamma))/pi
        #stddev_glhl_masking = (params['basis_masking'] + 1) * math.sqrt(1 + bound_error_masking) * math.sqrt(in_the_sqrt)
        stddev_glhl_masking = compute_gaussian_LHL_stddev(params['cp_secpar'], params['bound_on_Gaussian'] * params['stddev_masking'], params['Q'], params['basis_masking'], params['ell_masking']) 
        # Take the max from both
        #stddev_masking_gadget = max(stddev_lhl_masking, stddev_glhl_masking)
        if(stddev_lhl_masking >= stddev_glhl_masking):
            stddev_masking_gadget = stddev_lhl_masking
            print("For Masking: using stddev form the Classic Leftover Hash Lemma. stddev: " + str(stddev_masking_gadget))
        else:
            stddev_masking_gadget = stddev_glhl_masking
            print("For Masking: using stddev form the Gaussian Leftover Hash Lemma. stddev: " + str(stddev_masking_gadget))
        var_masking_gadget = stddev_masking_gadget**2
        var_masking = params['stddev_masking']**2
        B_masking = params['ell_masking'] * var_masking_gadget * var_masking
        # Put things together and return
        B_BR = B_BR + B_masking
        B_BR = math.sqrt(B_BR)
        return {'stddev':B_BR,'mean':0}
    elif params['version'] == 'flood':
        var_RGSW_gadget = params['basis_RGSW']**2
        # Blind rotation
        B_BR = 2 * n *  N  * params['ell_RGSW'] * var_RGSW_gadget  * var_br

        # Compute the masking error and standard deviation for the gadget  
        #bound_error_masking = params['bound_on_Gaussian'] * params['stddev_masking']
        # Compute the LHL
        #gamma = 2**(-params['cp_secpar']-4)
        #epsilon = 2 * gamma 
        #stddev_lhl_masking = 4 * ((1 - gamma) * epsilon**2)**(-1/params['ell_masking'])
        stddev_lhl_masking = compute_classic_LHL_stddev(params['cp_secpar'], params['ell_masking'])
        # Compute GLHL
        #in_the_sqrt = math.log(2 * params['ell_masking'] * (1 + 1/gamma))/pi
        #stddev_glhl_masking = (params['basis_masking'] + 1) * math.sqrt(1 + bound_error_masking) * math.sqrt(in_the_sqrt)
        stddev_glhl_masking = compute_gaussian_LHL_stddev(params['cp_secpar'], params['bound_on_Gaussian'] * params['stddev_masking'], params['Q'], params['basis_masking'], params['ell_masking'])
        # Take the max from both
        #stddev_masking_gadget = max(stddev_lhl_masking, stddev_glhl_masking)
        if(stddev_lhl_masking >= stddev_glhl_masking):
            stddev_masking_gadget = stddev_lhl_masking
            print("For Masking: using stddev form the Classic Leftover Hash Lemma. stddev: " + str(stddev_masking_gadget))
        else:
            stddev_masking_gadget = stddev_glhl_masking
            print("For Masking: using stddev form the Gaussian Leftover Hash Lemma. stddev: " + str(stddev_masking_gadget))
        var_masking_gadget = stddev_masking_gadget**2
        var_masking = params['stddev_masking']**2
        B_masking = params['ell_masking'] * var_masking_gadget * var_masking
        B_BR = B_BR + B_masking
        B_BR = math.sqrt(B_BR)

        # Compute the flooding error 
        bound_B_BR = params['bound_on_Gaussian'] * B_BR
        # Noise flooding error is a little bit different - we compute the bound B_flood =  2**(params['cp_secpar']/params['washing_cycles']) * bound_B_BR 
        # The washing cycles are the number of times we  run the bootstrapping algorithm
        # But the random variable is uniformly random, and not Gaussian distributed
        B_flood =  2**(params['cp_secpar']/params['washing_cycles']) * bound_B_BR 
        print("flooring error: " + str(math.log(B_flood, 2))) 
        return {'stddev':B_BR,'mean':B_flood}
    else:
        raise Exception("Version " + str(params['version']) + " not supported.")
          



def modulus_reduction(error_param, Q, q, hamming_weigth_s, variance_vec_s, mean_vec_s): 
    error_param['stddev'] = math.sqrt((q/Q * error_param['stddev'])**2 +  (hamming_weigth_s * variance_vec_s)/4)
    error_param['mean'] = q/Q * error_param['mean'] + (1 + hamming_weigth_s * mean_vec_s)/2
    return error_param




def print_ct_and_key_sizes(params):
    # How many bytes I need to store a single element in Q:
    bits_Q = math.ceil((math.log2(params['Q'])))
    bytes_Q = math.ceil(bits_Q/8)   
    # How many bytes I need to store a single element in q:
    q = 2 * params['N']
    bits_q = math.ceil((math.log2(q)))
    bytes_q = math.ceil(bits_q/8)  

    # Size of a RLWE the ciphertext: 
    ct_size_bytes = 2 * params['N'] * bytes_Q
    ct_size_kilo_bytes = ct_size_bytes/1000.0
    print("- RLWE Ciphertext size in KB: " + str(ct_size_kilo_bytes))
    #Size of a LWE ciphertexts:
    ct_LWE_size_bytes = (params['n']+1) * bytes_q
    ct_LWE_size_kilo_bytes = ct_LWE_size_bytes/1000.0
    print("- LWE Ciphertext size in KB: " + str(ct_LWE_size_kilo_bytes))

    ct_KSK_size_bytes = (params['n']+1) * bytes_Q
    ct_KSK_size_kilo_bytes = ct_KSK_size_bytes/1000.0
    ks_key_size_kilo_bytes = params['N'] *  params['ell_ksk'] * ct_KSK_size_kilo_bytes
    ks_key_size_mega_bytes = ks_key_size_kilo_bytes/1000.0
    print("- Key Switching Key is MB: " + str(ks_key_size_mega_bytes))
 
    # Size of a NTRU (Power of two) ciphertext: 
    ct_size_bytes_Q = 2 * params['N'] * bytes_Q 
    ct_size_kilo_bytes_Q = ct_size_bytes_Q/1000.0
    # Key Switching key (N * ell_keySwitchKey * LWE sample): 
    # Size of the blind rotation keys n *   (NTRU samples)  
    br_key_size_kilo_bytes = params['n'] * params['ell_RGSW'] * ct_size_kilo_bytes_Q
    br_key_size_mega_bytes = br_key_size_kilo_bytes/1000.0 
    print("- Blind Rotation Key size in MB: " + str(br_key_size_mega_bytes))

    # Sum of all keys
    sum_of_all_key_mega_bytes = br_key_size_mega_bytes + ks_key_size_mega_bytes

    if params['version'] == 'simul':
        masking_key_kilo_bytes = params['ell_masking'] * ct_size_kilo_bytes_Q
        masking_key_mega_bytes = masking_key_kilo_bytes/1000.0 
        print("- Masking Key size in MB: " + str(masking_key_mega_bytes))
        sum_of_all_key_mega_bytes = sum_of_all_key_mega_bytes + masking_key_mega_bytes
    elif params['version'] == 'flood':
        masking_key_kilo_bytes = params['ell_masking'] * ct_size_kilo_bytes_Q
        masking_key_mega_bytes = masking_key_kilo_bytes/1000.0 
        print("- Masking Key size in MB: " + str(masking_key_mega_bytes))
        sum_of_all_key_mega_bytes = sum_of_all_key_mega_bytes + masking_key_mega_bytes
    
    sum_of_all_key_giga_bytes = sum_of_all_key_mega_bytes/1000.0
    print("- Sum of all keys in MB: " + str(sum_of_all_key_mega_bytes))
    print("- Sum of all keys in GB: " + str(sum_of_all_key_giga_bytes)) 





def TFHE_Simplified(params):  
    ell = math.ceil(math.log(params['Q'], params['basis_RGSW']))
    if(params['basis_RGSW']**ell == params['Q']):
        print("RGSW_basis**(math.ceil(math.log(Q, base))) == Q: ")
    else:
        print("RGSW_basis**(math.ceil(math.log(Q, base))) != Q: ")

    # Blind rotation
    error_params = blind_rotation_error(params)

    # Keyswitching
    error_params = key_switching_error(error_params, params['N'], params['stddev_ksk'], params['basis_ksk'], params['ell_ksk'])
    print("Correctness after blind rotation:")
    print_correctness(error_params, params['Q'], params['t'])
  
    # Mode switch from Q to 2N 
    #TODO: Actually for circuit privacy I need to check the worst case: i.e. instead of variance_s it should be computed via a bound actually
    error_params_mod_switch = modulus_reduction(error_params, params['Q'], 2 * params['N'], params['hamming_weigth_s'], params['variance_vec_s'], params['expected_vec_s'])
    print("Correctness before blind rotation:")
    print_correctness(error_params_mod_switch, 2 * params['N'], params['t']) 
    contribution_of_the_secret_key(error_params_mod_switch['stddev'], params['hamming_weigth_s'], params['variance_vec_s'], params['expected_vec_s'])
 
    print("")
    print("--- Size ================") 
    print_ct_and_key_sizes(params)


def number_of_washing_machine_cycles(stat_dist, sec_par):
    # Remind that secpar is given s.t. we have 2**(secpar) of security 
    return math.log(2**(sec_par), stat_dist)

def set_bound_on_Gaussian(secpar):
    if secpar== 80:
        return 10.285
    elif secpar == 70:
        return 9.595
    elif secpar == 60:
        return 8.857 
    elif secpar == 50:
        return 8.042
    elif secpar == 40:
        return 7.19  
    elif secpar == 30:
        return 6.13 
    elif secpar == 20:
        return 4.99 
    elif secpar == 10:
        return 3.3
    elif secpar == 5:
        return 2.2
    else:
        raise Exception("We don't have a pre-estimated bound_on_Gaussian for cp_secpar " + str(params['cp_secpar']))

def tfhe_10(): 
    print("========= TFHE 10 ============")
    params = {} 
    params['t'] = 11
    #Ring Modulus  
    m = 48
    params['Q'] = 2**(m) - 16383
    # RLWE degree
    params['N'] = 2**(11) 
    params['basis_RGSW'] = 2**25
    params['ell_RGSW'] = math.ceil(math.log(params['Q'], params['basis_RGSW'])) 
    print('ell_RGSW: ' + str(params['ell_RGSW']))
    params['stddev_br'] = 3.2
    # LWE dimension
    params['n'] = 912
    params['basis_ksk'] = 2**(7)
    params['ell_ksk'] = math.ceil(math.log(params['Q'], params['basis_ksk']))
    print('ell_ksk: ' + str(params['ell_ksk']))
    params['stddev_ksk'] = 2**26
    # Variance of the LWE secret key.
    # This is used in mod switching from Q to 2N   (the data below is for binary secret keys)
    params['variance_vec_s'] = 1/3.0
    params['expected_vec_s'] = 0.5
    # At default hamming weigth is n
    params['hamming_weigth_s'] = params['n'] 
    
    #version = 'any'
    params['version'] = 'det'
    #params['version'] = 'simul'
    #params['version'] = 'flood'

    #if(params['version'] == 'det'):  
    # If signed decomposition - we assume to have uniform distribution in (-basis_RGSW/2, basis_RGSW/2]
    # Actually no new parameters are needed to be defined for this case 
    if(params['version'] == 'simul'):
        # If simul Gaussian Sampling (secpar is the statistical distance from uniform for the regularity lemma) 
        params['cp_secpar'] = 80
        params['basis_masking'] = 2**3
        params['ell_masking'] = math.ceil(math.log(params['Q'], params['basis_masking']))
        params['stddev_masking'] = params['stddev_br'] 
        params['bound_on_Gaussian'] = set_bound_on_Gaussian(params['cp_secpar'])
    elif(params['version'] == 'flood'):
        # Take security parameter, + masking parameters (also neededm but we only care whether the a parameter is close to uniform from the LHL right?)
        params['cp_secpar'] = 80
        params['washing_cycles'] = 16
        params['basis_masking'] = 2**3
        params['ell_masking'] = math.ceil(math.log(params['Q'], params['basis_masking']))
        params['stddev_masking'] = params['stddev_br']
        #params['bound_on_Gaussian'] = set_bound_on_Gaussian(params['cp_secpar'])
        params['bound_on_Gaussian'] = set_bound_on_Gaussian(params['cp_secpar'])

    # Here I set a precalculated constant bounding_constant, s.t. Pr[X < C * \sigma] = 2^(-secpar)
    # The value of bounding_constant is precalculated using the function bound_on_Gaussian 
    

    TFHE_Simplified(params) 


#bound_on_Gaussian(1, 8.042 * 1)


tfhe_10() 