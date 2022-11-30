import numpy as np
import itertools

# The three functions defined below create the p vector and the corresponding lambda vectors

def p_array(n):
    '''Generates a random p vector, where p represents the parameter and is the difference between the 
    probabilities of the state being in the ground or the excited state i.e. 0 or 1 state.
    The arguement is n, the number of RPQ's but total number of qubits is N = n + 1 as this includes the SPQ'''
    return np.random.random(n+1)

def uniform_p_array(n):
    '''Generates a random uniformly sampled p vector, where p represents the parameter and is the difference between the 
    probabilities of the state being in the ground or the excited state i.e. 0 or 1 state.
    The arguement is n, the number of RPQ's but total number of qubits is N = n + 1 as this includes the SPQ'''
    return np.random.random_sample(n)

def lambda0_array(array):
    '''This function takes in an array as an argument and works out the corresponding lambda 0 vector.
    So a random array is only generated at the beginning and in one instance, in the above code for p, thus making the code more efficient.
    I then assign that array to a variable and that variable goes into the rest of my functions below here
    which take an array as an argument.'''
    comps = []
    for x in array:
        comp = (1 + x) / 2
        comps.append(comp)
    return np.array(comps)

def lambda1_array(array):
    '''This function takes in an array as an argument and works out the corresponding lambda 1 vector.
    So a random array is only generated at the beginning and in one instance, in the above code for p, thus making the code more efficient.
    I then assign that array to a variable and that variable goes into the rest of my functions below here
    which take an array as an argument.'''
    comps = []
    for x in array:
        comp = (1 - x)/2
        comps.append(comp)
    return np.array(comps)

        
def string_permutations(n):
    '''This generates a string of binary digits of length n and all the unique permutations, where n is the number of RPQ's.
    Essentially generating the measurement results of the RPQ's.'''
    permutations = []
    for i in itertools.product([0,1], repeat = n):
        permutations.append(i)
    return np.array(permutations)

# The two below functions get fed a specific permutation and they calculate the numerators of the CFI formula and I will then loop over all permutations in the overall CFI function

def lambda1_productfunc(permutation, array):
    '''This lambda1 product function works out the lambda 1 product sum present in the numerator of the 
    final CFI formula but it only does so for a specifc permutation of the binary string j, which represents 
    a measurement outcome of the RPQ's'''
    product = lambda1_array(array)[0]
    for j, i, l in zip(permutation, lambda0_array(array)[1:], lambda1_array(array)[1:]):
        if j == 0:
            product *= l
        else:
            product *= i
    return product

def lambda0_productfunc(permutation, array):
    '''This lambda0 product function works out the lambda0 product sum present in the numerator of the 
    final CFI formula but it only does so for a specifc permutation of the binary string j, which represents 
    a specfic measurement outcome of the RPQ's'''
    product = lambda0_array(array)[0]
    for j, i, l in zip(permutation, lambda0_array(array)[1:], lambda1_array(array)[1:]):
        if j == 0:
            product *= i
        else:
            product *= l
    return product     
           
    
""" These product sum functions below are calculating the different products of the components of 
the sum of the kappa, eta_parallel and p vector with the power of the object depending on the state of the RPQ's measured; 
whether it is in the 0 or 1 state after a sigma_z measurement, this power is represented by the j in the 
below formula and is independent for each qubit. These below functions are calculating the products for a specific permutation
In the below function, n is the number of RPQ's"""

def product_sum1(permutation, kappa, eta_parallel, array):
    product = 1
    for j, p in zip(permutation, array[1:]):
        product *= 1 + (-1)**j * kappa + (-1)**j * eta_parallel * p
    return product

def product_sum2(permutation, kappa, eta_parallel, array):
    product = 1
    for j, p in zip(permutation, array[1:]):
        product *= 1 + (-1)**j * kappa + (-1)**(1-j) * eta_parallel * p
    return product

def product_sum3(permutation, kappa, eta_parallel, array):
    product = 1
    for j, p in zip(permutation, array[1:]):
        product *= 1 + (-1)**(1-j) * kappa + (-1)**j * eta_parallel * p
    return product

def product_sum4(permutation, kappa, eta_parallel, array):
    product = 1
    for j, p in zip(permutation, array[1:]):
        product *= 1 + (-1)**(1-j) * kappa + (-1)**(1-j) * eta_parallel * p
    return product 


# =============================================================================
# def CFI(kappa, eta_parallel, array):
#     '''This function calculates the Classical Fisher information over all permutations of the string_permutations object. 
#     The np.count_nonzero function returns the hamming weight for that specific permutation of binary digits. 
#     The functions inside the for loop depend on the permutations which we are looping over. n is the number of RPQ's.'''
#     Ap = 1 + kappa + eta_parallel
#     Am = 1 + kappa - eta_parallel
#     Bp = 1 - kappa + eta_parallel
#     Bm = 1 - kappa - eta_parallel
#     n = len(array) - 1
#     fisher = 0
#     for permutation in string_permutations(n):
#         fisher += ((2**(n+2) * (n + 1 - 2 * np.count_nonzero(permutation))**2 * (lambda1_productfunc(permutation, array) - lambda0_productfunc(permutation, array))**2) /
#                       ((lambda0_array(array)[0] * Ap + lambda1_array(array)[0] * Am) * product_sum1(permutation, kappa, eta_parallel, array) + 
#                        (lambda0_array(array)[0] * Am + lambda1_array(array)[0] * Ap) * product_sum2(permutation, kappa, eta_parallel, array) +
#                        (lambda0_array(array)[0] * Bp + lambda1_array(array)[0] * Bm) * product_sum3(permutation, kappa, eta_parallel, array) +
#                        (lambda0_array(array)[0] * Bm + lambda1_array(array)[0] * Bp) * product_sum4(permutation, kappa, eta_parallel, array)))
#     return fisher
# =============================================================================


# The below CFI formula is for when you include eta_perpendicular into the calculations and plots and if you want dephasing case or to ignore
# eta_perpendicular effects then just set it equal to 1. 
def CFI(kappa, eta_parallel, eta_perpendicular, array):
    '''This function calculates the Classical Fisher information over all permutations of the string_permutations object. 
    The np.count_nonzero function returns the hamming weight for that specific permutation of binary digits. 
    The functions inside the for loop depends on the permutations which we are looping over.'''
    Ap = 1 + kappa + eta_parallel
    Am = 1 + kappa - eta_parallel
    Bp = 1 - kappa + eta_parallel
    Bm = 1 - kappa - eta_parallel
    n = len(array) - 1
    fisher = 0
    for permutation in string_permutations(n):
        fisher += ((2**(n+2) * eta_perpendicular**(2*(n+1)) * (n + 1 - 2 * np.count_nonzero(permutation))**2 * (lambda1_productfunc(permutation, array) - lambda0_productfunc(permutation, array))**2) /
                      ((lambda0_array(array)[0] * Ap + lambda1_array(array)[0] * Am) * product_sum1(permutation, kappa, eta_parallel, array) + 
                       (lambda0_array(array)[0] * Am + lambda1_array(array)[0] * Ap) * product_sum2(permutation, kappa, eta_parallel, array) +
                       (lambda0_array(array)[0] * Bp + lambda1_array(array)[0] * Bm) * product_sum3(permutation, kappa, eta_parallel, array) +
                       (lambda0_array(array)[0] * Bm + lambda1_array(array)[0] * Bp) * product_sum4(permutation, kappa, eta_parallel, array)))
    return fisher


# =============================================================================
# need to confirm the Bp and Bm swap on the last two product functions in the code above. In my written calculations, I have
# Bp * lambda0 and Bp * lambda1 for the product_sum3 function and the opposite of that for the product_sum4 function. 
# Agnez's is opposite to this. this was resolved as we defined our Bp and Bm opposite to each other.
# =============================================================================