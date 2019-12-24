import csv


from itertools import combinations


def get_primitive_polynomial(power, k):
    """
    Retrieves a table of primitive
    polynomials from the given in
    the config file, and returns a
    corresponding to the
    parameters polynomial.
    :return: a primitive polynomial
    in a binary representation.
    """
    if power < 1 | power > 51:
        raise ValueError('The parameter n should be 1 <= power <= 51, power is {0}'.format(power))
    if k < 1 | k > 3:
        raise ValueError('The parameter k should be 1 <= k <= 3, k is {0}'.format(k))
    dictionary = {1: 1, 2: 2, 3: 5, 4: 10}
    with open('primitive-polynomials.csv', newline='') as csvfile:
        primitive_polynomials = csv.reader(csvfile, delimiter=',')
        for row in primitive_polynomials:
            ints = list(map(lambda x: 0 if x == '' else int(x), row))
            if ints[0] == power:
                polynomial_powers = ints[dictionary[k]:dictionary[k + 1]]
                polynomial_binary = 1 << power
                polynomial_binary |= 1
                for i in polynomial_powers:
                    polynomial_binary |= 1 << i
                return polynomial_binary


def build_logarithmic_table(power, primitive_polynomial):
    """
    Builds a logarithmic table where a
    logarithm is mapped to the binary
    representation of the corresponding
    polynomial.
    :return: the logarithmic table of
    the field.
    """
    if len(bin(primitive_polynomial)) - 3 != power:
        raise ValueError('The primitive polynomial {0:b} '
                         'is not of the specified power n = {1}'.
                         format(primitive_polynomial, power))
    logarithmic_table = {-1: 0}
    for i in range(power):
        logarithmic_table[i] = 1 << i
    logarithmic_table[power] = primitive_polynomial & ((2 ** power) - 1)
    for i in range(power + 1, 2 ** power - 1):
        multiplied_by_x_polynomial = logarithmic_table[i - 1] << 1
        if multiplied_by_x_polynomial & (2 ** power):
            multiplied_by_x_polynomial ^= logarithmic_table[power]
        logarithmic_table[i] = multiplied_by_x_polynomial & ((2 ** power) - 1)
    return logarithmic_table


def multiply_polynomials(polynomial1, polynomial2):
    """
                                      power
    Multiplies two polynomials in GF(2     ).
    :param polynomial1: 1st polynomial.
    :param polynomial2: 2nd polynomial.
    :return: the product of two polynomials.
    """
    result = 0
    for i in range(len(bin(polynomial2)) - 2):
        if polynomial2 & (1 << i):
            result ^= polynomial1 << i
    return result


def divide_polynomials(polynomial1, polynomial2):
    quotient = 0
    reminder = polynomial1
    while len(bin(reminder)) >= len(bin(polynomial2)):
        shift = len(bin(reminder)) - len(bin(polynomial2))
        reminder ^= polynomial2 << shift
        quotient ^= 1 << shift
    return quotient, reminder


def polynomial_of_argument_to_power(polynomial, power):
    """
    Calculates terms of a polynomial putting
    in corresponding power of the variable
    and reduces the resulting polynomial.
    :return: a reduced polynomial of the
    specified power of the variable.
    """
    length = len(bin(polynomial)) - 2
    result = 0
    for i in range(length):
        if polynomial >> i & 1 == 1:
            result |= 1 << i * power
    return result


def get_cyclotomic_cosets(power):
    """
    Fills a list of cyclotomic cosets.
    :return: a list of cyclotomic cosets
    in the binary representation.
    """
    cyclotomic_cosets = []
    all_cyclotomic_members = 1
    i = 0
    while all_cyclotomic_members < 2 ** (2 ** power - 2) - 1:
        cyclotomic_cosets.append(0)
        k = 0
        while True:
            if not 1 & (all_cyclotomic_members >> k):
                break
            k += 1
        while True:
            k = k % (2 ** power - 1)
            if 1 & (cyclotomic_cosets[i] >> k):
                break
            cyclotomic_cosets[i] ^= 1 << k
            k *= 2
        all_cyclotomic_members ^= cyclotomic_cosets[i]
        i += 1
    return cyclotomic_cosets


def get_polynomial_from_roots(roots, power, logarithmic_table):
    """
    Performs multiplication of a
    polynomial represented by its
    roots.
    :return: a binary vector
    represents a polynomial.
    """
    if roots == 0:
        return 0
    number_of_field_elements = 2 ** power - 1
    root_array = get_positions_of_binary_ones(roots)
    polynomial = 1 << len(root_array)
    for i in range(len(root_array)):
        coefficient = 0
        for combination in combinations(root_array, i + 1):
            coefficient ^= logarithmic_table[sum(combination) % number_of_field_elements]
        addition = coefficient << len(root_array) - i - 1
        polynomial ^= addition
    return polynomial


def get_positions_of_binary_ones(number):
    """
    Gets positions of 1s in a binary
    vector.
    :param number: a vector to
    analyze.
    :return: positions of 1s.
    """
    length = len(bin(number)) - 2
    result = []
    for i in range(0, length):
        if 1 & (number >> i):
            result.append(i)
    return result
