import math
import random

from finites_fields import get_primitive_polynomial, build_logarithmic_table, \
    get_cyclotomic_cosets, multiply_polynomials, \
    get_polynomial_from_roots, divide_polynomials, \
    get_positions_of_binary_ones, polynomial_of_argument_to_power


class BCH(object):
    def __init__(self, p, n):
        """
        Constructs a BCH code with the
        specified parameters.
        :param p: a parameter of a
        binary symmetric channel,
        the error probability.
        :param n: length of a code.
        """
        t, n, k, power = initiate(p, n)
        self.t = t
        self.n = n
        self.k = k
        self.power = power
        self.primitive_polynomial = get_primitive_polynomial(power, 1)
        self.cyclotomic_cosets = get_cyclotomic_cosets(power)
        self.logarithmic_table = build_logarithmic_table(power, self.primitive_polynomial)
        self.generator_polynomial = calculate_generator_polynomial(self.primitive_polynomial, self.cyclotomic_cosets,
                                                                   self.logarithmic_table, self.power, self.t)
        print("---------------------------------------")
        print("t: {0}, n: {1}, k: {2}".format(t, n, k))
        print("Generator polynomial: {0:b}".format(self.generator_polynomial))


def calculate_generator_polynomial(primitive_polynomial, cyclotomic_cosets, logarithmic_table, power, t):
    """
    Calculates a generator polynomial of
    a particular BCH code.
    :param primitive_polynomial: a primitive
    polynomial of a particular Galois field.
    :return: a generator polynomial which is
    a product of several polynomails including
    a primitive polynomial obtained from
    cyclotomic cosets.
    """
    generator_polynomial = primitive_polynomial
    for i in range(1, t):
        generator_polynomial = multiply_polynomials(generator_polynomial,
                                                    get_polynomial_from_roots(cyclotomic_cosets[i], power,
                                                                              logarithmic_table))
    return generator_polynomial


def encode(generator_polynomial, message, power, t):
    """
    Encodes a message by shifting
    it to the highest power and
    adding a divided by a generator
    polynomial message to it.
    """
    message = message << power * t
    return message ^ divide_polynomials(message, generator_polynomial)[1]


def get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, power, t):
    """
    Calculates syndromes based on a particular
    received message, in number of 2 * t.
    :param primitive_polynomial: a primitive
    polynomial, primitive in sense of the power given GF(2).
    :param received_message: a message which
    was received by the decoder.
    :param cyclotomic_cosets: cyclotomic cosets
    for building minimal polynomials.
    :param logarithmic_table: a convenient form
    of the field elements for multiplication.
    :param power: the power in GF(2).
    :param t: a number of errors to correct.
    :return: a list of powers of a primitive
    element a as shortcuts for polynomials.
    """
    length = t * 2
    syndromes = [0] * length
    flipped_logarithmic_table = dict((v, k) for k, v in logarithmic_table.items())
    for i in cyclotomic_cosets:
        for position in get_positions_of_binary_ones(number=i):
            if position - 1 < length:
                syndrome_polynomial = divide_polynomials(received_message,
                                                         get_polynomial_from_roots(i, power, logarithmic_table))[1]
                syndrome_polynomial_of_argument_to_power = polynomial_of_argument_to_power(syndrome_polynomial,
                                                                                           position)
                if len(bin(syndrome_polynomial_of_argument_to_power)) >= len(bin(primitive_polynomial)):
                    syndrome_polynomial_of_argument_to_power = divide_polynomials(syndrome_polynomial_of_argument_to_power,
                                                                                  primitive_polynomial)[1]
                syndrome = flipped_logarithmic_table[syndrome_polynomial_of_argument_to_power]
                syndromes[position - 1] = syndrome
    return syndromes


def berlekamp_massey_decode(syndromes, logarithmic_table, power, t):
    """
    Calculates an error locator polynomial using
    the Berlekamp-Massey algorithm.
    :param syndromes: a previously calculated
    array of syndromes corresponding to the
    received distorted message.
    :param logarithmic_table: a mapping of a
    power to a corresponding element of a particular
    Galois field.
    :param power: the power in a Galois field GF(2).
    :param t: a number of error to be corrected.
    :return: an array representing an error locator
    polynomial.
    """
    flipped_logarithmic_table = dict((v, k) for k, v in logarithmic_table.items())
    number_of_elements = 2 ** power - 1
    capital_lambdas = [[0] * (2 * t + 1)] * (2 * t + 1)
    capital_lambdas[-1][0] = 1
    capital_lambdas[0][0] = 1
    capital_ls = [0] * number_of_elements
    deltas = [0] * (2 * t + 1)
    deltas[0] = 1
    for k in range(1, 2 * t + 1):
        # calculates ∆
        for i in range(capital_ls[k - 1] + 1):
            if (k - 1 - i < len(syndromes)) and (syndromes[k - 1 - i] >= 0) and (capital_lambdas[k - 1][i] > 0):
                deltas[k] ^= logarithmic_table[
                    (flipped_logarithmic_table[capital_lambdas[k - 1][i]] +
                     syndromes[k - 1 - i]) % number_of_elements]
        # ends
        m_bag = []
        for i in range(k):
            if capital_ls[i] == capital_ls[k - 1]:
                m_bag.append(i)
        m = sorted(m_bag)[0]
        if deltas[k] == 0:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            capital_ls[k] = capital_ls[k - 1]
        else:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            multiplied_c_l = [0] * (len(capital_lambdas[m - 1]) + k - m)
            for i in range(len(capital_lambdas[m - 1])):
                if capital_lambdas[m - 1][i] > 0:
                    multiplied_c_l[i + k - m] = logarithmic_table[
                        (flipped_logarithmic_table[capital_lambdas[m - 1][i]] +
                         flipped_logarithmic_table[deltas[k]] +
                         number_of_elements - flipped_logarithmic_table[deltas[m]]
                         ) % number_of_elements
                        ]
            who = len(capital_lambdas[k]) - len(multiplied_c_l)
            if who <= 0:
                for i in range(who):
                    capital_lambdas[k][i].append(0)
            else:
                for i in range(who):
                    multiplied_c_l.append(0)
            for i in range(len(capital_lambdas[k])):
                capital_lambdas[k][i] = capital_lambdas[k][i] ^ multiplied_c_l[i]
            # ends
            capital_ls[k] = max(capital_ls[k - 1], capital_ls[m - 1] + k - m)
    result = []
    for i in capital_lambdas[2 * t]:
        result.append(flipped_logarithmic_table[i])
    return result


def find_roots_of_sigma(sigma, power, logarithmic_table):
    """
    Tries all the elements of a field to solve an
    error polynomial equls to zero.
    :param sigma: a sigma is an array representing
    an error locator polynomial.
    :param power: the power in a Galois field GF(2).
    :param logarithmic_table: a mapping of a power
    to a corresponding element of a particular
    Galois field.
    :return: an array of roots of a sigma.
    """

    def get_order_of_sigma():
        """
        Returns an order of an
        arror locator polynomial.
        """
        counter = 0
        for i in reversed(sigma):
            if i <= -1:
                counter += 1
                continue
            return len(sigma) - 1 - counter
        return 0

    roots = []
    for candidate in range(2 ** power - 1):
        result = logarithmic_table[sigma[0]]
        for polynomial_power in range(1, get_order_of_sigma() + 1):
            if sigma[polynomial_power] >= 0:
                result ^= logarithmic_table[(sigma[polynomial_power] + candidate * polynomial_power) % (2 ** power - 1)]
        if result == 0:
            roots.append(candidate)
    return roots


def decode(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, power, t):
    """
    Performs decoding of a received message.
    :param primitive_polynomial: a primitive
    polynomial over a Galois field.
    :return: a decoded message.
    """

    def get_error_positions(_roots, _power):
        """
        Flips the established roots of
        an error locator polynomial to
        the positions of errors.
        """
        positions = []
        for root in _roots:
            positions.append((2 ** _power - 1 - root) % (2 ** _power - 1))
        return positions

    syndromes = get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets,
                              logarithmic_table, power, t)
    sigma = berlekamp_massey_decode(syndromes, logarithmic_table, power, t)
    roots = find_roots_of_sigma(sigma, power, logarithmic_table)
    error_positions = get_error_positions(roots, power)
    for position in error_positions:
        received_message ^= 1 << position
    received_message >>= power * t
    return received_message


def get_random_number_of_hamming_weight(length, weight):
    """
    Returns a random number of a particular Hamming
    weight.
    :param length: a length of a number.
    :param weight: a weight of a number.
    :return: the number meets the requirements of
    length and Hamming weight.
    :raises: ValueError if the weight is greater
    than the length.
    """
    if weight > length:
        raise ValueError('The weight shouldn\'t be greater than the length: {0} > {1}'.format(weight, length))
    i = 0
    result = 0
    while True:
        if i == weight:
            return result
        shift = random.randrange(length)
        power_of_two = 1 << shift
        if power_of_two & result == power_of_two:
            continue
        result |= power_of_two
        i += 1


def translate_message_to_bits_and_split_on_blocks_of_length_k(message, k):
    """
    Translates a message to bit representation and splits bits
    on blocks of length k.
    """
    bits = int.from_bytes(message.encode('utf-8', 'custompassword'), 'big')
    blocks = []
    while bits != 0:
        block = bits & (2 ** k - 1)
        blocks.append(block)
        bits >>= k
    return blocks


def translate_bits_to_message_and_glue_blocks_of_length_k(blocks, k):
    """
    Translates an array of bit vectors to string.
    :param blocks: an array of bit vectors.
    :param k: length of a bit vector.
    :return: a string of united translated vectors.
    """
    result = 0
    for block in reversed(blocks):
        result <<= k
        result ^= block
    return result, result.to_bytes((result.bit_length() + 7) // 8, 'big').decode('utf-8', 'custompassword') or '\0'


def initiate(p, n):
    """
    Initiates a process of data transmission.
    :return: a number of errors which can be
    corrected by a BCH code, p, a corrected n,
    length of a message as an input to the channel, power
    and the power of GF(2) which corresponds to n.
    """

    def is_power_of_two(num):
        """
        Checks the power of 2.
        :return: true if the number has only one 1
        and all the remaining 0, false otherwise.
        """
        counter = 0
        for i in (range(len(bin(num)) - 2)):
            if num & 1:
                counter += 1
            num >>= 1
            i += 1
        return counter == 1

    if (p > 1 / 3) | (p < 0):
        raise ValueError('The parameter p should be 0 <= p <= 1/3, p is {0}'.format(p))
    power = 0
    if not is_power_of_two(n + 1):
        power = len(bin(n)) - 3
        n = 2 ** power - 1
    else:
        power = len(bin(n + 1)) - 3
    while p * n > power - 1:
        n = (n + 1) / 2 - 1
        power -= 1
    d = 2 * math.ceil(p * n) + 1
    t = int((d - 1) / 2)
    k = 2 ** power - 1 - t * power
    return t, int(n), k, power
