from pyqubo import *
from pyqubo import Binary, And, Xor, Not


def initialize_q_vars(num_amino, qubits_per_amino):
    # Initialize the q variables as Boolean symbols
    q_vars = {}
    for t in range(num_amino):
        for q in range(qubits_per_amino):
            q_vars[(t, q)] = Binary(f'q_{t}{chr(ord('a') + q)}')
    return q_vars

def half_adder_array(bit_list):
    for bit in range(len(bit_list) - 1, 0, -1):
        x = bit_list[bit]
        y = bit_list[bit - 1]
        bit_list[bit] = And(x, y)
        bit_list[bit - 1] = Xor(x, y)
    return bit_list

def half_adder_loader(bit_list):
    sum_bit_list = []
    size = len(bit_list)

    for _ in range(size - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list

def sum_of_directions(dir_func, amino_start, amino_end):
    bit_list = []
    for i in range(amino_end - amino_start):
        bit_list.append(dir_func(amino_start + i))

    return half_adder_loader(bit_list)


def sum_of_directions_plus_one(dir_func, amino_start, amino_end):
    bit_list = [Num(1)]
    for i in range(amino_end - amino_start):
        bit_list.append(dir_func(amino_start + i))

    return half_adder_loader(bit_list)


def sum_of_y(dir_func1, dir_func2, amino_start, amino_end):
    bit_list = []
    for i in range(amino_end - amino_start):
        bit_list.append(dir_func1(amino_start + i))
        bit_list.append(dir_func2(amino_start + i))

    return half_adder_loader(bit_list)


def sum_of_y_plus_one(dir_func1, dir_func2, amino_start, amino_end):
    bit_list = [Num(1)]
    for i in range(amino_end - amino_start):
        bit_list.append(dir_func1(amino_start + i))
        bit_list.append(dir_func2(amino_start + i))

    return half_adder_loader(bit_list)


def sum_of_y_plus_two(dir_func1, dir_func2, amino_start, amino_end):
    bit_list = [Num(1), Num(1)]
    for i in range(amino_end - amino_start):
        bit_list.append(dir_func1(amino_start + i))
        bit_list.append(dir_func2(amino_start + i))

    return half_adder_loader(bit_list)


def Xnor(a, b):
    return 1 - a - b + 2 * a * b
