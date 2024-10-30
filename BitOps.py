def sum_half_adder(a, b):
    return "-(a) * (b) + (a) * -(b)".replace('a', a).replace('b', b)

def carry_half_adder(a, b):
    return "(a) * (b)".replace('a', a).replace('b', b)

def xnor(a, b):
    return "1 + -(a) + -(b) + 2((a) * (b))".replace('a', a).replace('b', b)


def half_adder_array(bit_list):
    for bit in range(len(bit_list) - 1, 0, -1):
        x = bit_list[bit]
        y = bit_list[bit - 1]
        bit_list[bit] = carry_half_adder(x, y)
        bit_list[bit - 1] = sum_half_adder(x, y)
    return bit_list


def sum_of_directions_helper(dir_func, amino_start, amino_end):
    counter = amino_end - amino_start
    bit_list = [None] * counter
    for i in range(counter):
        bit_list[i] = dir_func(amino_start + i)

    sum_bit_list = []

    for branch in range(counter - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list


def sum_of_directions_plus_one_helper(dir_func, amino_start, amino_end):
    counter = amino_end - amino_start + 1
    bit_list = [None] * counter
    for i in range(1, counter):
        bit_list[i] = dir_func(amino_start + i)
    bit_list[0] = '1'

    sum_bit_list = []

    for branch in range(counter - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list

