"""
Direction & Movement Vector & Integer Approximation & Encoding
North (N)           & $[0,, 1,, 0]$                       & $[0,, 2,, 0]$   & $111$ 
South (S)           & $[0,,-1,, 0]$                      & $[0,,-2,, 0]$  & $000$ 
Northeast (NE)      & $left[ tfrac{sqrt{3}}{2},, tfrac{1}{2},, 0 right]$ & $[1,, 1,, 0]$   & $011$ 
Southwest (SW)      & $left[ -tfrac{sqrt{3}}{2},, -tfrac{1}{2},, 0 right]$ & $[-1,,-1,, 0]$  & $100$ 
Northwest (NW)      & $left[ -tfrac{sqrt{3}}{2},, tfrac{1}{2},, 0 right]$ & $[-1,, 1,, 0]$  & $101$ 
Southeast (SE)      & $left[ tfrac{sqrt{3}}{2},, -tfrac{1}{2},, 0 right]$ & $[1,,-1,, 0]$   & $010$ 
Up (U)              & $[0,, 0,, 1]$                       & $[0,, 0,, 1]$   & $001$ 
Down (D)            & $[0,, 0,,-1]$                      & $[0,, 0,,-1]$  & $110$ 
"""

import math
from BitOps import sum_half_adder, carry_half_adder, xnor, half_adder_array, sum_of_directions_plus_one_helper, \
    sum_of_directions_helper
from Energy import get_energy_matrix


def create_energy_function(sequence, energy_model):
    # Create the energy function
    # Return the energy function
    num_amino = len(sequence)

    overlap = create_overlap_constraint(num_amino)
    backup = create_back_constraint(num_amino)

    energy_matrix = get_energy_matrix(sequence, energy_model)
    interactions = create_interactions(sequence, energy_matrix)

    return ((f"OVERLAP_PENALTY\t({overlap}) + \n"
            f"BACKUP_PENATLY\t({backup}) + \n"
            f"{interactions}")
            .replace('q^0_1', '0').replace('q^0_2', '0'))

def dx_plus(t):
    return "-q^t_1 * q^t_2".replace('t', str(t))

def dx_minus(t):
    return "q^t_1 * -q^t_2".replace('t', str(t))

def dy_plus(t):
    return "(q^t_1 * q^t_2 * q^t_3) + (-q^t_1 * q^t_2 * q^t_3) + (q^t_1 * -q^t_2 * q^t_3)".replace('t', str(t))

def dy_minus(t):
    return '(q^t_1 * q^t_2 * q^t_3) + (-q^t_1 * q^t_2 * q^t_3) + (q^t_1 * -q^t_2 * q^t_3)'.replace('t', str(t))

def dy_plus_plus(t):
    return "q^t_1 * q^t_2 * q^t_3".replace('t', str(t))

def dy_minus_minus(t):
    return '-q^t_1 * -q^t_2 * -q^t_3'.replace('t', str(t))

def dz_plus(t):
    return "-q^t_1 * -q^t_2 * q^t_3".replace('t', str(t))

def dz_minus(t):
    return 'q^t_1 * q^t_2 * -q^t_3'.replace('t', str(t))


def sum_of_y_plus_two(amino_start, amino_end):
    counter = amino_end - amino_start
    bit_list = [None] * ((counter * 2) + 2)
    for i in range(2, counter):
        bit_list[i] = dy_plus(amino_start + i)
        bit_list[i + amino_end - amino_start] = dy_plus_plus(amino_start + i)

    bit_list[0] = '1'
    bit_list[1] = '1'

    sum_bit_list = []

    for branch in range((counter * 2) - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list


def sum_of_y_minus_two(amino_start, amino_end):
    counter = amino_end - amino_start
    bit_list = [None] * ((counter * 2) + 2)
    for i in range(2, counter):
        bit_list[i] = dy_minus(amino_start + i)
        bit_list[i + counter] = dy_minus_minus(amino_start + i)

    bit_list[0] = '1'
    bit_list[1] = '1'

    sum_bit_list = []

    for branch in range((counter * 2) - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list


def sum_of_directions_plus_one(direction, amino_start, amino_end):
    # Map the direction string to the corresponding function
    direction_funcs = {
        'dx_plus': dx_plus,
        'dx_minus': dx_minus,
        'dy_plus': dy_plus,
        'dy_minus': dy_minus,
        'dz_plus': dz_plus,
        'dz_minus': dz_minus
    }
    dir_func = direction_funcs[direction]

    if direction == 'dy_plus':
        return sum_of_y_plus_two(amino_start, amino_end)
    elif direction == 'dy_minus':
        return sum_of_y_minus_two(amino_start, amino_end)

    return sum_of_directions_plus_one_helper(dir_func, amino_start, amino_end)


def sum_of_y_plus(amino_start, amino_end):
    counter = amino_end - amino_start
    bit_list = [None] * (counter * 2)
    for i in range(counter):
        bit_list[i] = dy_plus(amino_start + i)
        bit_list[i + counter] = dy_plus_plus(amino_start + i)

    sum_bit_list = []

    for branch in range((counter * 2) - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list


def sum_of_y_minus(amino_start, amino_end):
    counter = amino_end - amino_start
    bit_list = [None] * (counter * 2)
    for i in range(counter):
        bit_list[i] = dy_minus(amino_start + i)
        bit_list[i + counter] = dy_minus_minus(amino_start + i)

    sum_bit_list = []

    for branch in range((counter * 2) - 1):
        bit_list = half_adder_array(bit_list)
        sum_bit_list.append(bit_list[0])
        bit_list = bit_list[1:]

    sum_bit_list.append(bit_list[0])
    return sum_bit_list


def sum_of_directions(direction, amino_start, amino_end):
    # Map the direction string to the corresponding function
    direction_funcs = {
        'dx_plus': dx_plus,
        'dx_minus': dx_minus,
        'dy_plus': dy_plus,
        'dy_minus': dy_minus,
        'dy_plus_plus': dy_plus_plus,
        'dy_minus_minus': dy_minus_minus,
        'dz_plus': dz_plus,
        'dz_minus': dz_minus
    }
    dir_func = direction_funcs[direction]

    if direction == 'dy_plus':
        return sum_of_y_plus(amino_start, amino_end)
    elif direction == 'dy_minus':
        return sum_of_y_minus(amino_start, amino_end)

    return sum_of_directions_helper(dir_func, amino_start, amino_end)


def create_overlap_constraint(num_amino):
    # Create the overlap constraint
    # Return the overlap constraint
    # returns a string representing an overlap constraint
    overlap = ""

    for amino1 in range(num_amino):
        for amino2 in range(num_amino):
            if amino1 + 1 >= amino2:
                continue

            dx_plus_sum = sum_of_directions('dx_plus', amino1, amino2)
            dx_minus_sum = sum_of_directions('dx_minus', amino1, amino2)
            dy_plus_sum = sum_of_directions('dy_plus', amino1, amino2)
            dy_minus_sum = sum_of_directions('dy_minus', amino1, amino2)
            dz_plus_sum = sum_of_directions('dz_plus', amino1, amino2)
            dz_minus_sum = sum_of_directions('dz_minus', amino1, amino2)

            for bit in range(math.ceil(math.log2(amino2 - amino1))):
                overlap += xnor(dx_plus_sum[bit], dx_minus_sum[bit]) + " + "
                overlap += xnor(dz_plus_sum[bit], dz_minus_sum[bit]) + " + "

            for bit in range(math.ceil(math.log2((amino2 - amino1 + 2))) ):
                overlap += xnor(dy_plus_sum[bit], dy_minus_sum[bit]) + " + "

    return overlap[:-3]


def create_interactions(sequence, energy_matrix):
    # Create the interactions
    # Return the interactions
    # returns a string representing the interactions
    interactions = ""

    for amino1 in range(len(sequence)):
        for amino2 in range(len(sequence)):
            if amino1 + 1 >= amino2:
                continue

            energy_value = energy_matrix[amino1, amino2]
            if energy_value == 0:
                continue

            interactions += f'INTERACTION({amino1}-{amino2})\t' + str(energy_value) + ' * (' + adjacency_indicator(amino1, amino2) + ')\n'

    return interactions


def create_back_constraint(num_amino):
    # Create the back constraint
    # Return the back constraint string
    # returns a string representing the back constraint

    def backwards_north(t):
        return "(q^a_1 * q^a_2 * q^a_3 * -q^b_1 * -q^b_2 * -q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_south(t):
        return "(-q^a_1 * -q^a_2 * -q^a_3 * q^b_1 * q^b_2 * q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_northeast(t):
        return "(-q^a_1 * q^a_2 * q^a_3 * q^b_1 * -q^b_2 * -q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_southwest(t):
        return "(q^a_1 * -q^a_2 * -q^a_3 * -q^b_1 * q^b_2 * q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_northwest(t):
        return "(q^a_1 * -q^a_2 * q^a_3 * -q^b_1 * q^b_2 * -q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_southeast(t):
        return "(-q^a_1 * q^a_2 * -q^a_3 * q^b_1 * -q^b_2 * q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_up(t):
        return "(-q^a_1 * -q^a_2 * q^a_3 * q^b_1 * q^b_2 * -q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_down(t):
        return "(q^a_1 * q^a_2 * -q^a_3 * -q^b_1 * -q^b_2 * q^b_3)".replace('a', str(t)).replace('b', str(t + 1))

    back = ""
    for turn in range(num_amino - 2):
        back += backwards_north(turn) + " + "
        back += backwards_south(turn) + " + "
        back += backwards_northeast(turn) + " + "
        back += backwards_southwest(turn) + " + "
        back += backwards_northwest(turn) + " + "
        back += backwards_southeast(turn) + " + "
        back += backwards_up(turn) + " + "
        back += backwards_down(turn) + " + "

    return back[:-3]


def adjacency_indicator(amino1, amino2):
    # amino1 is smaller than amino2
    # returns a string representing the adjacency indicator
    dx_plus_sum = sum_of_directions('dx_plus', amino1, amino2)
    dx_minus_sum = sum_of_directions('dx_minus', amino1, amino2)
    dy_plus_sum = sum_of_directions('dy_plus', amino1, amino2)
    dy_minus_sum = sum_of_directions('dy_minus', amino1, amino2)
    dz_plus_sum = sum_of_directions('dz_plus', amino1, amino2)
    dz_minus_sum = sum_of_directions('dz_minus', amino1, amino2)

    dx_plus_sum_plus_one = sum_of_directions_plus_one('dx_plus', amino1, amino2)
    dx_minus_sum_plus_one = sum_of_directions_plus_one('dx_minus', amino1, amino2)
    dy_plus_sum_plus_two = sum_of_directions_plus_one('dy_plus', amino1, amino2)
    dy_minus_sum_plus_two = sum_of_directions_plus_one('dy_minus', amino1, amino2)
    dz_plus_sum_plus_one = sum_of_directions_plus_one('dz_plus', amino1, amino2)
    dz_minus_sum_plus_one = sum_of_directions_plus_one('dz_minus', amino1, amino2)

    x_equal = ''
    y_equal = ''
    z_equal = ''
    x_plus_one = ''
    x_minus_one = ''
    y_plus_two = ''
    y_minus_two = ''
    z_plus_one = ''
    z_minus_one = ''

    for bit in range(math.ceil(math.log2(amino2 - amino1))):
        x_equal += xnor(dx_plus_sum[bit], dx_minus_sum[bit]) + " + "
        y_equal += xnor(dy_plus_sum[bit], dy_minus_sum[bit]) + " + "
        z_equal += xnor(dz_plus_sum[bit], dz_minus_sum[bit]) + " + "
        x_plus_one += xnor(dx_plus_sum_plus_one[bit], dx_minus_sum[bit]) + " + "
        x_minus_one += xnor(dx_plus_sum[bit], dx_minus_sum_plus_one[bit]) + " + "
        z_plus_one += xnor(dz_plus_sum_plus_one[bit], dz_minus_sum[bit]) + " + "
        z_minus_one = xnor(dz_plus_sum[bit], dz_minus_sum_plus_one[bit]) + " + "

    for bit in range(math.ceil(math.log2((amino2 - amino1) + 2))):
        y_plus_two += xnor(dy_plus_sum_plus_two[bit], dy_minus_sum[bit]) + " + "
        y_minus_two += xnor(dy_plus_sum[bit], dy_minus_sum_plus_two[bit]) + " + "

    # x equal and z equal while y is +2 or -2
    x_equal_z_equal_y_offset = '((' + x_equal[:-3] + ') * ((' + z_equal[:-3] + ') * ((' + y_plus_two[:-3] + ') + (' + y_minus_two[:-3] + '))))'

    # x equal and y equal while z is +1 or -1
    xy_equal_z_offset = '((' + x_equal[:-3] + ') * (' + y_equal[:-3] + ') * ((' + z_plus_one[:-3] + ') + (' + z_minus_one[:-3] + ')))'

    # z equal and while x is +1 or -1 and y is +1 or -1
    z_equal_xy_offset = '((' + z_equal[:-3] + ') * ((' + x_plus_one[:-3] + ') + (' + x_minus_one[:-3] + ') * ((' + y_plus_two[:-3] + ') + (' + y_minus_two[:-3] + ')))'

    return x_equal_z_equal_y_offset + ' + ' + xy_equal_z_offset + ' + ' + z_equal_xy_offset
