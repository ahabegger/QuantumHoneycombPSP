"""
Direction & Movement Vector & Integer Approximation & Encoding $(q^t_1 q^t_2 q^t_3 q^t_4)$ 
Northeast (NE)          & $[ {sqrt{2}}{2},, {sqrt{2}}{2},, 0 ]$       & $[1,,1,,0]$     & $1011$ 
Northwest (NW)          & $[ -{sqrt{2}}{2},, {sqrt{2}}{2},, 0 ]$      & $[-1,,1,,0]$    & $1111$ 
Southeast (SE)          & $[ {sqrt{2}}{2},, -{sqrt{2}}{2},, 0 ]$      & $[1,,-1,,0]$    & $1010$ 
Southwest (SW)          & $[ -{sqrt{2}}{2},, -{sqrt{2}}{2},, 0 ]$     & $[-1,,-1,,0]$   & $1110$ 
Up-North (UN)           & $[ 0,, {sqrt{2}}{2},, {sqrt{2}}{2} ]$       & $[0,,1,,1]$     & $0111$ 
Up-South (US)           & $[ 0,,-{sqrt{2}}{2},, {sqrt{2}}{2} ]$       & $[0,,-1,,1]$    & $0101$ 
Down-North (DN)         & $[ 0,, {sqrt{2}}{2},, -{sqrt{2}}{2} ]$      & $[0,,1,,-1]$    & $0110$ 
Down-South (DS)         & $[ 0,,-{sqrt{2}}{2},, -{sqrt{2}}{2} ]$      & $[0,,-1,,-1]$   & $0100$ 
Up-East (UE)            & $[ {sqrt{2}}{2},, 0,, {sqrt{2}}{2} ]$       & $[1,,0,,1]$     & $1001$ 
Up-West (UW)            & $[ -{sqrt{2}}{2},, 0,, {sqrt{2}}{2} ]$      & $[-1,,0,,1]$    & $1101$ 
Down-East (DE)          & $[ {sqrt{2}}{2},, 0,,-{sqrt{2}}{2} ]$       & $[1,,0,,-1]$    & $1000$ 
Down-West (DW)          & $[ -{sqrt{2}}{2},, 0,,-{sqrt{2}}{2} ]$      & $[-1,,0,,-1]$   & $1100$ 
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
    redun = create_redundancy_constraint(num_amino)

    energy_matrix = get_energy_matrix(sequence, energy_model)
    interactions = create_interactions(sequence, energy_matrix)

    return ((f"REDUNDUNCY_PENALTY\t({redun})\n"
            f"OVERLAP_PENALTY\t({overlap})\n"
            f"BACKUP_PENATLY\t({backup})\n"
            f"{interactions}")
            .replace('q^0_1', '1').replace('q^0_2', '1')
            .replace('q^0_3', '1').replace('q^0_4', '0')
            .replace('q^1_1', '1').replace('q^1_4', '0'))

def dx_plus(t):
    return "q^t_1 * -q^t_2".replace('t', str(t))

def dx_minus(t):
    return "q^t_1 * q^t_2".replace('t', str(t))

def dy_plus(t):
    return "q^t_1 * q^t_3 * q^t_4 + -q^t_1 * q^t_2 * q^t_3".replace('t', str(t))

def dy_minus(t):
    return 'q^t_1 * q^t_3 * -q^t_4 + -q^t_1 * q^t_2 * -q^t_3'.replace('t', str(t))

def dz_plus(t):
    return "q^t_1 * -q^t_3 * q^t_4 + -q^t_1 * q^t_2 * q^t_4".replace('t', str(t))

def dz_minus(t):
    return 'q^t_1 * -q^t_3 * -q^t_4 + -q^t_1 * q^t_2 * -q^t_4'.replace('t', str(t))

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

    return sum_of_directions_plus_one_helper(dir_func, amino_start, amino_end)


def sum_of_directions(direction, amino_start, amino_end):
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

    return sum_of_directions_helper(dir_func, amino_start, amino_end)


def create_redundancy_constraint(num_amino):
    # Create the redundancy constraint
    # Return the redundancy constraint
    # returns a string representing the redundancy constraint
    redundancy = ''

    #0000 0011 0001 0010 invalid directions
    for t in range(num_amino - 1):
        redundancy += f'(-q^{t}_1 * -q^{t}_2) + '
    return redundancy[:-3]


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
                overlap += xnor(dy_plus_sum[bit], dy_minus_sum[bit]) + " + "
                overlap += xnor(dz_plus_sum[bit], dz_minus_sum[bit]) + " + "

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

    def backwards_ne(t):
        return "(q^a_1 * -q^a_2 * q^a_3 * q^a_4 * q^b_1 * q^b_2 * q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_sw(t):
        return "(q^a_1 * q^a_2 * q^a_3 * -q^a_4 * q^b_1 * -q^b_2 * q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_nw(t):
        return "(q^a_1 * q^a_2 * q^a_3 * q^a_4 * q^b_1 * -q^b_2 * q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_se(t):
        return "(q^a_1 * -q^a_2 * q^a_3 * -q^a_4 * q^b_1 * q^b_2 * q^b_3 * q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_un(t):
        return "(-q^a_1 * q^a_2 * q^a_3 * q^a_4 * -q^b_1 * q^b_2 * -q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_ds(t):
        return "(-q^a_1 * q^a_2 * -q^a_3 * -q^a_4 * -q^b_1 * q^b_2 * q^b_3 * q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_us(t):
        return "(-q^a_1 * q^a_2 * -q^a_3 * q^a_4 * -q^b_1 * q^b_2 * q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_dn(t):
        return "(-q^a_1 * q^a_2 * q^a_3 * -q^a_4 * -q^b_1 * q^b_2 * -q^b_3 * q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_ue(t):
        return "(q^a_1 * -q^a_2 * -q^a_3 * q^a_4 * q^b_1 * q^b_2 * -q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_dw(t):
        return "(q^a_1 * q^a_2 * -q^a_3 * -q^a_4 * q^b_1 * -q^b_2 * -q^b_3 * q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_uw(t):
        return "(q^a_1 * q^a_2 * -q^a_3 * q^a_4 * q^b_1 * -q^b_2 * -q^b_3 * -q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    def backwards_de(t):
        return "(q^a_1 * -q^a_2 * -q^a_3 * -q^a_4 * q^b_1 * q^b_2 * -q^b_3 * q^b_4)".replace('a', str(t)).replace('b', str(t + 1))

    back = ""
    for turn in range(num_amino - 2):
        back += backwards_ne(turn) + " + "
        back += backwards_sw(turn) + " + "
        back += backwards_nw(turn) + " + "
        back += backwards_se(turn) + " + "
        back += backwards_un(turn) + " + "
        back += backwards_ds(turn) + " + "
        back += backwards_us(turn) + " + "
        back += backwards_dn(turn) + " + "
        back += backwards_ue(turn) + " + "
        back += backwards_dw(turn) + " + "
        back += backwards_uw(turn) + " + "
        back += backwards_de(turn) + " + "

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
    dy_plus_sum_plus_one = sum_of_directions_plus_one('dy_plus', amino1, amino2)
    dy_minus_sum_plus_one = sum_of_directions_plus_one('dy_minus', amino1, amino2)
    dz_plus_sum_plus_one = sum_of_directions_plus_one('dz_plus', amino1, amino2)
    dz_minus_sum_plus_one = sum_of_directions_plus_one('dz_minus', amino1, amino2)

    x_equal = ''
    y_equal = ''
    z_equal = ''
    x_plus_one = ''
    x_minus_one = ''
    y_plus_one = ''
    y_minus_one = ''
    z_plus_one = ''
    z_minus_one = ''

    for bit in range(math.ceil(math.log2(amino2 - amino1))):
        x_equal += xnor(dx_plus_sum[bit], dx_minus_sum[bit]) + " + "
        y_equal += xnor(dy_plus_sum[bit], dy_minus_sum[bit]) + " + "
        z_equal += xnor(dz_plus_sum[bit], dz_minus_sum[bit]) + " + "
        x_plus_one += xnor(dx_plus_sum_plus_one[bit], dx_minus_sum[bit]) + " + "
        x_minus_one += xnor(dx_plus_sum[bit], dx_minus_sum_plus_one[bit]) + " + "
        y_plus_one += xnor(dy_plus_sum_plus_one[bit], dy_minus_sum[bit]) + " + "
        y_minus_one += xnor(dy_plus_sum[bit], dy_minus_sum_plus_one[bit]) + " + "
        z_plus_one += xnor(dz_plus_sum_plus_one[bit], dz_minus_sum[bit]) + " + "
        z_minus_one = xnor(dz_plus_sum[bit], dz_minus_sum_plus_one[bit]) + " + "

    # x equal and while y is +1 or -1 and while z is +1 or -1
    x_equal_yz_offset =  '((' + x_equal[:-3] + ') * ((' + y_plus_one[:-3] + ') + (' + y_minus_one[:-3] + ')) * ((' + z_plus_one[:-3] + ') + (' + z_minus_one[:-3] + ')))'

    # y equal while x is +1 or -1 and while z is +1 or -1
    y_equal_xz_offset = '((' + y_equal[:-3] + ') * ((' + x_plus_one[:-3] + ') + (' + x_minus_one[:-3] + ')) * ((' + z_plus_one[:-3] + ') + (' + z_minus_one[:-3] + ')))'

    # z equal while x is +1 or -1 and while y is +1 or -1
    z_equal_xy_offset = '((' + z_equal[:-3] + ') * ((' + x_plus_one[:-3] + ') + (' + x_minus_one[:-3] + ')) * ((' + y_plus_one[:-3] + ') + (' + y_minus_one[:-3] + ')))'

    return x_equal_yz_offset + ' + ' + y_equal_xz_offset + ' + ' + z_equal_xy_offset
