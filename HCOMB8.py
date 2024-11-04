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

Presets : 00x S or U
"""

from BitOps import sum_of_y, sum_of_y_plus_two, sum_of_y_plus_one
import math
from Energy import get_energy_matrix
from sympy import *
from BitOps import sum_of_directions, sum_of_directions_plus_one, initialize_q_vars
from sympy.logic.boolalg import Or, And, Not, true, false, Xnor

def create_energy_function(sequence, energy_model):
    num_amino = len(sequence)
    global q_vars
    q_vars = initialize_q_vars(num_amino, 3)

    overlap = create_overlap_constraint(num_amino)
    overlap = set_default_and_dnf(overlap)
    print(f"overlap: {overlap}")

    energy_matrix = get_energy_matrix(sequence, energy_model)
    create_interactions(sequence, energy_matrix)


def set_default_and_dnf(expression):
    try:
        new_expression = expression.subs(q_vars[(0, 0)], False)
        new_expression = new_expression.subs(q_vars[(0, 1)], False)
        new_expression = to_dnf(new_expression, True, True)
        return new_expression
    except:
        pass


def dx_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    return And(Not(q_ta), q_tb)

def dx_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    return And(q_ta, Not(q_tb))

def dy_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return Or(And(q_ta, q_tb, q_tc), And(Not(q_ta), q_tb, q_tc), And(q_ta, Not(q_tb), q_tc))

def dy_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return Or(And(Not(q_ta), Not(q_tb), Not(q_tc)), And(q_ta, Not(q_tb), Not(q_tc)), And(Not(q_ta), q_tb, Not(q_tc)))

def dy_plus_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(q_ta, q_tb, q_tc)

def dy_minus_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(Not(q_ta), Not(q_tb), Not(q_tc))

def dz_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(Not(q_ta), Not(q_tb), q_tc)

def dz_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(q_ta, q_tb, Not(q_tc))


def create_overlap_constraint(num_amino):
    # Initialize overlap constraint as False
    overlap = false

    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 2, num_amino):
            dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
            dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
            dy_plus_sum = sum_of_y(dy_plus, dy_plus_plus, amino1, amino2)
            dy_minus_sum = sum_of_y(dy_minus, dy_minus_minus, amino1, amino2)
            dz_plus_sum = sum_of_directions(dz_plus, amino1, amino2)
            dz_minus_sum = sum_of_directions(dz_minus, amino1, amino2)

            bits = math.ceil(math.log2(amino2 - amino1))
            amino_overlap = true
            for bit in range(bits):
                bit_overlap_flag = Xnor(dx_plus_sum[bit], dx_minus_sum[bit]) # Xnor is equivalent to '=='
                bit_overlap_flag = And(bit_overlap_flag, Xnor(dz_plus_sum[bit], dz_minus_sum[bit]))
                amino_overlap = And(amino_overlap, bit_overlap_flag)

            bits_2 = math.ceil(math.log2(amino2 - amino1)) + 1
            for bit in range(bits_2):
                bit_overlap_flag = Xnor(dy_plus_sum[bit], dy_minus_sum[bit])
                amino_overlap = And(amino_overlap, bit_overlap_flag)

            overlap = Or(overlap, amino_overlap)

    return overlap

def adjacency_indicator(amino1, amino2):
    dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
    dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
    dy_plus_sum = sum_of_y(dy_plus, dy_plus_plus, amino1, amino2)
    dy_minus_sum = sum_of_y(dy_minus, dy_minus_minus, amino1, amino2)
    dz_plus_sum = sum_of_directions(dz_plus, amino1, amino2)
    dz_minus_sum = sum_of_directions(dz_minus, amino1, amino2)

    dx_plus_sum_plus_one = sum_of_directions_plus_one(dx_plus, amino1, amino2)
    dx_minus_sum_plus_one = sum_of_directions_plus_one(dx_minus, amino1, amino2)
    dy_plus_sum_plus_one = sum_of_y_plus_one(dy_plus, dy_plus_plus, amino1, amino2)
    dy_minus_sum_plus_one = sum_of_y_plus_one(dy_minus, dy_minus_minus, amino1, amino2)
    dz_plus_sum_plus_one = sum_of_directions_plus_one(dz_plus, amino1, amino2)
    dz_minus_sum_plus_one = sum_of_directions_plus_one(dz_minus, amino1, amino2)

    dy_plus_sum_plus_two = sum_of_y_plus_two(dy_plus, dy_plus_plus, amino1, amino2)
    dy_minus_sum_plus_two = sum_of_y_plus_two(dy_minus, dy_minus_minus, amino1, amino2)

    x_equal = true
    y_equal = true
    z_equal = true

    x_plus_one = true
    x_minus_one = true
    y_plus_one = true
    y_minus_one = true
    z_plus_one = true
    z_minus_one = true

    y_plus_two = true
    y_minus_two = true

    bits = math.ceil(math.log2(amino2 - amino1))

    for bit in range(bits):
        x_equal = And(x_equal, Xnor(dx_plus_sum[bit], dx_minus_sum[bit]))
        z_equal = And(z_equal, Xnor(dz_plus_sum[bit], dz_minus_sum[bit]))

    bits = math.ceil(math.log2(amino2 - amino1)) + 1
    for bit in range(bits):
        y_equal = And(y_equal, Xnor(dy_plus_sum[bit], dy_minus_sum[bit]))

    bits = math.ceil(math.log2(amino2 - amino1 + 1))
    for bit in range(bits):
        x_plus_one = And(x_plus_one, Xnor(dx_plus_sum_plus_one[bit], dx_minus_sum[bit]))
        x_minus_one = And(x_minus_one, Xnor(dx_plus_sum[bit], dx_minus_sum_plus_one[bit]))
        z_plus_one = And(z_plus_one, Xnor(dz_plus_sum_plus_one[bit], dz_minus_sum[bit]))
        z_minus_one = And(z_minus_one, Xnor(dz_plus_sum[bit], dz_minus_sum_plus_one[bit]))

    bits = math.ceil(math.log2(amino2 - amino1 + 1)) + 1
    for bit in range(bits):
        y_plus_one = And(y_plus_one, Xnor(dy_plus_sum_plus_one[bit], dy_minus_sum[bit]))
        y_minus_one = And(y_minus_one, Xnor(dy_plus_sum[bit], dy_minus_sum_plus_one[bit]))

    bits = math.ceil(math.log2(amino2 - amino1 + 2)) + 1
    for bit in range(bits):
        y_plus_two = And(y_plus_two, Xnor(dy_plus_sum_plus_two[bit], dy_minus_sum[bit]))
        y_minus_two = And(y_minus_two, Xnor(dy_plus_sum[bit], dy_minus_sum_plus_two[bit]))

    # xz equal and y offset by 2
    xz_equal_y_offset = And(x_equal, z_equal, Or(y_plus_two, y_minus_two))
    # xy equal and z offset by 1
    xy_equal_z_offset = And(x_equal, y_equal, Or(z_plus_one, z_minus_one))
    # z equal and x offset by 1 and y offset by 1
    z_equal_xy_offset = And(z_equal, Or(x_plus_one, x_minus_one), Or(y_plus_one, y_minus_one))

    return Or(xz_equal_y_offset, xy_equal_z_offset, z_equal_xy_offset)


def create_interactions(sequence, energy_matrix):
    num_amino = len(sequence)
    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 2, num_amino):
            energy_value = energy_matrix[amino1][amino2]
            if energy_value != 0:
                interaction = adjacency_indicator(amino1, amino2)
                interaction = set_default_and_dnf(interaction)
                print(f"interaction ({amino1} - {amino2}): {energy_value} * {interaction}")


if __name__ == '__main__':
    sequence = 'GAGA'
    energy_model = 'HP'  # Placeholder
    create_energy_function(sequence, energy_model)
