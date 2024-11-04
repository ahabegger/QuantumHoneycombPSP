"""
Direction & Movement Vector & Encoding $(q^t_1 q^t_2)$
East (E)  & $[1, 0]$   & $01$
West (W)  & $[-1, 0]$  & $10$
North (N) & $[0, 1]$   & $11$
South (S) & $[0,-1]$   & $00$

Presets : 01 E , x1 E or N
"""

import math
from Energy import get_energy_matrix
from sympy import *
from BitOps import sum_of_directions, sum_of_directions_plus_one, initialize_q_vars
from sympy.logic.boolalg import Or, And, Not, true, false, Xnor


def create_energy_function(sequence, energy_model):
    num_amino = len(sequence)
    global q_vars
    q_vars = initialize_q_vars(num_amino, 2)

    overlap = create_overlap_constraint(num_amino)
    overlap = set_default_and_dnf(overlap)
    print(f"overlap: {overlap}")

    backup = create_back_constraint(num_amino)
    backup = set_default_and_dnf(backup)
    print(f"backup: {backup}")

    energy_matrix = get_energy_matrix(sequence, energy_model)
    create_interactions(sequence, energy_matrix)


def set_default_and_dnf(expression):
    try:
        new_expression = expression.subs(q_vars[(0, 0)], False)
        new_expression = new_expression.subs(q_vars[(0, 1)], True)
        new_expression = new_expression.subs(q_vars[(1, 1)], True)
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
    return And(q_ta, q_tb)


def dy_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    return And(Not(q_ta), Not(q_tb))


def create_overlap_constraint(num_amino):
    # Initialize overlap constraint as False
    overlap = false

    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 3, num_amino):
            dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
            dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
            dy_plus_sum = sum_of_directions(dy_plus, amino1, amino2)
            dy_minus_sum = sum_of_directions(dy_minus, amino1, amino2)

            bits = math.ceil(math.log2(amino2 - amino1))
            amino_overlap = true
            for bit in range(bits):
                bit_overlap_flag = Xnor(dx_plus_sum[bit], dx_minus_sum[bit]) # Xnor is equivalent to '=='
                bit_overlap_flag = And(bit_overlap_flag, Xnor(dy_plus_sum[bit], dy_minus_sum[bit]))
                amino_overlap = And(amino_overlap, bit_overlap_flag)
            overlap = Or(overlap, amino_overlap)

    return overlap


def create_back_constraint(num_amino):
    # Initialize back constraint as False
    back = false
    for t in range(num_amino - 2):
        back = Or(back, And(dx_plus(t), dx_minus(t + 1)))
        back = Or(back, And(dx_minus(t), dx_plus(t + 1)))
        back = Or(back, And(dy_plus(t), dy_minus(t + 1)))
        back = Or(back, And(dy_minus(t), dy_plus(t + 1)))
    return back


def adjacency_indicator(amino1, amino2):
    dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
    dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
    dy_plus_sum = sum_of_directions(dy_plus, amino1, amino2)
    dy_minus_sum = sum_of_directions(dy_minus, amino1, amino2)
    dx_plus_sum_plus_one = sum_of_directions_plus_one(dx_plus, amino1, amino2)
    dx_minus_sum_plus_one = sum_of_directions_plus_one(dx_minus, amino1, amino2)
    dy_plus_sum_plus_one = sum_of_directions_plus_one(dy_plus, amino1, amino2)
    dy_minus_sum_plus_one = sum_of_directions_plus_one(dy_minus, amino1, amino2)

    bits = math.ceil(math.log2(amino2 - amino1))

    x_equal = true
    y_equal = true
    x_plus_one = true
    x_minus_one = true
    y_plus_one = true
    y_minus_one = true

    for bit in range(bits):
        x_equal = And(x_equal, Xnor(dx_plus_sum[bit], dx_minus_sum[bit]))
        y_equal = And(y_equal, Xnor(dy_plus_sum[bit], dy_minus_sum[bit]))

    bits_2 = math.ceil(math.log2(amino2 - amino1 + 1))

    for bit in range(bits_2):
        x_plus_one = And(x_plus_one, Xnor(dx_plus_sum_plus_one[bit], dx_minus_sum[bit]))
        x_minus_one = And(x_minus_one, Xnor(dx_plus_sum[bit], dx_minus_sum_plus_one[bit]))
        y_plus_one = And(y_plus_one, Xnor(dy_plus_sum_plus_one[bit], dy_minus_sum[bit]))
        y_minus_one = And(y_minus_one, Xnor(dy_plus_sum[bit], dy_minus_sum_plus_one[bit]))

    x_equal_y_offset = And(x_equal, Or(y_plus_one, y_minus_one))
    y_equal_x_offset = And(y_equal, Or(x_plus_one, x_minus_one))

    return Or(x_equal_y_offset, y_equal_x_offset)


def create_interactions(sequence, energy_matrix):
    num_amino = len(sequence)
    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 3, num_amino):
            energy_value = energy_matrix[amino1][amino2]
            if energy_value != 0:
                interaction = adjacency_indicator(amino1, amino2)
                interaction = set_default_and_dnf(interaction)
                print(f"interaction ({amino1} - {amino2}): {energy_value} * {interaction}")



if __name__ == '__main__':
    sequence = 'GACAA'
    energy_model = 'HP'  # Placeholder
    create_energy_function(sequence, energy_model)
