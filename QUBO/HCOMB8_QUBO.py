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

import math
from Energy import get_energy_matrix
from pyqubo import *
from QUBO.BitOps_QUBO import (sum_of_directions, sum_of_directions_plus_one,
                         initialize_q_vars, sum_of_y, sum_of_y_plus_one, sum_of_y_plus_two)
from pprint import pprint
from QUBO.BitOps_QUBO import Xnor

def create_energy_function(sequence, energy_model):
    num_amino = len(sequence)
    global q_vars
    q_vars = initialize_q_vars(num_amino, 3)
    q_vars = set_default(q_vars)

    energy_matrix = get_energy_matrix(sequence, energy_model)
    interactions, energy_values = create_interactions(sequence, energy_matrix)
    total_interaction_energy = Num(0)
    for i in range(len(interactions)):
        interaction_value = interactions[i] * energy_values[i]
        total_interaction_energy = total_interaction_energy + interaction_value

    penalty = 40

    back = create_back_constraint(num_amino)
    back = Num(penalty) * back

    overlap = create_overlap_constraint(num_amino)
    overlap = Num(penalty) * overlap

    model = total_interaction_energy + overlap + back
    model = model.compile(1)

    bqm = model.to_bqm()
    qubo = model.to_qubo()
    ising = model.to_ising()

    return model, bqm, qubo, ising

def set_default(vars):
    vars[(0, 0)] = Num(0)
    vars[(0, 1)] = Num(0)
    return vars

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
    return Or(And(q_ta, And(q_tb, q_tc)),
              Or(And(Not(q_ta), And(q_tb, q_tc)),
                 And(q_ta, And(Not(q_tb), q_tc))))

def dy_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return Or(And(Not(q_ta), And(Not(q_tb), Not(q_tc))),
              Or(And(q_ta, And(Not(q_tb), Not(q_tc))),
                 And(Not(q_ta), And(q_tb, Not(q_tc)))))

def dy_plus_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(q_ta, And(q_tb, q_tc))

def dy_minus_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(Not(q_ta), And(Not(q_tb), Not(q_tc)))

def dz_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(Not(q_ta), And(Not(q_tb), q_tc))

def dz_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    return And(q_ta, And(q_tb, Not(q_tc)))


def create_back_constraint(num_amino):
    # Initialize back constraint as False
    back = Num(0)

    for t in range(num_amino - 2):
        q_ta = q_vars[(t, 0)]
        q_tb = q_vars[(t, 1)]
        q_tc = q_vars[(t, 2)]
        q_t2a = q_vars[(t + 1, 0)]
        q_t2b = q_vars[(t + 1, 1)]
        q_t2c = q_vars[(t + 1, 2)]

        back = Or(back, And(Not(q_ta), And(q_tb, And(q_tc, And(q_t2a, And(Not(q_t2b), Not(q_t2c))))))) # ne
        back = Or(back, And(q_ta, And(Not(q_tb), And(Not(q_tc), And(Not(q_t2a), And(q_t2b, q_t2c)))))) # sw

        back = Or(back, And(q_ta, And(Not(q_tb), And(q_tc, And(Not(q_t2a), And(q_t2b, Not(q_t2c))))))) # nw
        back = Or(back, And(Not(q_ta), And(q_tb, And(Not(q_tc), And(q_t2a, And(Not(q_t2b), q_t2c)))))) # se

        back = Or(back, And(dy_plus_plus(t), dy_minus_minus(t + 1))) # n
        back = Or(back, And(dy_minus_minus(t), dy_plus_plus(t + 1))) # s
        back = Or(back, And(dz_plus(t), dz_minus(t + 1))) # u
        back = Or(back, And(dz_minus(t), dz_plus(t + 1))) # d

    return back

def create_overlap_constraint(num_amino):
    # Initialize overlap constraint as False
    overlap = Num(0)

    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 3, num_amino):
            dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
            dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
            dy_plus_sum = sum_of_y(dy_plus, dy_plus_plus, amino1, amino2)
            dy_minus_sum = sum_of_y(dy_minus, dy_minus_minus, amino1, amino2)
            dz_plus_sum = sum_of_directions(dz_plus, amino1, amino2)
            dz_minus_sum = sum_of_directions(dz_minus, amino1, amino2)

            bits = math.ceil(math.log2(amino2 - amino1))
            amino_overlap = Num(1)
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

    x_equal = Num(1)
    y_equal = Num(1)
    z_equal = Num(1)

    x_plus_one = Num(1)
    x_minus_one = Num(1)
    y_plus_one = Num(1)
    y_minus_one = Num(1)
    z_plus_one = Num(1)
    z_minus_one = Num(1)

    y_plus_two = Num(1)
    y_minus_two = Num(1)

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
    xz_equal_y_offset = And(And(x_equal, z_equal), Or(y_plus_two, y_minus_two))
    # xy equal and z offset by 1
    xy_equal_z_offset = And(And(x_equal, y_equal), Or(z_plus_one, z_minus_one))
    # z equal and x offset by 1 and y offset by 1
    z_equal_xy_offset = And(And(z_equal, Or(x_plus_one, x_minus_one)), Or(y_plus_one, y_minus_one))

    return Or(xz_equal_y_offset, Or(xy_equal_z_offset, z_equal_xy_offset))


def create_interactions(sequence, energy_matrix):
    interations = []
    energy_values = []
    num_amino = len(sequence)

    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 2, num_amino):
            energy_value = energy_matrix[amino1][amino2]
            if energy_value != 0:
                interaction = adjacency_indicator(amino1, amino2)
                interations.append(interaction)
                energy_values.append(energy_value)

    return interations, energy_values


if __name__ == '__main__':
    sequence = 'GAAA'
    energy_model = 'HP'  # Placeholder
    bqm, qubo, ising = create_energy_function(sequence, energy_model)

    pprint(bqm)
    pprint(qubo)
    pprint(ising)
