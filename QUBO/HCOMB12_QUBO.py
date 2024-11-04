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

Presets : 1011 NE , 1xx1 NE or NW or UE or UW
"""

import math
from Energy import get_energy_matrix
from pyqubo import *
from QUBO.BitOps_QUBO import sum_of_directions, sum_of_directions_plus_one, initialize_q_vars
from pprint import pprint
from QUBO.BitOps_QUBO import Xnor

def create_energy_function(sequence, energy_model):
    num_amino = len(sequence)
    global q_vars
    q_vars = initialize_q_vars(num_amino, 4)
    q_vars = set_default(q_vars)

    energy_matrix = get_energy_matrix(sequence, energy_model)
    interactions, energy_values = create_interactions(sequence, energy_matrix)
    total_interaction_energy = Num(0)
    for i in range(len(interactions)):
        interaction_value = interactions[i] * energy_values[i]
        total_interaction_energy = total_interaction_energy + interaction_value\

    penalty = (sum(energy_values) * -1) + 1

    redundancy = create_redundancy_constraint(num_amino)
    redundancy = Num(penalty) * redundancy

    overlap = create_overlap_constraint(num_amino)
    overlap = Num(penalty) * overlap

    backup = create_back_constraint(num_amino)
    backup = Num(penalty) * backup

    model = total_interaction_energy + overlap + backup + redundancy
    model = model.compile(10)

    bqm = model.to_bqm()
    qubo = model.to_qubo()
    ising = model.to_ising()

    return model, bqm, qubo, ising


def set_default(vars):
    vars[(0, 0)] = Num(1)
    vars[(0, 1)] = Num(0)
    vars[(0, 2)] = Num(1)
    vars[(0, 3)] = Num(1)
    vars[(1, 0)] = Num(1)
    vars[(1, 3)] = Num(1)
    return vars

def dx_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    return And(q_ta, Not(q_tb))

def dx_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    return And(q_ta, q_tb)

def dy_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    q_td = q_vars[(t, 3)]
    return Or(And(q_ta, And(q_tc, q_td)), And(Not(q_ta), And(q_tb, q_tc)))

def dy_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    q_td = q_vars[(t, 3)]
    return Or(And(q_ta, And(q_tc, Not(q_td))), And(Not(q_ta), And(q_tb, Not(q_tc))))

def dz_plus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    q_td = q_vars[(t, 3)]
    return Or(And(q_ta, And(Not(q_tc), q_td)), And(Not(q_ta), And(q_tb, q_td)))

def dz_minus(t):
    q_ta = q_vars[(t, 0)]
    q_tb = q_vars[(t, 1)]
    q_tc = q_vars[(t, 2)]
    q_td = q_vars[(t, 3)]
    return Or(And(q_ta, And(Not(q_tc), Not(q_td))), And(Not(q_ta), And(q_tb, Not(q_td))))


def create_redundancy_constraint(num_amino):
    redundancy = Num(0)

    #0000 0011 0001 0010 invalid directions
    for t in range(num_amino - 1):
        q_ta = q_vars[(t, 0)]
        q_tb = q_vars[(t, 1)]
        redundancy = Or(redundancy, And(Not(q_ta), Not(q_tb)))

    return redundancy


def create_back_constraint(num_amino):
    # Initialize back constraint as False
    back = Num(0)
    for t in range(num_amino - 2):
        back = Or(back, And(dx_plus(t), And(dy_plus(t), And(dx_minus(t + 1), dy_minus(t + 1)))))
        back = Or(back, And(dx_minus(t), And(dy_minus(t), And(dx_plus(t + 1), dy_plus(t + 1)))))
        back = Or(back, And(dy_plus(t), And(dz_plus(t), And(dy_minus(t + 1), dz_minus(t + 1)))))
        back = Or(back, And(dy_minus(t), And(dz_minus(t), And(dy_plus(t + 1), dz_plus(t + 1)))))
        back = Or(back, And(dz_plus(t), And(dx_plus(t), And(dz_minus(t + 1), dx_minus(t + 1)))))
        back = Or(back, And(dz_minus(t), And(dx_minus(t), And(dz_plus(t + 1), dx_plus(t + 1)))))
        back = Or(back, And(dx_plus(t), And(dy_minus(t), And(dx_minus(t + 1), dy_plus(t + 1)))))
        back = Or(back, And(dx_minus(t), And(dy_plus(t), And(dx_plus(t + 1), dy_minus(t + 1)))))
        back = Or(back, And(dy_plus(t), And(dz_minus(t), And(dy_minus(t + 1), dz_plus(t + 1)))))
        back = Or(back, And(dy_minus(t), And(dz_plus(t), And(dy_plus(t + 1), dz_minus(t + 1)))))
        back = Or(back, And(dz_plus(t), And(dx_minus(t), And(dz_minus(t + 1), dx_plus(t + 1)))))
        back = Or(back, And(dz_minus(t), And(dx_plus(t), And(dz_plus(t + 1), dx_minus(t + 1)))))

    return back


def create_overlap_constraint(num_amino):
    # Initialize overlap constraint as False
    overlap = Num(0)

    for amino1 in range(num_amino):
        for amino2 in range(amino1 + 3, num_amino):
            dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
            dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
            dy_plus_sum = sum_of_directions(dy_plus, amino1, amino2)
            dy_minus_sum = sum_of_directions(dy_minus, amino1, amino2)
            dz_plus_sum = sum_of_directions(dz_plus, amino1, amino2)
            dz_minus_sum = sum_of_directions(dz_minus, amino1, amino2)

            bits = math.ceil(math.log2(amino2 - amino1))
            amino_overlap = Num(1)
            for bit in range(bits):
                bit_overlap_flag = Xnor(dx_plus_sum[bit], dx_minus_sum[bit]) # Xnor is equivalent to '=='
                bit_overlap_flag = And(bit_overlap_flag, Xnor(dy_plus_sum[bit], dy_minus_sum[bit]))
                bit_overlap_flag = And(bit_overlap_flag, Xnor(dz_plus_sum[bit], dz_minus_sum[bit]))
                amino_overlap = And(amino_overlap, bit_overlap_flag)

            overlap = Or(overlap, amino_overlap)

    return overlap


def adjacency_indicator(amino1, amino2):
    dx_plus_sum = sum_of_directions(dx_plus, amino1, amino2)
    dx_minus_sum = sum_of_directions(dx_minus, amino1, amino2)
    dy_plus_sum = sum_of_directions(dy_plus, amino1, amino2)
    dy_minus_sum = sum_of_directions(dy_minus, amino1, amino2)
    dz_plus_sum = sum_of_directions(dz_plus, amino1, amino2)
    dz_minus_sum = sum_of_directions(dz_minus, amino1, amino2)

    dx_plus_sum_plus_one = sum_of_directions_plus_one(dx_plus, amino1, amino2)
    dx_minus_sum_plus_one = sum_of_directions_plus_one(dx_minus, amino1, amino2)
    dy_plus_sum_plus_one = sum_of_directions_plus_one(dy_plus, amino1, amino2)
    dy_minus_sum_plus_one = sum_of_directions_plus_one(dy_minus, amino1, amino2)
    dz_plus_sum_plus_one = sum_of_directions_plus_one(dz_plus, amino1, amino2)
    dz_minus_sum_plus_one = sum_of_directions_plus_one(dz_minus, amino1, amino2)

    x_equal = Num(1)
    y_equal = Num(1)
    z_equal = Num(1)

    x_plus_one = Num(1)
    x_minus_one = Num(1)
    y_plus_one = Num(1)
    y_minus_one = Num(1)
    z_plus_one = Num(1)
    z_minus_one = Num(1)

    bits = math.ceil(math.log2(amino2 - amino1))

    for bit in range(bits):
        x_equal = And(x_equal, Xnor(dx_plus_sum[bit], dx_minus_sum[bit]))
        y_equal = And(y_equal, Xnor(dy_plus_sum[bit], dy_minus_sum[bit]))
        z_equal = And(z_equal, Xnor(dz_plus_sum[bit], dz_minus_sum[bit]))

    bits_2 = math.ceil(math.log2(amino2 - amino1 + 1))

    for bit in range(bits_2):
        x_plus_one = And(x_plus_one, Xnor(dx_plus_sum_plus_one[bit], dx_minus_sum[bit]))
        x_minus_one = And(x_minus_one, Xnor(dx_plus_sum[bit], dx_minus_sum_plus_one[bit]))
        y_plus_one = And(y_plus_one, Xnor(dy_plus_sum_plus_one[bit], dy_minus_sum[bit]))
        y_minus_one = And(y_minus_one, Xnor(dy_plus_sum[bit], dy_minus_sum_plus_one[bit]))
        z_plus_one = And(z_plus_one, Xnor(dz_plus_sum_plus_one[bit], dz_minus_sum[bit]))
        z_minus_one = And(z_minus_one, Xnor(dz_plus_sum[bit], dz_minus_sum_plus_one[bit]))

    # x equal yz offset
    x_equal_yz_offset = And(x_equal, And(Or(y_plus_one, y_minus_one), Or(z_plus_one, z_minus_one)))

    # y equal xz offset
    y_equal_xz_offset = And(y_equal, And(Or(x_plus_one, x_minus_one), Or(z_plus_one, z_minus_one)))

    # z equal xy offset
    z_equal_xy_offset = And(z_equal, And(Or(x_plus_one, x_minus_one), Or(y_plus_one, y_minus_one)))

    return Or(x_equal_yz_offset, Or(y_equal_xz_offset, z_equal_xy_offset))


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