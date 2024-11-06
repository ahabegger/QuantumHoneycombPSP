import argparse

from Energy import encode_hp, encode_hpab, get_energy_matrix
from Binary import HCOMB4, HCOMB6, HCOMB8, HCOMB12
from QUBO import HCOMB4_QUBO, HCOMB6_QUBO, HCOMB8_QUBO, HCOMB12_QUBO


def main(sequence, energy_model, lattice_type, binary):
    if energy_model == 'HP':
        encoded_sequence = encode_hp(sequence)
    elif energy_model == 'HPAB':
        encoded_sequence = encode_hpab(sequence)
    elif energy_model == 'WHPAB':
        encoded_sequence = encode_hpab(sequence)
    else:
        encoded_sequence = sequence

    length = len(sequence)

    print(f'SEQUENCE:\t\t{sequence}')
    print(f'LATTICE_TYPE:\t\t{lattice_type}')
    print(f'ENERGY_MODEL:\t\t{energy_model}')
    print(f'ENCODED_SEQUENCE:\t\t{encoded_sequence}')
    print(f"LENGTH:\t\t{length} acids")

    interaction_matrix = get_energy_matrix(sequence, energy_model)

    print(f'INTERACTION_MATRIX:\n{interaction_matrix}')

    if binary:
        if lattice_type == 4:
            energy_function = HCOMB4.create_energy_function(sequence, energy_model)
        elif lattice_type == 6:
            energy_function = HCOMB6.create_energy_function(sequence, energy_model)
        elif lattice_type == 8:
            energy_function = HCOMB8.create_energy_function(sequence, energy_model)
        else:
            energy_function = HCOMB12.create_energy_function(sequence, energy_model)

        exit(0)

    if lattice_type == 4:
        qubits_per_amino = 2
        model, bqm, qubo, ising = HCOMB4_QUBO.create_energy_function(sequence, energy_model)
    elif lattice_type == 6:
        qubits_per_amino = 3
        model, bqm, qubo, ising = HCOMB6_QUBO.create_energy_function(sequence, energy_model)
    elif lattice_type == 8:
        qubits_per_amino = 3
        model, bqm, qubo, ising = HCOMB8_QUBO.create_energy_function(sequence, energy_model)
    else:
        qubits_per_amino = 4
        model, bqm, qubo, ising = HCOMB12_QUBO.create_energy_function(sequence, energy_model)

    print(f'QUBO:\n{qubo}')
    print(f'BQM:\n{bqm}')
    print(f'ISING:\n{ising}')

    qubo = qubo[0]

    variables = []
    for t in range(length-1):
        for q in range(qubits_per_amino):
            variables.append(f'q_{t}{chr(ord('a') + q)}')

    if lattice_type == 4:
        variables.remove('q_0a')
        variables.remove('q_0b')
        variables.remove('q_1b')
    elif lattice_type == 6:
        variables.remove('q_0a')
        variables.remove('q_0b')
        variables.remove('q_0c')
        variables.remove('q_1b')
        variables.remove('q_1c')
    elif lattice_type == 8:
        variables.remove('q_0a')
        variables.remove('q_0b')
    elif lattice_type == 12:
        variables.remove('q_0a')
        variables.remove('q_0b')
        variables.remove('q_0c')
        variables.remove('q_0d')
        variables.remove('q_1a')
        variables.remove('q_1d')


    qubo_vars = {}
    counter = 0
    for item in variables:
        qubo_vars[item] = f"q({counter})"
        counter += 1

    print(f'QUBO_VARS:\n{qubo_vars}')

    def replace_old(expression, dictionary):
        for key, value in dictionary.items():
            expression = expression.replace(key, value)
        return expression.replace(' ', '')

    new_qubo = {}
    for key, value in qubo.items():
        new_key = (replace_old(key[0], qubo_vars),
                   replace_old(key[1], qubo_vars))
        new_qubo[new_key] = value

    print(f'NEW_QUBO:\n{new_qubo}')

    new_eq = ''
    for key, value in new_qubo.items():
        if key[0] == key[1]:
            new_eq += f'{value}*{key[0]}+'
        else:
            new_eq += f'{value}*{key[0]}*{key[1]}+'

    new_eq = new_eq[:-1]
    new_eq = new_eq.replace('+-', '-')
    new_eq = new_eq.replace(' ', '')

    print(f'NEW_EQ:\n{new_eq}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='QuantumHoneycombPSP CLI')
    parser.add_argument('sequence', type=str, help='Protein sequence')
    parser.add_argument('energy_model', type=str, choices=['HP', 'HPAB', 'WHPAB', 'MJ'], help='Energy model')
    parser.add_argument('lattice_type', type=int, choices=[4, 6, 8, 12], help='Lattice type')
    parser.add_argument('--binary', action='store_true', help='Use Binary model')

    args = parser.parse_args()
    main(args.sequence, args.energy_model, args.lattice_type, args.binary)