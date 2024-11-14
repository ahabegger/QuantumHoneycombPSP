import argparse
import time
import pickle

from Annealer import annealer
from Energy import encode_hp, encode_hpab, get_energy_matrix
from Binary import HCOMB4, HCOMB6, HCOMB8, HCOMB12
from QUBO import HCOMB4_QUBO, HCOMB6_QUBO, HCOMB8_QUBO, HCOMB12_QUBO
from Sample_Analysis import sample_analysis


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
        print(f'ENERGY_FUNCTION:\n{energy_function}')
        exit(0)

    if lattice_type == 4:
        model, bqm, qubo, ising = HCOMB4_QUBO.create_energy_function(sequence, energy_model)
    elif lattice_type == 6:
        model, bqm, qubo, ising = HCOMB6_QUBO.create_energy_function(sequence, energy_model)
    elif lattice_type == 8:
        model, bqm, qubo, ising = HCOMB8_QUBO.create_energy_function(sequence, energy_model)
    else:
        model, bqm, qubo, ising = HCOMB12_QUBO.create_energy_function(sequence, energy_model)

    # Pickle the quadratic problem
    with open(f'Results/{sequence}_{lattice_type}_{energy_model}_BQM_{int(time.time())}.pkl', 'wb') as f:
        pickle.dump(bqm, f)

    # RUN D Wave Annealer
    samples = annealer(bqm)

    # Run Sample Analysis
    sample_analysis(samples)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='QuantumHoneycombPSP CLI')
    parser.add_argument('sequence', type=str, help='Protein sequence')
    parser.add_argument('energy_model', type=str, choices=['HP', 'HPAB', 'WHPAB', 'MJ'], help='Energy model')
    parser.add_argument('lattice_type', type=int, choices=[4, 6, 8, 12], help='Lattice type')
    parser.add_argument('--binary', action='store_true', help='Use Binary model')

    args = parser.parse_args()
    main(args.sequence, args.energy_model, args.lattice_type, args.binary)
