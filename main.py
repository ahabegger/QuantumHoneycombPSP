import HCOMB12
import HCOMB6
import HCOMB8
from Energy import encode_hp, encode_hpab, get_energy_matrix
import HCOMB4


sequence = 'GAGA'
energy_model = 'HP'
lattice_type = 4

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

if lattice_type == 4:
    energy_function = HCOMB4.create_energy_function(sequence, energy_model)
elif lattice_type == 6:
    energy_function = HCOMB6.create_energy_function(sequence, energy_model)
elif lattice_type == 8:
    energy_function = HCOMB8.create_energy_function(sequence, energy_model)
else:
    energy_function = HCOMB12.create_energy_function(sequence, energy_model)

print(f'ENERGY_FUNCTION:\n{energy_function}')
