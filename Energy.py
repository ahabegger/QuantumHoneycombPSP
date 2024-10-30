import numpy as np

def get_energy_matrix(sequence, energy_model):
    if energy_model == 'HP':
        energy_table = get_energy_hp_table()
        encoded_sequence = encode_hp(sequence)
    elif energy_model == 'HPAB':
        energy_table = get_energy_hpab_table()
        encoded_sequence = encode_hpab(sequence)
    elif energy_model == 'WHPAB':
        energy_table = get_energy_whpab_table()
        encoded_sequence = encode_hpab(sequence)
    elif energy_model == 'MJ':
        energy_table = get_energy_mj_table()
        encoded_sequence = sequence
    else:
        raise ValueError('Invalid energy model.')

    # create a df to store the energy matrix
    energy_matrix = np.zeros((len(sequence), len(sequence)))
    energy_values = []
    for i in range(len(encoded_sequence)):
        for j in range(len(encoded_sequence)):
            if energy_model == 'MJ':
                energy_value = energy_table[encoded_sequence[i]][encoded_sequence[j]]
                energy_matrix[i][j] = energy_value
                energy_values.append(energy_value)
            else:
                energy_matrix[i][j] = energy_table[encoded_sequence[i] + encoded_sequence[j]]

    if energy_model == 'MJ':
        energies = list(set(energy_values))
        energies.sort(reverse = True)
        ranked_energy_values = []
        for i in range(len(energies)):
            ranked_energy_values.append(-1 - i)

        energy_dict = dict(zip(energies, ranked_energy_values))

        for i in range(len(encoded_sequence)):
            for j in range(len(encoded_sequence)):
                energy_matrix[i][j] = energy_dict[energy_matrix[i][j]]


    return energy_matrix


def encode_hp(sequence):
    """
    Encodes a given amino acid sequence into the HP model.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to be encoded.

    Returns
    -------
    str
        Encoded sequence in the HP model.
    """
    hp_sequence = ''
    for amino_acid in sequence:
        if amino_acid in 'AGILMFPWV':
            hp_sequence += 'H'
        else:
            hp_sequence += 'P'
    return hp_sequence


def encode_hpab(sequence):
    """
    Encodes a given amino acid sequence into the HPAB model.

    Parameters
    ----------
    sequence : str
        Amino acid sequence to be encoded.

    Returns
    -------
    str
        Encoded sequence in the HPAB model.
    """
    hpab_sequence = ''
    for amino_acid in sequence:
        if amino_acid in 'AGILMFPWV':
            hpab_sequence += 'H'
        elif amino_acid in 'NQSTY':
            hpab_sequence += 'P'
        elif amino_acid in 'DE':
            hpab_sequence += 'A'
        else:
            hpab_sequence += 'B'
    return hpab_sequence

def get_energy_hp_table():
    """
    Get the energy table for the HP model.

    Returns
    -------
    dict
        Energy table for the HP model.
    """
    return {'HH': 1, 'HP': 0, 'PH': 0, 'PP': 0}

def get_energy_hpab_table():
    """
    Get the energy table for the HPAB model.

    Returns
    -------
    dict
        Energy table for the HPAB model.
    """

    return {'HH': -1, 'HP': 0, 'HA': 0, 'HB': 0, 'PH': 0, 'PP': 0, 'PA': 0, 'PB': 0, 'AH': 0, 'AP': 0, 'AA': 1, 'AB': -1, 'BH': 0, 'BP': 0, 'BA': -1, 'BB': 1}

def get_energy_whpab_table():
    """
    Get the energy table for the weighted HPAB model.

    Returns
    -------
    dict
        Energy table for the weighted HPAB model.
    """
    return {'HH': -4, 'HP': 0, 'HA': -1, 'HB': -1, 'PH': 0, 'PP': 0, 'PA': 0, 'PB': 0, 'AH': -1, 'AP': 0, 'AA': 2, 'AB': -2, 'BH': -1, 'BP': 0, 'BA': -2, 'BB': 2}

def get_energy_mj_table():
    """
    Get the energy table for the Miyazawa-Jernigan model.
    :return:
    """
    mj_contact_energies = {
        "C": {"C": -5.44, "M": -4.99, "F": -5.80, "I": -5.50, "L": -5.83, "V": -4.96, "W": -4.95,
              "Y": -4.16, "A": -3.57, "G": -3.16,
              "T": -3.11, "S": -2.86, "N": -2.59, "Q": -2.85, "D": -2.41, "E": -2.27, "H": -3.60,
              "R": -2.57, "K": -1.95, "P": -3.07},
        "M": {"C": -4.99, "M": -5.46, "F": -6.56, "I": -6.02, "L": -6.41, "V": -5.32, "W": -5.55,
              "Y": -4.91, "A": -3.94, "G": -3.39,
              "T": -3.51, "S": -3.03, "N": -2.95, "Q": -3.30, "D": -2.57, "E": -2.89, "H": -3.98,
              "R": -3.12, "K": -2.48, "P": -3.45},
        "F": {"C": -5.80, "M": -6.56, "F": -7.26, "I": -6.84, "L": -7.28, "V": -6.29, "W": -6.16,
              "Y": -5.66, "A": -4.81, "G": -4.13,
              "T": -4.28, "S": -4.02, "N": -3.75, "Q": -4.10, "D": -3.48, "E": -3.56, "H": -4.77,
              "R": -3.98, "K": -3.36, "P": -4.25},
        "I": {"C": -5.50, "M": -6.02, "F": -6.84, "I": -6.54, "L": -7.04, "V": -6.05, "W": -5.78,
              "Y": -5.25, "A": -4.58, "G": -3.78,
              "T": -4.03, "S": -3.52, "N": -3.24, "Q": -3.67, "D": -3.17, "E": -3.27, "H": -4.14,
              "R": -3.63, "K": -3.01, "P": -3.76},
        "L": {"C": -5.83, "M": -6.41, "F": -7.28, "I": -7.04, "L": -7.37, "V": -6.48, "W": -6.14,
              "Y": -5.67, "A": -4.91, "G": -4.16,
              "T": -4.34, "S": -3.92, "N": -3.74, "Q": -4.04, "D": -3.40, "E": -3.59, "H": -4.54,
              "R": -4.03, "K": -3.37, "P": -4.20},
        "V": {"C": -4.96, "M": -5.32, "F": -6.29, "I": -6.05, "L": -6.48, "V": -5.52, "W": -5.18,
              "Y": -4.62, "A": -4.04, "G": -3.38,
              "T": -3.46, "S": -3.05, "N": -2.83, "Q": -3.07, "D": -2.48, "E": -2.67, "H": -3.58,
              "R": -3.07, "K": -2.49, "P": -3.32},
        "W": {"C": -4.95, "M": -5.55, "F": -6.16, "I": -5.78, "L": -6.14, "V": -5.18, "W": -5.06,
              "Y": -4.66, "A": -3.82, "G": -3.42,
              "T": -3.22, "S": -2.99, "N": -3.07, "Q": -3.11, "D": -2.84, "E": -2.99, "H": -3.98,
              "R": -3.41, "K": -2.69, "P": -3.73},
        "Y": {"C": -4.16, "M": -4.91, "F": -5.66, "I": -5.25, "L": -5.67, "V": -4.62, "W": -4.66,
              "Y": -4.17, "A": -3.36, "G": -3.01,
              "T": -3.01, "S": -2.78, "N": -2.76, "Q": -2.97, "D": -2.76, "E": -2.79, "H": -3.52,
              "R": -3.16, "K": -2.60, "P": -3.19},
        "A": {"C": -3.57, "M": -3.94, "F": -4.81, "I": -4.58, "L": -4.91, "V": -4.04, "W": -3.82,
              "Y": -3.36, "A": -2.72, "G": -2.31,
              "T": -2.32, "S": -2.01, "N": -1.84, "Q": -1.89, "D": -1.70, "E": -1.51, "H": -2.41,
              "R": -1.83, "K": -1.31, "P": -2.03},
        "G": {"C": -3.16, "M": -3.39, "F": -4.13, "I": -3.78, "L": -4.16, "V": -3.38, "W": -3.42,
              "Y": -3.01, "A": -2.31, "G": -2.24,
              "T": -2.08, "S": -1.82, "N": -1.74, "Q": -1.66, "D": -1.59, "E": -1.22, "H": -2.15,
              "R": -1.72, "K": -1.15, "P": -1.87},
        "T": {"C": -3.11, "M": -3.51, "F": -4.28, "I": -4.03, "L": -4.34, "V": -3.46, "W": -3.22,
              "Y": -3.01, "A": -2.32, "G": -2.08,
              "T": -2.12, "S": -1.96, "N": -1.88, "Q": -1.90, "D": -1.80, "E": -1.74, "H": -2.42,
              "R": -1.90, "K": -1.31, "P": -1.90},
        "S": {"C": -2.86, "M": -3.03, "F": -4.02, "I": -3.52, "L": -3.92, "V": -3.05, "W": -2.99,
              "Y": -2.78, "A": -2.01, "G": -1.82,
              "T": -1.67, "S": -1.67, "N": -1.58, "Q": -1.49, "D": -1.63, "E": -1.48, "H": -2.11,
              "R": -1.62, "K": -1.05, "P": -1.57},
        "N": {"C": -2.59, "M": -2.95, "F": -3.75, "I": -3.24, "L": -3.74, "V": -2.83, "W": -3.07,
              "Y": -2.76, "A": -1.84, "G": -1.74,
              "T": -1.88, "S": -1.58, "N": -1.68, "Q": -1.71, "D": -1.68, "E": -1.51, "H": -2.08,
              "R": -1.64, "K": -1.21, "P": -1.53},
        "Q": {"C": -2.85, "M": -3.30, "F": -4.10, "I": -3.67, "L": -4.04, "V": -3.07, "W": -3.11,
              "Y": -2.97, "A": -1.89, "G": -1.66,
              "T": -1.90, "S": -1.49, "N": -1.71, "Q": -1.54, "D": -1.46, "E": -1.42, "H": -1.98,
              "R": -1.80, "K": -1.29, "P": -1.73},
        "D": {"C": -2.41, "M": -2.57, "F": -3.48, "I": -3.17, "L": -3.40, "V": -2.48, "W": -2.84,
              "Y": -2.76, "A": -1.70, "G": -1.59,
              "T": -1.80, "S": -1.63, "N": -1.68, "Q": -1.21, "D": -1.21, "E": -1.02, "H": -2.32,
              "R": -2.29, "K": -1.68, "P": -1.33},
        "E": {"C": -2.27, "M": -2.89, "F": -3.56, "I": -3.27, "L": -3.59, "V": -2.67, "W": -2.99,
              "Y": -2.79, "A": -1.51, "G": -1.22,
              "T": -1.74, "S": -1.48, "N": -1.51, "Q": -1.42, "D": -1.02, "E": -0.91, "H": -2.15,
              "R": -2.27, "K": -1.80, "P": -1.26},
        "H": {"C": -3.60, "M": -3.98, "F": -4.77, "I": -4.14, "L": -4.54, "V": -3.58, "W": -3.98,
              "Y": -3.52, "A": -2.41, "G": -2.15,
              "T": -2.42, "S": -2.11, "N": -2.08, "Q": -1.98, "D": -2.32, "E": -2.15, "H": -3.05,
              "R": -2.16, "K": -1.35, "P": -2.25},
        "R": {"C": -2.57, "M": -3.12, "F": -3.98, "I": -3.63, "L": -4.03, "V": -3.07, "W": -3.41,
              "Y": -3.16, "A": -1.83, "G": -1.72,
              "T": -1.90, "S": -1.62, "N": -1.64, "Q": -1.80, "D": -2.29, "E": -2.27, "H": -2.16,
              "R": -1.55, "K": -0.59, "P": -1.70},
        "K": {"C": -1.95, "M": -2.48, "F": -3.36, "I": -3.01, "L": -3.37, "V": -2.49, "W": -2.69,
              "Y": -2.60, "A": -1.31, "G": -1.15,
              "T": -1.31, "S": -1.05, "N": -1.21, "Q": -1.29, "D": -1.68, "E": -1.80, "H": -1.35,
              "R": -0.59, "K": -0.12, "P": -0.97},
        "P": {"C": -3.07, "M": -3.45, "F": -4.25, "I": -3.76, "L": -4.20, "V": -3.32, "W": -3.73,
              "Y": -3.19, "A": -2.03, "G": -1.87,
              "T": -1.90, "S": -1.57, "N": -1.53, "Q": -1.73, "D": -1.33, "E": -1.26, "H": -2.25,
              "R": -1.70, "K": -0.97, "P": -1.75}
    }

    return mj_contact_energies