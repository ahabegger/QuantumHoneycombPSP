import argparse
import math
from Energy import encode_hp, encode_hpab
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def visualize(sequence, energy_model, lattice_type, binary_output):
    # Uses Matplotlib to visualize the protein sequence on a lattice

    if lattice_type == 4:
        visualize_4(sequence, energy_model, binary_output)
    elif lattice_type == 6:
        visualize_6(sequence, energy_model, binary_output)
    elif lattice_type == 8:
        visualize_8(sequence, energy_model, binary_output)
    else:
        visualize_12(sequence, energy_model, binary_output)


def get_interaction_coords(coordinates):
    # cycle through the coordinates and return pairs of coordinates that are 1 unit apart
    interaction_coords = []
    connection_coords = []

    for i in range(len(coordinates)):
        for j in range(i+2, len(coordinates)):
            if close_to_one(math.dist(coordinates[i], coordinates[j])):
                interaction_coords.append((coordinates[i], coordinates[j]))

    for i in range(len(coordinates) - 1):
        j = i + 1
        if close_to_one(math.dist(coordinates[i], coordinates[j])):
            connection_coords.append((coordinates[i], coordinates[j]))

    return (interaction_coords, connection_coords)


def close_to_one(num):
    if 0.98 < num < 1.02:
        return True
    return False


def visualize_4(sequence, energy_model, binary_output):
    # Visualizes the protein sequence on a square lattice
    moves = {
        '01': [1, 0],
        '10': [-1, 0],
        '11': [0, 1],
        '00': [0, -1]
    }

    binary = '01' + binary_output[0] + '1' + binary_output[1:]
    coordinates = generate_carteisan_coordinates(binary, moves, 2)
    color_sequences = encoded_sequence(energy_model, sequence)
    interaction_coords, connection_coords = get_interaction_coords(coordinates)
    plot_lattice(coordinates, color_sequences, interaction_coords, connection_coords)


def visualize_6(sequence, energy_model, binary_output):
    # Visualizes the protein sequence on a cubic lattice
    moves = {
        '000': [1, 0, 0],
        '001': [-1, 0, 0],
        '010': [0, 1, 0],
        '011': [0, -1, 0],
        '100': [0, 0, 1],
        '101': [0, 0, -1]
    }

    binary = '000' + binary_output[0] + '00' + binary_output[1:]
    coordinates = generate_carteisan_coordinates(binary, moves, 3)
    color_sequences = encoded_sequence(energy_model, sequence)
    interaction_coords, connection_coords = get_interaction_coords(coordinates)
    plot_lattice(coordinates, color_sequences, interaction_coords, connection_coords)


def visualize_8(sequence, energy_model, binary_output):
    # Visualizes the protein sequence on a trianguler prismatic lattice
    sqrt3_2 = math.sqrt(3) / 2
    moves = {
        '111': [0, 1, 0],
        '000': [0, -1, 0],
        '011': [sqrt3_2, 1/2, 0],
        '100': [-1 * sqrt3_2, -1/2, 0],
        '101': [-1 * sqrt3_2, 1/2, 0],
        '010': [sqrt3_2, -1/2, 0],
        '001': [0, 0, 1],
        '110': [0, 0, -1]
    }

    binary = '00' + binary_output
    coordinates = generate_carteisan_coordinates(binary, moves, 3)
    color_sequences = encoded_sequence(energy_model, sequence)
    interaction_coords, connection_coords = get_interaction_coords(coordinates)
    plot_lattice(coordinates, color_sequences, interaction_coords, connection_coords)


def visualize_12(sequence, energy_model, binary_output):
    # Visualizes the protein sequence on a FCC lattice
    sqrt2_2 = math.sqrt(2) / 2
    moves = {
        '1011': [sqrt2_2, sqrt2_2, 0],
        '1111': [-1 * sqrt2_2, sqrt2_2, 0],
        '1010': [sqrt2_2, -1 * sqrt2_2, 0],
        '1110': [-1 * sqrt2_2, -1 * sqrt2_2, 0],
        '0111': [0, sqrt2_2, sqrt2_2],
        '0101': [0, -1 * sqrt2_2, sqrt2_2],
        '0110': [0, sqrt2_2, -1 * sqrt2_2],
        '0100': [0, -1 * sqrt2_2, -1 * sqrt2_2],
        '1001': [sqrt2_2, 0, sqrt2_2],
        '1101': [-1 * sqrt2_2, 0, sqrt2_2],
        '1000': [sqrt2_2, 0, -1 * sqrt2_2],
        '1100': [-1 * sqrt2_2, 0, -1 * sqrt2_2]
    }

    binary = '10111' + binary_output[0] + binary_output[1] + '1' + binary_output[2:]
    coordinates = generate_carteisan_coordinates(binary, moves, 4)
    color_sequences = encoded_sequence(energy_model, sequence)
    interaction_coords, connection_coords = get_interaction_coords(coordinates)
    plot_lattice(coordinates, color_sequences, interaction_coords, connection_coords)


def generate_carteisan_coordinates(binary, moves, bits_per_move):
    if bits_per_move == 2:
        # Generates the cartesian coordinates of the protein sequence
        x, y = 0, 0
        coordinates = [(x, y)]

        for i in range(0, len(binary), 2):
            move = binary[i:i + 2]
            x += moves[move][0]
            y += moves[move][1]
            coordinates.append((x, y))
    else:
        # Generates the cartesian coordinates of the protein sequence
        x, y, z = 0, 0, 0
        coordinates = [(x, y, z)]

        for i in range(0, len(binary), bits_per_move):
            move = binary[i:i+bits_per_move]
            x += moves[move][0]
            y += moves[move][1]
            z += moves[move][2]
            coordinates.append((x, y, z))

    return coordinates


def encoded_sequence(energy_model, sequence):
    if energy_model == 'HP':
        encode_sequence = encode_hp(sequence)
    elif energy_model == 'HPAB':
        encode_sequence = encode_hpab(sequence)
    elif energy_model == 'WHPAB':
        encode_sequence = encode_hpab(sequence)
    else:
        encode_sequence = sequence

    return encode_sequence


def plot_lattice(coordinates, color_sequences, interaction_coords, connection_coords):
    num_of_dimensions = len(coordinates[0])

    # Map unique characters to colors
    unique_chars = sorted(list(set(color_sequences)))
    predefined_colors = ['red', 'blue', 'green', 'yellow', 'cyan']
    color_map = {}

    for idx, char in enumerate(unique_chars):
        color_map[char] = predefined_colors[idx % len(predefined_colors)]

    # Create a list of colors corresponding to color_sequences
    node_colors = [color_map[char] for char in color_sequences]

    if num_of_dimensions == 2:
        # Create a 2D plot
        fig, ax = plt.subplots()
        # Unzip the coordinates
        x, y = zip(*coordinates)
        # Plot the nodes with colors
        ax.scatter(x, y, s=100, c=node_colors, zorder=3)

        # Plot interactions first
        for coord_pair in interaction_coords:
            x_pair, y_pair = zip(*coord_pair)
            ax.plot(x_pair, y_pair, c='orange', linewidth=3, linestyle='--', zorder=2, label='Interaction')

        # Plot connections after interactions
        for coord_pair in connection_coords:
            x_pair, y_pair = zip(*coord_pair)
            ax.plot(x_pair, y_pair, c='grey', linewidth=3, zorder=1)

        # Set labels and title
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_title(color_sequences + ' Sequence on 2D Lattice')

        # Add legend for node types
        if len(unique_chars) <= 5:
            handles = [mpatches.Patch(color=color_map[char], label=char) for char in unique_chars]
            ax.legend(handles=handles, title='Amino Acid Types')

    elif num_of_dimensions == 3:
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Unzip the coordinates
        x, y, z = zip(*coordinates)
        # Plot the nodes with colors
        ax.scatter(x, y, z, s=100, c=node_colors, zorder=3)

        # Plot interactions first
        for coord_pair in interaction_coords:
            x_pair, y_pair, z_pair = zip(*coord_pair)
            ax.plot(x_pair, y_pair, z_pair, c='orange', linewidth=3, linestyle='--', zorder=2, label='Interaction')

        # Plot connections after interactions
        for coord_pair in connection_coords:
            x_pair, y_pair, z_pair = zip(*coord_pair)
            ax.plot(x_pair, y_pair, z_pair, c='grey', linewidth=3, zorder=1)

        # Set labels and title
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        ax.set_title(color_sequences + ' Sequence on 3D Lattice')

        # Add legend for node types
        if len(unique_chars) <= 5:
            handles = [mpatches.Patch(color=color_map[char], label=char) for char in unique_chars]
            ax.legend(handles=handles, title='Amino Acid Types')

    else:
        raise ValueError("Coordinates must be either 2D or 3D")

    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='QuantumHoneycombPSP CLI')
    parser.add_argument('sequence', type=str, help='Protein sequence')
    parser.add_argument('energy_model', type=str, choices=['HP', 'HPAB', 'WHPAB', 'MJ'], help='Energy model')
    parser.add_argument('lattice_type', type=int, choices=[4, 6, 8, 12], help='Lattice type')
    parser.add_argument('binary_output', type=str, help='Binary move output from QUBO')


    args = parser.parse_args()
    visualize(args.sequence, args.energy_model, args.lattice_type, args.binary_output)
