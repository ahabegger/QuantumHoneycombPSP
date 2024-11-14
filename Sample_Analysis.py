import pickle
import numpy as np

def sample_analysis(samples):
    bit_samples = samples.record.sample

    print(len(samples.variables))
    columns_to_delete = []

    for x in range(len(samples.variables)):
        if '*' in samples.variables[x]:
            columns_to_delete.append(x)

    print(len(columns_to_delete))

    for element in reversed(columns_to_delete):
        bit_samples = np.delete(bit_samples, element, 1)

    bitstring_samples = []
    for bit_sample in bit_samples:
        bitstring_samples.append(''.join([str(int(bit)) for bit in bit_sample]))

    # Assuming bit_samples, samples.record.num_occurrences, and samples.record.energy are lists or arrays of the same length
    for bit_sample, occurrence, energy in zip(bitstring_samples, samples.record.num_occurrences, samples.record.energy):
        print(bit_sample, occurrence, energy)

    # Create Matplot lib graph x-axis is bitstring_samples, y-axis is samples.record.num_occurrences
    import matplotlib.pyplot as plt

    # Plot for Number of Occurrences
    plt.bar(bitstring_samples, samples.record.num_occurrences)
    plt.xlabel('Bitstring Samples')
    plt.ylabel('Number of Occurrences')
    plt.title('Bitstring Samples vs Number of Occurrences')
    plt.xticks(rotation=90, fontsize=8)  # Rotate labels and adjust font size
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.show()

    # Plot for Energy
    plt.bar(bitstring_samples, samples.record.energy)
    plt.xlabel('Bitstring Samples')
    plt.ylabel('Energy')
    plt.title('Bitstring Samples vs Energy')
    plt.xticks(rotation=90, fontsize=8)  # Rotate labels and adjust font size
    plt.tight_layout()  # Adjust layout to prevent clipping
    plt.show()



def unpickle_file(file_path):
  with open(file_path, 'rb') as f:
    return pickle.load(f)

if __name__ == '__main__':
    samples = unpickle_file('Results/Samples_1731334129.pkl')

    sample_analysis(samples)

