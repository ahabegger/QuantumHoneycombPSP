# QuantumHoneycombPSP

QuantumHoneycombPSP is a Python-based project designed to model protein structure prediction (PSP) using quantum computing techniques. The project includes various modules to encode sequences, generate energy matrices, and create energy functions for different lattice types. It also supports the conversion of these energy functions into Quadratic Unconstrained Binary Optimization (QUBO) and Ising models.

## Features

- **Sequence Encoding**: Encode protein sequences using different models such as HP, HPAB, and WHPAB.
- **Energy Matrix Generation**: Generate energy matrices based on the encoded sequences and selected energy models.
- **Energy Function Creation**: Create energy functions for different lattice types (4, 6, 8, 12) and convert them into QUBO and Ising models.
- **Interaction Calculation**: Calculate interactions between amino acids in the sequence and generate corresponding energy values.
- **QUBO Variable Mapping**: Map QUBO variables to a new format for further processing.
- **Structure Imager**: Generate images of protein structures based on the encoded sequences.

## Installation

To install the required dependencies, run:

```bash
pip install -r requirements.txt
```

## Usage

### CLI

The main.py script can be used as a command-line interface (CLI) to specify the sequence, energy model, lattice type, and whether to use the binary model.
Command-Line Arguments

- sequence: Protein sequence (e.g., 'GAAGA')
- energy_model: Energy model (choices: 'HP', 'HPAB', 'WHPAB', 'MJ')
- lattice_type: Lattice type (choices: 4, 6, 8, 12)
- --binary: Use the binary model (optional flag)

### Example Command

```bash
python main.py GAAGA HP 4 --binary
```

## Output

The script will output the following information:

- Encoded sequence
- Interaction matrix
- QUBO, BQM, and Ising models (if --binary is not used)
- Mapped QUBO variables
- Final QUBO equation

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
