simulated_annealing.py will take sequences of 'p' and 'h' and optimize them with simulated annealing 20 times, then save a png file with the 20 structures obtained for each sequence. The files
are named with number, starting at 1.
To run the script, use the sequences you wish to optimize as command line arguments. An example command is:
python simulated_annealing.py hhhhhhhhhhhhhhhhhhhhppppppppppppppppppppppppppppppppp hpppppphhhhhhhhhhhhppphhhhhhhpppppppphhhhhhpppp hppppppppphhhhhhhhhhhhhhhppppppppphhhhhhpppphhhh
Pseudocode:

Function calculate_energy(sequence, path)
    Initialize energy to 0
    Initialize position as a list with two zeros
    Initialize positions as a list containing the position list
    Initialize directions as a list of tuples representing four directions

    Define a dictionary mapping from path to direction
    Set direction as the first element in directions

    For each residue and move in sequence and path
        Update the direction based on the move and current direction
        Update the position based on the direction

        If the updated position is already in positions
            Increase energy by 10^10
        Add the updated position to positions

        For each direction in directions
            Calculate the neighbor position
            If the neighbor position is in positions
                Get the index of the neighbor position in positions
                If the sequence at the index is 'h'
                    Decrease energy by 1 if residue is 'h'
                Else if the sequence at the index is 'p'
                    Increase energy by 1 if residue is 'h'

    Return energy

Function optimize_structure(sequences)
    Initialize energies and paths as empty lists

    For each sequence in sequences
        Initialize path as a list of 'f' with length of sequence minus 1

        For each iteration in range 30000
            Calculate temperature T
            Calculate initial energy Ei using calculate_energy function
            Copy path to new_path
            Randomly select a position in new_path and change its direction
            Calculate new energy Ej using calculate_energy function
            If Ej is less than Ei or a random number is less than the exponential of negative (Ej - Ei) divided by T
                Update path to new_path

        Add path to paths
        Add Ei to energies

    Return paths and energies

Function main
    Get sequences from command line arguments

    For each sequence in sequences
        Initialize all_paths and all_energies as empty lists

        For each iteration in range 20
            Get paths and energies by calling optimize_structure function with sequence
            Add paths to all_paths
            Add energies to all_energies

        Call render_structures function with all_paths, sequence repeated 20 times, all_energies, and index

If the script is run directly
    Call main function
