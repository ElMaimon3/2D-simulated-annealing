import sys, os
import getopt
import math
import random
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

def render_structures(paths, sequences, energies, i):
    # Create a big figure
    fig, axs = plt.subplots(5, 4, figsize=(20, 20))  # Adjust as needed
    axs = axs.flatten()

    for j, (path, sequence) in enumerate(zip(paths, sequences)):
        # Initialize the position of the protein
        position = [0, 0]
        positions = [list(position)]
        directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]

        # Define the mapping from path to direction
        path_to_direction = {'f': 0, 'r': 1, 'l': 3}
        direction = directions[0]
        # Iterate over the sequence and the path
        for move in path:
            # Update the position
            direction = directions[(path_to_direction[move] + directions.index(direction)) % 4]
            position[0] += direction[0]
            position[1] += direction[1]
            positions.append(list(position))

        # Determine the size of the grid
        x_positions = [pos[0] for pos in positions]
        y_positions = [pos[1] for pos in positions]
        grid_size = max(max(x_positions) - min(x_positions), max(y_positions) - min(y_positions)) + 3  # Add 3 to leave a space around the protein

        # Draw the protein on the grid
        for pos, residue in zip(positions, sequence):
            x, y = pos[0] - min(x_positions) + 1, pos[1] - min(y_positions) + 1  # Add 1 to leave a space around the protein
            color = 'black' if residue == 'h' else 'red'  # Black for hydrophobic residues, Red for polar residues
            circle = patches.Circle((x, y), radius=0.5, facecolor=color)
            axs[j].add_patch(circle)

        # Mark the initial and final residues with a white dot
        circle = patches.Circle((positions[0][0] - min(x_positions) + 1, positions[0][1] - min(y_positions) + 1), radius=0.2, facecolor='white')
        axs[j].add_patch(circle)
        circle = patches.Circle((positions[-1][0] - min(x_positions) + 1, positions[-1][1] - min(y_positions) + 1), radius=0.2, facecolor='white')
        axs[j].add_patch(circle)

        # Set the limits and aspect ratio of the plot
        axs[j].set_xlim(0, grid_size)
        axs[j].set_ylim(0, grid_size)
        axs[j].set_aspect('equal', 'box')

        # Add a grid
        axs[j].grid(True, which='both', color='grey', linewidth=4)

        # Hide the axes
        axs[j].axis('off')

        # Add the sequence and path under the figure
        path = ''.join(path)
        axs[j].set_title(f"Path: {path}\nEnergy: {energies[j]}")

    # Save the figure
    plt.savefig(f'T10_9998_seq{i+1}.png', bbox_inches='tight')
    plt.close()

def calculate_energy(sequence, path):
    # Initialize the energy and the position of the protein
    energy = 0
    position = [0, 0]
    positions = [list(position)]
    directions = [(1, 0), (0, 1), (-1, 0), (0, -1)]

    # Define the mapping from path to direction
    path_to_direction = {'f': 0, 'r': 1, 'l': 3}
    direction = directions[0]
    # Iterate over the sequence and the path
    for residue, move in zip(sequence[:-1], path):
        # Update the position
        direction = directions[(path_to_direction[move] + directions.index(direction)) % 4]
        position[0] += direction[0]
        position[1] += direction[1]

        # Check for self-crossing
        if position in positions:
            energy += 10**10
        positions.append(list(position))

        # Calculate the energy
        for dx, dy in directions:
            neighbor = [position[0] + dx, position[1] + dy]
            if neighbor in positions:
                index = positions.index(neighbor)
                if sequence[index] == 'h':
                    energy -= 1 if residue == 'h' else 0
                elif sequence[index] == 'p':
                    energy += 1 if residue == 'h' else 0

    return energy


def optimize_structure(sequences):
    # Takes a list of sequences of Hs and Ps and returns the path reached with simulated annealing as a set of instructions
    energies = []
    paths = []
    for seq in sequences:
        path = ['f'] * (len(seq) - 1)  # start with a straight line path
        for _ in range(30000):  # number of iterations
            T = 10*0.9998**_
            Ei = calculate_energy(seq, path)
            new_path = path.copy()
            change_pos = random.randint(0, len(path) - 1)
            new_path[change_pos] = random.choice(['l', 'r', 'f'])
            Ej = calculate_energy(seq, new_path)
            if Ej < Ei or random.random() < math.exp(-(Ej - Ei) / T):
                path = new_path
        paths.append(path)
        energies.append(Ei)
    return paths,energies


def main():
    sequences = sys.argv[1:]
    for i, sequence in enumerate(sequences):
        all_paths = []
        all_energies = []
        for _ in range(20):
            paths, energies = optimize_structure([sequence])
            all_paths.extend(paths)
            all_energies.extend(energies)
        render_structures(all_paths, [sequence] * 20, all_energies, i)


if __name__=='__main__':
    main()
