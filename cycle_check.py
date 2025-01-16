import sys
import numpy as np

# Function to parse the input file
def parse_input_file(file_path):
    permutations = {}
    cycles = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith("Permutations["):
                name, values = line.split("=")
                name = name.strip()
                #print(f'The name is {name}')
                values = list(map(int, values.strip().split()))
                permutations[name] = values

            elif line.startswith("cycles["):
                _, values = line.split("=")
                values = list(map(int, values.strip().split()))
                cycles.append(values)

    return permutations, cycles

# Function to validate cycles against permutations
def validate_cycles(permutations, cycles):
    permutation_keys = list(permutations.keys())
    num_permutations = len(permutation_keys)
    invalid_cycles = []
    invalid_cycles_sum = []

    print(f'The number of permutations is {num_permutations}')
    print(' ')

    for i, cycle in enumerate(cycles):
        if len(cycle) != num_permutations:
            print(f"Error: Cycle[{i}] length does not match the number of permutations.")
            continue

        final_array = [0] * len(permutations[permutation_keys[0]])

        for index, multiplier in enumerate(cycle):
            permutation_array = permutations[permutation_keys[index]]
            final_array = [
                final_array[j] + multiplier * permutation_array[j]
                for j in range(len(permutation_array))
            ]

        if any(final_array):
            invalid_cycles.append(i)
            invalid_cycles_sum.append(np.sum(final_array))

    return invalid_cycles , invalid_cycles_sum

def get_cycle_indices(cycle):
    perm_indices = []
    for i in range(len(cycle)):
        if abs(cycle[i]) > 1E-7:
            perm_indices.append(i)
    
    return perm_indices

# Main script
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]

    # Parse the input file
    permutations, cycles = parse_input_file(input_file)

    # print(f'Cycles are {cycles}')
    print(f'The permutation keys are {permutations.keys()}')
    # Validate the cycles
    invalid_cycles , invalid_cycles_sum = validate_cycles(permutations, cycles)
    cycle_indices = []

    for inv_cycle_num in invalid_cycles:
        cycle_indices.append(get_cycle_indices(cycles[inv_cycle_num]))

    # Report results
    if invalid_cycles:
        print(f"The following cycles did not sum to all zeros: {invalid_cycles}")
        print(f"The sum of the invalid cycles is: {invalid_cycles_sum}")
        print(' ')
        print(f'The perm indices of the invalid cycles are {cycle_indices}')
    else:
        print("All cycles summed to all zeros.")
