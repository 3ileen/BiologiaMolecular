
sequences = [
    "VGA--HAGEY",
    "V----NVDEV",
    "VEA--DVAGH",
    "VKG------D",
    "VYS--TYETS",
    "FNA--NIPKH",
    "IAGADNGAGY"
]

# calcular la cantidad de estados (columnas) donde más de la mitad de los caracteres son residuos
num_columns = len(sequences[0])
num_sequences = len(sequences)
min_residues = num_sequences // 2 + 1  # Más de la mitad
states_count = sum(1 for i in range(num_columns) if sum(s[i] != '-' for s in sequences) >= min_residues)

# Calculamos las probabilidades de emisión
amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # Lista estándar de aminoácidos
emission_probabilities = {aa: [0] * num_columns for aa in amino_acids}

# frecuencias
for col in range(num_columns):
    total_residues = 0
    for seq in sequences:
        aa = seq[col]
        if aa != '-':
            emission_probabilities[aa][col] += 1
            total_residues += 1
    # Normalizar las frecuencias para convertirlas en probabilidades
    if total_residues > 0:
        for aa in amino_acids:
            emission_probabilities[aa][col] /= total_residues

# Probabilidades de transición y background
# distribución equiprobable para el background
background_probabilities = {aa: 1/len(amino_acids) for aa in amino_acids}

# Probabilidades de Transición
# Inicializamos los contadores para las transiciones
transition_counts = {
    "M-M": 0, "M-I": 0, "M-D": 0,
    "I-M": 0, "I-I": 0,
    "D-M": 0, "D-D": 0
}

# Función para determinar el estado de una columna (Match, Insert, Delete)
def column_state(column):
    if '-' in column and any(residue != '-' for residue in column):
        return "I"  # Insert
    elif all(residue == '-' for residue in column):
        return "D"  # Delete
    else:
        return "M"  # Match

# Calcular las transiciones entre columnas
for col in range(num_columns - 1):
    current_state = column_state([seq[col] for seq in sequences])
    next_state = column_state([seq[col + 1] for seq in sequences])

    transition_key = f"{current_state}-{next_state}"
    transition_counts[transition_key] += 1

# Aplicamos suavizado para manejar frecuencias cero
for key in transition_counts:
    transition_counts[key] += 1

# Calcular probabilidades de transición
total_transitions = sum(transition_counts.values())
transition_probabilities = {key: count / total_transitions for key, count in transition_counts.items()}




print("\nCantidad de Estados:")
print(states_count)

print("\nProbabilidades de Emisión:")
for aa in amino_acids:
    print(f"{aa}: {['{:.3f}'.format(prob) for prob in emission_probabilities[aa]]}")

print("\nProbabilidades de Transición:")
for transition, prob in transition_probabilities.items():
    print(f"{transition}: {prob:.3f}")

print("\nVector de Probabilidades de Fondo:")
for aa, prob in background_probabilities.items():
    print(f"{aa}: {prob:.3f}")