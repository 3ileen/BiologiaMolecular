def find_consensus_sequence(sequences):
    max_length = 0
    for sequence in sequences:
        max_length = max(max_length, len(sequence))

    consensus = ''
    for i in range(max_length):
        nucleotide_count = {}
        for sequence in sequences:
            if i < len(sequence):
                nucleotide_count[sequence[i]] = nucleotide_count.get(sequence[i], 0) + 1

        most_frequent = ' '
        max_count = 0
        for nucleotide, count in nucleotide_count.items():
            if count > max_count:
                max_count = count
                most_frequent = nucleotide

        consensus += most_frequent
    return consensus


def reverse_complement(sequence):
    complement = ''
    for nucleotide in sequence:
        if nucleotide == 'A':
            complement += 'T'
        elif nucleotide == 'T':
            complement += 'A'
        elif nucleotide == 'C':
            complement += 'G'
        elif nucleotide == 'G':
            complement += 'C'
    return complement[::-1]


def hamiltonian_paths(sequences, min_overlap, length):
    combined_sequences = []
    for index in range(len(sequences)):
        combined_sequence = sequences[index]
        contributing_sequences = [sequences[index]]
        for i in range(len(sequences)):
            if i == index:
                continue
            overlap = 0
            len1 = len(combined_sequence)
            len2 = len(sequences[i])

            for j in range(min_overlap, min(len1, len2) + 1):
                if combined_sequence[len1 - j:] == sequences[i][:j]:
                    overlap = j
            if overlap < min_overlap:
                continue
            combined_sequence += sequences[i][overlap:]
            contributing_sequences.append(sequences[i])
        print("\nCadena", index, ": ")
        for sequence in contributing_sequences:
            print(sequence, "->", end=' ')
        contributing_sequences.append(combined_sequence)
        combined_sequences.append(contributing_sequences)

    best_index = 0
    for i in range(1, len(combined_sequences)):
        if len(combined_sequences[best_index][-1]) < len(combined_sequences[i][-1]):
            best_index = i

    return combined_sequences[best_index]


sequences = [
    "ATCCGTTGAAGCCGCGGGC",
    "TTAACTCGAGG",
    "TTAAGTACTGCCCG",
    "ATCTGTGTCGGG",
    "CGACTCCCGACACA",
    "CACAGATCCGTTGAAGCCGCGGG",
    "CTCGAGTTAAGTA",
    "CGCGGGCAGTACTT"
]

consensus_sequence = find_consensus_sequence(sequences)
print("Secuencia de Consenso:", consensus_sequence)
print("------------------------------------------------------------")

all_sequences = sequences + [reverse_complement(sequence) for sequence in sequences]

desired_length = 55
min_overlap = 1
results = hamiltonian_paths(all_sequences, min_overlap, desired_length)

print("\n------------------------------------------------------------")
print("Secuencia proxima a longitud", desired_length, ":", results[-1])
print("Tamaño:", len(results[-1]))
print("Camino Hamiltoneano:")
for sequence in results:
    print(sequence, "->", end=' ')
print()
