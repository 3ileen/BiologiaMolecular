import numpy as np
def scores(m_scores, m_dir, i, j, seq_max, seq_min):
    if seq_max[i - 1] == seq_min[j - 1]:
        a = [1 + m_scores[i - 1][j - 1], "d"]
    else:
        a = [m_scores[i - 1][j - 1] - 1,"d"]
    b = [m_scores[i - 1][j] - 2,"l"] # GAP
    c = [m_scores[i][j - 1] - 2, "u"] # GAP
    list_score = [a, b, c]
    max_score = max(list_score, key=lambda x: x[0])
    m_scores[i][j]= max_score[0]
    list_dir=[]
    for score, direction in list_score:
        if score == max_score[0]:
            list_dir.append(direction)
    m_dir[i][j] = list_dir


def get_alignments(seq_max, seq_min, m_dir, max_len, min_len):
    alignments = [[[] for j in range(min_len+1)] for i in range(max_len+1)]
    alignments[0][0] = [("", "")]

    for i in range(1, max_len+1):
        alignments[i][0] = [(seq_max[:i], "-" * i)]

    for j in range(1, min_len+1):
        alignments[0][j] = [("-" * j, seq_min[:j])]

    for i in range(1, max_len+1):
        for j in range(1, min_len+1):
            if "d" in m_dir[i][j]:
                for alignment in alignments[i-1][j-1]:
                    alignments[i][j].append((alignment[0] + seq_max[i-1], alignment[1] + seq_min[j-1]))
            if "l" in m_dir[i][j]:
                for alignment in alignments[i-1][j]:
                    alignments[i][j].append((alignment[0] + seq_max[i-1], alignment[1] + "-"))
            if "u" in m_dir[i][j]:
                for alignment in alignments[i][j-1]:
                    alignments[i][j].append((alignment[0] + "-", alignment[1] + seq_min[j-1]))

    return alignments[max_len][min_len]

    
def get_results_global(seq1, seq2):
    # Comparar cual secuencia es la mas larga
    if len(seq1) > len(seq2):
        max_len = len(seq1)
        min_len = len(seq2)
        seq_max = seq1
        seq_min = seq2
    else:
        max_len = len(seq2)
        min_len = len(seq1)
        seq_max = seq2
        seq_min = seq1

    # Matriz de scores
    m_scores = np.zeros((max_len + 1, min_len + 1)) - np.inf

    # Completar la primera columna y la primera fila de acuerdo al valor del GAP=-2
    m_scores[0, :] = np.arange(0, -(min_len * 2) - 2, -2)
    m_scores[:, 0] = np.arange(0, -(max_len * 2) - 2, -2)

    # Matriz de direcciones
    m_dir = np.empty((len(seq_max) + 2, len(seq_min) + 2), dtype=object)
    # Completar la primera columna y la primera fila
    for i in range(1, max_len + 1):
        m_dir[i][0] = "l"
    for j in range(1, min_len + 1):
        m_dir[0][j] = "u"

    # Completar la matriz de scores y direcciones
    for i in range(1, max_len + 1):
        for j in range(1, min_len + 1):
            scores(m_scores, m_dir, i, j, seq_max, seq_min)

    # Obtener los alineamientos como strings
    alignments = get_alignments(seq_max, seq_min, m_dir, max_len, min_len)

    return m_scores[max_len][min_len], len(alignments), alignments, m_scores, m_dir

def create_match_matrix(seq1, seq2):
    match_matrix = []
    for i in range(len(seq1)):
        row = []
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                row.append('*')
            else:
                row.append(' ')
        match_matrix.append(row)
    return match_matrix



#seq1 = "AAAC"
#seq2 = "AGC"
seq1 = "tcaagcgttagagaagtcattatgtgataa"
seq2 = "atggaagcaatatcactgatgactatacta "


match_matrix = create_match_matrix(seq1, seq2)
score, num_alignments, alignments, m_scores, m_dir = get_results_global(seq1, seq2)


print("Secuencia 1:", seq1)
print("Secuencia 2:", seq2)
print("Matriz de Scores:")
print(m_scores)


print("Matriz de coincidencias:")
for row in match_matrix:
    print(' '.join(row))

print("\nMatriz de Direcciones:")
print(m_dir)

print("\nScore:", score)
print("Número de alineamientos:", num_alignments)
print("Alineamientos:")
for alignment in alignments:
    print(alignment)
