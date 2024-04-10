from enum import Enum
import numpy as np
from functools import cmp_to_key

class Direction(Enum):
    UP = 1
    DIAGONAL = 2
    LEFT = 3

class Alignment:
    def __init__(self, sequenceA, sequenceB, score, peso):
        self.sequenceA = sequenceA
        self.sequenceB = sequenceB
        self.score = score
        self.peso = peso

def compareByPeso(a, b):
    return a.peso - b.peso

def Max(arriba, esquina, izquierda):
    return max(arriba, esquina, izquierda)

def Needleman(CadenaA, tamA, CadenaB, tamB):
    tamA += 1
    tamB += 1

    Matriz = np.zeros((tamA, tamB), dtype=int)
    Direcciones = [[[] for _ in range(tamB)] for _ in range(tamA)]

    menos_hor = -2
    menos_vert = -2

    for i in range(tamA):
        for j in range(tamB):
            if i == 0 and j == 0:
                Matriz[i][j] = 0
            elif i == 0 and j != 0:
                Matriz[i][j] = menos_hor
                menos_hor -= 2
                Direcciones[i][j].append(Direction.LEFT)
            elif i != 0 and j == 0:
                Matriz[i][j] = menos_vert
                menos_vert -= 2
                Direcciones[i][j].append(Direction.UP)

    for i in range(1, tamA):
        for j in range(1, tamB):
            match = 1 if CadenaA[i-1] == CadenaB[j-1] else -1
            diagonalScore = Matriz[i - 1][j - 1] + match
            upScore = Matriz[i - 1][j] - 2
            leftScore = Matriz[i][j - 1] - 2

            maxScore = Max(upScore, diagonalScore, leftScore)

            if maxScore == upScore:
                Direcciones[i][j].append(Direction.UP)
            if maxScore == diagonalScore:
                Direcciones[i][j].append(Direction.DIAGONAL)
            if maxScore == leftScore:
                Direcciones[i][j].append(Direction.LEFT)

            Matriz[i][j] = maxScore

    alignments = []
    
    def backtrack(seqA, seqB, row, col, newSeqA, newSeqB, score, peso):
        nonlocal alignments
        if row == 0 and col == 0:
            alignments.append(Alignment(seqA[::-1], seqB[::-1], score, peso))
            return

        for dir in Direcciones[row][col]:
            updatedSeqA = newSeqA
            updatedSeqB = newSeqB
            newRow = row
            newCol = col
            newScore = score
            newPeso = peso

            if dir == Direction.DIAGONAL:
                updatedSeqA += CadenaA[row-1]
                updatedSeqB += CadenaB[col-1]
                newRow -= 1
                newCol -= 1
                newScore += 1 if CadenaA[row-1] == CadenaB[col-1] else -1
            elif dir == Direction.UP:
                updatedSeqA += CadenaA[row-1]
                updatedSeqB += '-'
                newRow -= 1
                newScore -= 2
            elif dir == Direction.LEFT:
                updatedSeqA += '-'
                updatedSeqB += CadenaB[col-1]
                newCol -= 1
                newScore -= 2

            newPeso += Matriz[row][col]

            backtrack(seqA + (updatedSeqA[-1] if updatedSeqA[-1] != '-' else ""), 
                      seqB + (updatedSeqB[-1] if updatedSeqB[-1] != '-' else ""),
                      newRow, newCol, updatedSeqA, updatedSeqB, newScore, newPeso)

    backtrack("", "", tamA-1, tamB-1, "", "", 0, 0)

    alignments.sort(key=cmp_to_key(compareByPeso))

    return alignments[0]

def MSA_Estrella(cadenas):
    n = len(cadenas)
    Matriz_MSA = np.zeros((n, n), dtype=int)

    for i in range(n-1):
        for j in range(i+1, n):
            secuenciasAlineadas = Needleman(cadenas[i], len(cadenas[i]), cadenas[j], len(cadenas[j]))
            Matriz_MSA[i][j] = secuenciasAlineadas.score
            Matriz_MSA[j][i] = secuenciasAlineadas.score

    print("Matriz de scores:")
    print(Matriz_MSA)

    print("\nSuma de scores:")
    scores_suma = [sum(row) for row in Matriz_MSA]
    max_score_idx = scores_suma.index(max(scores_suma))
    print(scores_suma)
    print(max(scores_suma))
    print("Secuencia con mayor score:", cadenas[max_score_idx])

    SecuenciaOptimaLarga = ""

    print("\nCombinaciones finales MSA:")
    for j in range(n):
        if j != max_score_idx:
            secuenciasAlineadas = Needleman(cadenas[max_score_idx], len(cadenas[max_score_idx]), cadenas[j], len(cadenas[j]))
            print("Sec " + str(max_score_idx+1) + ": " + secuenciasAlineadas.sequenceA)
            print("Sec " + str(j+1) + ": " + secuenciasAlineadas.sequenceB)

            if len(SecuenciaOptimaLarga) < len(secuenciasAlineadas.sequenceA):
                SecuenciaOptimaLarga = secuenciasAlineadas.sequenceA

    print("\nSecuencia mas larga:")
    print(SecuenciaOptimaLarga)

    print("\nSec " + str(max_score_idx+1) + ": " + SecuenciaOptimaLarga)

    print("\nAlineaciones finales:")
    for j in range(n):
        if j != max_score_idx:
            secuenciasAlineadas = Needleman(SecuenciaOptimaLarga, len(SecuenciaOptimaLarga), cadenas[j], len(cadenas[j]))
            print("Sec " + str(j+1) + ": " + secuenciasAlineadas.sequenceB)

if __name__ == "__main__":
    secuencias = [
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    ]

    MSA_Estrella(secuencias)
