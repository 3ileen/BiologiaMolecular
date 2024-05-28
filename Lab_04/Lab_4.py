
import numpy as np
import pandas as pd

def couple(pair):
    # Return True if RNA nucleotides are Watson-Crick base pairs

    pairs = {"A": "U", "U": "A", "G": "C", "C": "G"}

    # check if pair is couplable
    if pair in pairs.items():
        return True

    return False


def fill(nm, rna):

    # Fill the matrix as per the Nussinov algorithm

    minimal_loop_length = 0

    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k

            if j - i >= minimal_loop_length:
                down = nm[i + 1][j]  # 1st rule
                left = nm[i][j - 1]  # 2nd rule
                diag = nm[i + 1][j - 1] + couple((rna[i], rna[j]))  # 3rd rule

                rc = max([nm[i][t] + nm[t + 1][j] for t in range(i, j)])  # 4th rule

                nm[i][j] = max(down, left, diag, rc)  # max of all

            else:
                nm[i][j] = 0

    return nm


def traceback(nm, rna, fold, i, L):

    j = L

    if i < j:
        if nm[i][j] == nm[i + 1][j]:  # 1st rule
            traceback(nm, rna, fold, i + 1, j)
        elif nm[i][j] == nm[i][j - 1]:  # 2nd rule
            traceback(nm, rna, fold, i, j - 1)
        elif nm[i][j] == nm[i + 1][j - 1] + couple((rna[i], rna[j])):  # 3rd rule
            fold.append((i, j))
            traceback(nm, rna, fold, i + 1, j - 1)
        else:
            for k in range(i + 1, j - 1):
                if nm[i][j] == nm[i, k] + nm[k + 1][j]:  # 4th rule
                    traceback(nm, rna, fold, i, k)
                    traceback(nm, rna, fold, k + 1, j)
                    break

    return fold


def dot_write(rna, fold):
    dot = ["." for i in range(len(rna))]

    for s in fold:
        # print(min(s), max(s))
        dot[min(s)] = "("
        dot[max(s)] = ")"

    return "".join(dot)


def init_matrix(rna):
    M = len(rna)

    # init matrix
    nm = np.empty([M, M])
    nm[:] = np.NAN

    # init diaganols to 0
    nm[range(M), range(M)] = 0
    nm[range(1, len(rna)), range(len(rna) - 1)] = 0

    return nm


if __name__ == "__main__":
    #rna = "ACUCGAUUCCGAG"
    rna = "GGAAAUCC"

    nm = init_matrix(rna)
    nm = fill(nm, rna)

    fold = []
    sec = traceback(nm, rna, fold, 0, len(rna) - 1)

    res = dot_write(rna, fold)

    names = [_ for _ in rna]
    df = pd.DataFrame(nm, index=names, columns=names)
    print(df, "\n", rna, res)