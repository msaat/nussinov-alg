import numpy as np

base_index = {
    "A": 0,
    "C": 1,
    "G": 2,
    "U": 3
}

bond_scores = np.array([[-10, -10, -10, 2],
                       [-10, -10, 3, -10],
                       [-10, 3, -10, 1],
                       [2, -10, 1, -10]])

seq = input("Enter the RNA sequence\n")

seq_len = len(seq)

alignment_matrix = np.zeros((seq_len, seq_len), dtype=int)

traceback_indeces = []

def score(index_i, index_j):
    score = bond_scores[base_index[seq[index_i]], base_index[seq[index_j]]] + alignment_matrix[index_i + 1, index_j - 1]

    k = index_i
    while k < index_j:
        value = alignment_matrix[index_i, k] + alignment_matrix[k + 1, index_j]
        if value > score:
            score = value
        k += 1

    return score

def traceback(matrix, seq_len, index_i, index_j):
    traceback_indeces.append([index_i, index_j])
    if index_j > index_i:
        if matrix[index_i, index_j] == matrix[index_i, index_j - 1]:
            traceback(matrix, seq_len, index_i, index_j - 1)
        elif matrix[index_i, index_j] == matrix[index_i + 1, index_j]:
            traceback(matrix, seq_len, index_i + 1, index_j)
        else:
            traceback(matrix, seq_len, index_i + 1, index_j - 1)
            

#iterate diagonal by diagonal
i = 1
while i < seq_len:
    j = 0
    while j < seq_len - i:
        # calculate algorithm score
        alignment_matrix[j, i+j] = score(j, i+j)
        j += 1
    i += 1


print("The final score matrix is:")
print(alignment_matrix)
traceback(alignment_matrix, seq_len, 0, seq_len - 1)
print("The indeces for the traceback are:")
print(traceback_indeces)




