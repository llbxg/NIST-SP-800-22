import numpy as np

# Swap the rows to lock 1 to the left edge.
def swap(mat, start, m, M, Q):
    for q in range(start, Q):
        if mat[q][m]==1:
            if q != 0:
                mat[start], mat[q] = mat[q], mat[start]
            return True, mat
    return False, mat

# Clean to make the lower triangular matrix zero.
def clean(mat, start, m, M, Q):
    p_m = mat[start]
    for i in range(start+1, Q):
        if mat[i][m] == 1:
            mat[i] = p_m ^ mat[i]
    return mat

# Delete unnecessary 1 as much as possible.
def delete(mat, M, Q):
    for q in range(Q):
        last_mat = None
        for m in reversed(range(M)):
            if mat[m][q]==1:
                if last_mat is None:
                    last_mat = mat[m]
                else:
                    mat[m] = last_mat ^ mat[m]
    return mat

# Adjust to make a step.
def adjust(mat, M, Q):
    for q in reversed(range(Q)):
        for m in list(reversed(range(M)))[Q-q:]:
            if mat[m][q] == 1:
                mat[m], mat[q] = mat[q], mat[m]
                break
    return mat

# Compute the rank of a matrix of 0s and 1s.
def rank(mat, M = 32, Q = 32):
    for m in range(M):
        b, mat = swap(mat, m, m, M, Q)
        if b:
            mat = clean(mat, m, m, M, Q)

    mat = delete(mat, M, Q)
    mat = adjust(mat, M,Q)

    s = sum([1 if all([i == 0 for i in i]) else 0 for i in mat ])

    return  M-s