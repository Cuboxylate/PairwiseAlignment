from platform import python_version
print("python", python_version())

import numpy as np
print("numpy", np.version.version)

import pandas
print("pandas", pandas.__version__)

import numpy as np

# F matrix, sequence A in row, B in col
def get_F_matrix(A, B, match=1.0, mismatch=-1.0, gap=-2.0):
    # The dimension of the matrix
    rown = len(A) + 1
    coln = len(B) + 1

    # init 2D F matrix, shape=(len(A)+1, len(B)+1). The +1 is because we will have a row and column past the end
    # of the sequences, to represent the score for if the last bases are match/mismatch/gap
    F = np.zeros(shape=(rown, coln))

    # go through each row in the outer loop, then each column in each row in the inner loop
    for i in range(rown):
        for j in range(coln):

            if i == 0:
                # Fill the i == 0 column, representing leading gaps in B
                if j != 0:
                    # [0][0] is already initialised to 0
                    F[i][j] = F[i][j-1] + gap
                continue
            elif j == 0:
                # fill j == 0 row, representing leading gaps in A
                F[i][j] = F[i-1][j] + gap
                continue

            # Then fill each remaining cell with the highest score out of the three options:
            # - Match or Mismatch score + F[i-1][j-1]
            # - Inserting a gap in A
            # - Inserting a gap in B

            basesMatch = A[i-1] == B[j-1]
            previousCellScore = F[i-1][j-1]
            diagonalScore = previousCellScore + match if basesMatch else previousCellScore + mismatch

            gapInA = F[i][j-1] + gap
            gapInB = F[i-1][j] + gap

            F[i][j] = max(diagonalScore, gapInA, gapInB)

    # the starting position is bottom right element in F matrix
    return F

# F matrix, A in the row headers (i.e. down the left side), B in column headers (i.e. across the top)
def trace_back(F, A, B, match=1.0, mismatch=-1.0, gap=-2.0):
    # i indexing row, j indexing col in F
    print(F)
    i, j = len(A), len(B)
    if i == 0 or j == 0 or i >= np.shape(F)[0] or j >= np.shape(F)[1]:
        raise ValueError('Invalid position [%s, %s] !' % (i, j))

    # aligned sequence A
    AlignmentA = ""
    # aligned sequence B
    AlignmentB = ""

    """
    Trace back in F matrix following the possible moves by recalculating the score.
    The starting position is bottom right element in F matrix.
    """

    while (i != 0 and j != 0) :
        currentCellScore = F[i][j]

        # Find which direction to trace back through. The order of these conditions
        # enforces taking the high road when there are equally scored directions

        if (currentCellScore == F[i-1][j] + gap) :
            # if this score came from a previous gap in B (vertical cell):
            AlignmentA = A[i-1] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i = i-1
        elif ((A[i-1] == B[j-1] and currentCellScore == F[i-1][j-1] + match) or
              (A[i-1] != B[j-1] and currentCellScore == F[i-1][j-1] + mismatch)) :
            # if this score came from the diagonal cell:
            AlignmentA = A[i-1] + AlignmentA
            AlignmentB = B[j-1] + AlignmentB
            i = i-1
            j = j-1
        else :
            # if this score came from a previous gap in A (horizontal cell):
            AlignmentA = "-" + AlignmentA
            AlignmentB = B[j-1]  + AlignmentB
            j = j-1


    # string type
    return AlignmentA, AlignmentB


# https://en.wikipedia.org/wiki/Needleman–Wunsch_algorithm
x = "GATTACA"
y = "GCATGCU"
mat = get_F_matrix(x, y, gap=-1)
alignment = trace_back(mat, x, y, gap=-1)
print("y: ", alignment[1], "\nx: ", alignment[0])
print(mat)

assert alignment[1] == "GCATG-CU"
assert alignment[0] == "G-ATTACA"
assert np.array_equal(mat,
                      np.reshape((0, -1, -2, -3, -4, -5, -6, -7,
                                  -1,  1,  0, -1, -2, -3, -4, -5,
                                  -2,  0,  0,  1,  0, -1, -2, -3,
                                  -3, -1, -1,  0,  2,  1,  0, -1,
                                  -4, -2, -2, -1,  1,  1,  0, -1,
                                  -5, -3, -3, -1,  0,  0,  0, -1,
                                  -6, -4, -2, -2, -1, -1,  1,  0,
                                  -7, -5, -3, -1, -2, -2,  0,  0), (8, 8)) )
