import pandas
import numpy as np

def readBLOSUM62():
    text_file = open("blosum62.txt", "r")

    # This method assumes the same number of tabs in each row of the .txt
    BLOSUM_df = pandas.read_csv(text_file, delimiter="\t", header=0, index_col=0)

    # pandas.DataFrame
    return BLOSUM_df

BLOSUM62 = readBLOSUM62()
print(list(BLOSUM62['A']))

# Find the overlap of greatest score. Will choose an overlap of A over B or B over A depending on which
# scores the highest.
def overlap(A,B,S,d=-2.0):

    numberRows = len(A)+1
    numberColumns = len(B)+1

    # init 2D F matrix, shape=(len(A)+1, len(B)+1)
    F = np.zeros((numberRows,numberColumns))

    # the highest value in the final row or column, and its position
    max_score = 0
    max_pos = None

    ''' 
    Returns a pairwise sequence alignment between sequences A and B,
    score is the the highest value at the F(i,m) or F(n,j) cell to start.
    S is the score matrix, and d is the penalty. 
    '''

    # go through each row in the outer loop, then each column in each row in the inner loop
    for i in range(numberRows):
        for j in range(numberColumns):

            if i == 0 or j == 0:
                # This column and row are already filled by 0s
                continue

            # Then fill each remaining cell with the highest score out of the three options:
            # - Match or mismatch with score from S. S is indexed by the amino acid's letter, taken from A and B
            # - Inserting a gap in A
            # - Inserting a gap in B

            previousCellScore = F[i-1][j-1]
            diagonalScore = previousCellScore + S.loc[A[i-1], B[j-1]]

            gapInA = F[i][j-1] + d
            gapInB = F[i-1][j] + d

            F[i][j] = max(diagonalScore, gapInA, gapInB)

            # Keep track of the max score calculated in the final row or column, and where it was
            if (i == len(A) or j == len(B)) and F[i][j] >= max_score :
                max_score = F[i][j]
                max_pos = (i, j)

    AlignmentA, AlignmentB = trace_back_overlap(F, A, B, S, max_pos, d)

    # AlignmentA or AlignmentB is string, score is a real number
    return AlignmentA, AlignmentB, max_score, F


# F is matrix of scores, A and B are the sequences, S is scoring matrix indexed by base, gap is gap penalty
def trace_back_overlap(F, A, B, S, startingCoords, gap):
    i, j = startingCoords[0], startingCoords[1]

    if i >= np.shape(F)[0] or j >= np.shape(F)[1]:
        raise ValueError('Invalid position [%s, %s] !' % (i, j))

    AlignmentA = ""
    AlignmentB = ""
    A_finished_first = False;

    if i == len(A):
        # This means A finishes before B in the alignment
        # Therefore B starts as the remainder of B left at j, and A as that many gaps
        AlignmentA = "-" * (len(B) - j)
        AlignmentB = B[j:]
        A_finished_first = True
    else :
        # B finishes before A. We initialise the alignment in the opposite way
        AlignmentA = A[i:]
        AlignmentB = "-" * (len(A) - i)

    # Loop until either i or j hits 0, based on whether B or A finishes first in the alignment
    while ((A_finished_first and j != 0) or i != 0) :
        currentCellScore = F[i][j]

        # Find which direction to trace back through. The order of these conditions
        # enforces taking the high road when there are equally scored directions

        if (currentCellScore == F[i-1][j] + gap):
            # if this score came from a previous gap in B (vertical cell):
            AlignmentA = A[i-1] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i = i-1
        elif (currentCellScore == F[i-1][j-1] + S.loc[A[i-1], B[j-1]]) :
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

        currentCellScore = F[i][j]

    # Now we have the end and overlapping parts of our alignment done, we just need to fill in the beginning
    if A_finished_first :
        # We need to put the rest of A in front of its alignment, and finish B with leading gaps
        AlignmentA = A[:i] + AlignmentA
        AlignmentB = "-" * i + AlignmentB
    else :
        # If B finishes first, we need to put leading gaps on AlightmentA and the front of B on AlignmentB
        AlignmentA = "-" * j + AlignmentA
        AlignmentB = B[:j] + AlignmentB

    # string type
    return AlignmentA, AlignmentB

x = 'PAWHEAE'
y = 'HEAGAWGHEE'

AlignmentA,AlignmentB,score,F = overlap(x, y, BLOSUM62, -8)
print(AlignmentB)
print(AlignmentA)
print(score)
print(F)

assert AlignmentB == "HEAGAWGHEE-"
assert AlignmentA == "---PAW-HEAE"
assert score == 17.0