import argparse

UNIT_COST_MATRIX = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
BASE_LIST = ["A", "C", "G", "T"]
string = "sss"

def main():
    parser = argparse.ArgumentParser(description='Calculate multiple alignment using the center star approximation')
    parser.add_argument('--input', metavar='file', type=str,
                        help='Input fasta file')
    parser.add_argument('--dist',metavar='file', action='store_const',
                        const=sum, default="unit costs",
                        help='File containing the distance matrix. If left out unit costs will be used.')

    args = parser.parse_args()


# calculates the edit distance for two DNA sequences and a
# score matrix, if no cost matrix is given, unit costs will be used
def edit_distance(seq1, seq2, mat=None):
    if mat is None:
        mat = UNIT_COST_MATRIX
    # initialize the first row of matrix
    row1 = [i for i in range(len(seq1) + 1)]
    # only two rows are saved of the matrix to save computational space
    row2 = [len(seq2) + 1]

    for i in range(len(seq2)):
        i = i + 1
        # initialize fist entry for recursive calculation
        row2[0] = i
        for j in range(len(seq1)):
            j = j + 1
            # implements the recursive definition, view sequence analysis script
            row2[j] = max(row1[j - 1] + get_score(seq1[j - 1], seq2[i - 1], mat), row1[j] + 1, row2[j - 1] + 1)
        row1 = row2
    return row1[len(row1) - 1]


# calculates score for 2 bases and a score matrix
def get_score(base1, base2, mat):
    index1 = BASE_LIST.index(base1)
    index2 = BASE_LIST.index(base2)
    return mat[index1][index2]


def add_sequence_to_alignment(alignment, seq):
    row1 = [i for i in range(len(alignment) + 1)]
    row2 = [len(alignment) + 1]
    hirschbergRowOne = get_hirschberg(row1, row2, alignment, seq)
    hirschbergRowTwo = get_hirschberg(row1.reverse(), row2, alignment.reverse(), seq.reverse())

def get_hirschberg(row1, row2, alignment, seq):
    limit = len()


if __name__ == '__main__':
    main()
