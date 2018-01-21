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
def edit_distance(seq1, seq2, mat=None, indel_cost = 1):
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
            row2[j] = max(row1[j - 1] + get_score(seq1[j - 1], seq2[i - 1], mat), row1[j] + indel_cost, row2[j - 1] + indel_cost)
        row1 = row2
    return row1[len(row1) - 1]


# calculates score for 2 bases and a score matrix
def get_score(base1, base2, mat):
    index1 = BASE_LIST.index(base1)
    index2 = BASE_LIST.index(base2)
    return mat[index1][index2]


def get_score_alignment(alignmentColumn, base, mat):
    totalScore = 0
    alignmentColumn += [base]
    for i in range(len(alignmentColumn)):
        for j in alignmentColumn:
            totalScore += get_score(j, base, mat)
        alignmentColumn.remove(alignmentColumn[0])
    return totalScore


def calculate_alignment(center_star, sequences, mat):
    multiple_alignment = [c for c in center_star]
    for i in range(len(sequences)):
        append_next_sequence(multiple_alignment, sequences[i], mat)
    return multiple_alignment


def append_next_sequence(multiple_alignment, seq, mat):
    new_multiple_alignment = [multiple_alignment]
    while True:
        coordinates = get_next_iteration_hirschberg(multiple_alignment, seq, mat)


def get_next_iteration_hirschberg(alignment, seq, mat):
    row1 = [i for i in range(len(alignment) + 1)]
    row2 = [len(alignment) + 1]
    hirschberg_row_one = get_hirschberg(row1, row2, alignment, seq, mat)
    row1 = [i for i in range(len(alignment) + 1)]
    row2 = [len(alignment) + 1]
    hirschberg_row_two = get_hirschberg(row1.reverse(), row2, alignment.reverse(), seq.reverse(), mat)
    row1 = [i for i in range(len(alignment) + 1)]
    row2 = [len(alignment) + 1]
    edit_dist = get_edit_distance_multiple(row1, row2, alignment, seq, mat)
    i_coordinate = 0
    j_coordinate = (len(seq)+1)/2 + ((len(seq) + 1) % 2)
    for i in range(len(hirschberg_row_one)):
        if (hirschberg_row_one[i] + hirschberg_row_two[i]) - edit_dist == 0:
            i_coordinate = i
            break
    return i_coordinate, j_coordinate


def get_edit_distance_multiple(row1, row2, alignment, seq, mat):
    for i in range(seq):
        i = i + 1
        row2[0] = i
        for j in range(len(alignment)):
            j = j + 1
            # implements the recursive definition, view sequence analysis script
            row2[j] = max(row1[j - 1] + get_score_alignment(alignment[j - 1], seq[i - 1], mat), row1[j] + 1, row2[j - 1] + 1)
        row1 = row2
        return row1[len(row1)-1]


def get_hirschberg(row1, row2, alignment, seq, mat):
    limit = int((len(seq)+1)/2) + ((len(seq) + 1) % 2)
    for i in range(limit-1):
        i = i + 1
        row2[0] = i
        for j in range(len(alignment)):
            j = j + 1
            # implements the recursive definition, view sequence analysis script
            row2[j] = max(row1[j - 1] + get_score_alignment(alignment[j - 1], seq[i - 1], mat), row1[j] + 1, row2[j - 1] + 1)
        row1 = row2
    return row1


if __name__ == '__main__':
    main()
