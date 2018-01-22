import argparse
import os.path

class DNASequence(object):
    def __init__(self, title, sequence):
        self.title = title
        self.sequence = sequence



UNIT_COST_MATRIX = [[0, 1, 1, 1, 1], [1, 0, 1, 1, 1], [1, 1, 0, 1, 1], [1, 1, 1, 0, 1], [1, 1, 1, 1]]
BASE_LIST = ["A", "C", "G", "T", "-"]


def main():
    parser = argparse.ArgumentParser(description='Calculate multiple alignment using the center star approximation')
    parser.add_argument('--input', metavar='file_name', type=str,
                        help='Input fasta file. Sequences should only contain DNA bases (A, C, G, T).')
    parser.add_argument('--dist', metavar='file_name',  default="unit costs",
                        help='File containing the cost matrix. If left out unit costs will be used. '
                             'For information about the format please read the man.')
    parser.add_argument('--indel', metavar='indel_cost', default="unit_costs",
                        help='Int to define the indel costs. If left out unit costs will be used.')
    parser.add_argument('--output', metavar='file_name', default="print_on_terminal",
                        help='Filename where the output will be written. If left out it will print in terminal')
    args = parser.parse_args()

    os.path.isfile(name)


def fasta_parser(filename):
    sep = ""
    dna_sequences_list = []
    dna_seq = DNASequence("", "")
    with open(filename) as f:
        for line in f:
            if line.startswith(">"):
                if dna_seq.sequence != "":
                    dna_sequences_list.append(dna_seq)
                dna_seq.title = (sep + line.strip() + ' ')
                dna_seq.sequence = ""
            else:
                dna_seq.sequence += line.strip()
    return dna_sequences_list


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
    multiple_alignment = list(center_star)
    for i in range(len(sequences)):
        multiple_alignment = append_next_sequence(multiple_alignment, sequences[i], mat)
    return multiple_alignment


def append_next_sequence(multiple_alignment, seq, mat):
    new_alignment = [(multiple_alignment, seq)]
    index = 0
    number_rows = len(multiple_alignment[0])
    while index < max(len(multiple_alignment), len(seq)):
        if isinstance(new_alignment[index], tuple):
            alignment = new_alignment[index][0]
            sequence = new_alignment[index][1]
            row1 = [i for i in range(len(alignment) + 1)]
            row2 = [len(alignment) + 1]
            hirschberg_row_one = get_hirschberg(row1, row2, alignment, sequence, mat)
            row1 = [i for i in range(len(alignment) + 1)]
            edit_dist = get_edit_distance_multiple(row1, row2, alignment, sequence, mat)
            row1 = [i for i in range(len(alignment) + 1)]
            row1.reverse()
            alignment.reverse()
            sequence.reverse()
            hirschberg_row_two = get_hirschberg(row1, row2, alignment, sequence, mat)
            alignment.reverse()
            sequence.reverse()
            i_coordinate = 0
            j_coordinate = (int((len(sequence) + 1) / 2) + ((len(sequence) + 1) % 2)) - 1
            for i in range(len(hirschberg_row_one)):
                if (hirschberg_row_one[i] + hirschberg_row_two[i]) - edit_dist == 0:
                    i_coordinate = i - 1
                    break
            if len(sequence) == 1 and len(alignment) == 1:
                alignment[0].append(sequence[0])
                new_alignment.insert(index, alignment[0])
                index += 1
            if len(sequence) == 0 and len(alignment) == 1:
                alignment[0].append('-')
                new_alignment.insert(index, alignment[0])
                index += 1
            if len(sequence) == 1 and len(alignment) == 0:
                new_column = list(number_rows * '-')
                new_column.append(sequence[0])
                new_alignment.insert(index, new_column)
                index += 1
            del new_alignment[index]
            new_alignment.insert(index, ([i[:i_coordinate] for i in alignment], sequence[:j_coordinate]))
            new_alignment.insert(index + 1, ([i[i_coordinate:] for i in alignment], sequence[j_coordinate:]))
    return new_alignment

def get_first_coordinate(coordinates_list):
    coordinate = 0
    for i in coordinates_list:
        if i != 0:
            coordinate = i
        else:
            return coordinate


def get_second_coordinate(coordinates_list):
    last_i = -1
    for i in coordinates_list:
        if last_i == 0 and i != 0:
            return i
        last_i = i


def get_next_iteration_hirschberg(alignment, sequence, mat):
    row1 = [i for i in range(len(alignment) + 1)]
    row2 = [len(alignment) + 1]
    hirschberg_row_one = get_hirschberg(row1, row2, alignment, sequence, mat)
    row1 = [i for i in range(len(alignment) + 1)]
    edit_dist = get_edit_distance_multiple(row1, row2, alignment, sequence, mat)
    row1 = [i for i in range(len(alignment) + 1)]
    row1.reverse()
    alignment.reverse()
    sequence.reverse()
    hirschberg_row_two = get_hirschberg(row1, row2, alignment, sequence, mat)
    alignment.reverse()
    sequence.reverse()
    i_coordinate = 0
    j_coordinate = int((len(sequence) + 1) / 2) + ((len(sequence) + 1) % 2)
    for i in range(len(hirschberg_row_one)):
        if (hirschberg_row_one[i] + hirschberg_row_two[i]) - edit_dist == 0:
            i_coordinate = i
            break
    if len(sequence) == 1 and len(alignment) == 1:
        return alignment[0], sequence[0]
    if len(sequence) == 0 and len(alignment) == 1:
        return alignment[0], '-'
    if len(sequence) == 1 and len(alignment) == 0:
        return sequence[0], '-'
    return get_next_iteration_hirschberg([i[:i_coordinate] for i in alignment], sequence[:j_coordinate]), get_next_iteration_hirschberg([i[i_coordinate:] for i in alignment], sequence[j_coordinate:])


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
