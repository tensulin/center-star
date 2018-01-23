import argparse
import os.path


class DNASequence(object):
    def __init__(self, title, sequence):
        self.title = title
        self.sequence = sequence


UNIT_COST_MATRIX = [[0, 1, 1, 1], [1, 0, 1, 1], [1, 1, 0, 1], [1, 1, 1, 0]]
BASE_LIST = ["A", "C", "G", "T"]
# line length for printing in terminal and file
LINE_LENGTH = 50


def main():
    parser = argparse.ArgumentParser(description='Calculate multiple alignment using the center star approximation')
    parser.add_argument('input', action="store", metavar='file_name', type=str,
                        help='Input fasta file. Sequences should only contain DNA bases (A, C, G, T).')
    parser.add_argument('--dist', action="store", dest="mat_file", metavar='file_name',  default="unit_costs",
                        help='File containing the cost matrix. If left out unit costs will be used. '
                             'For information about the format please read the man.')
    parser.add_argument('--indel', action="store", dest="indel_cost", type=int, metavar='indel_cost', default=1,
                        help='Int to define the indel costs. If left out unit costs will be used.')
    parser.add_argument('--output', action="store", dest="out_file_name", metavar='file_name', default="print_on_terminal",
                        help='Filename where the output will be written. If left out it will print in terminal')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    file_name = args.input
    if not os.path.isfile(file_name):
        print("Error: Input file does not exist")
        return
    dna_sequences_list = fasta_parser(file_name)
    if not dna_sequences_list:
        print("Error is file is no valid fasta file.")

    if args.mat_file != "unit_costs":
        if not os.path.isfile(file_name):
            print("Error: Matrix file does not exist")
            return
        else:
            try:
                mat_file = open(args.mat_file, "r")
                line = mat_file.readline().strip()
                mat_list = line[1:len(line)-1].split(";")
                mat = [parse_list(i) for i in mat_list]
            except:
                print("Error: Format of matrix file is wrong")
    else:
        mat = UNIT_COST_MATRIX
    center_star = find_center_star(dna_sequences_list, mat, args.indel_cost)
    dna_sequences_list.remove(center_star)
    alignment = calculate_alignment(center_star, dna_sequences_list, mat, args.indel_cost)
    if args.out_file_name == "print_on_terminal":
        print("Alignment Order:")
        print(center_star.title)
        for dna_seq in dna_sequences_list:
            print(dna_seq.title)
        print()
        for i in range(int(len(dna_sequences_list[0].sequence) / LINE_LENGTH) + 1):
            for seq in alignment:
                pos = LINE_LENGTH * i
                line  = ""
                for base in seq[pos:pos+LINE_LENGTH]:
                    line += base
                print(line)
            print()
    else:
        out_file = open(args.out_file_name, "w")
        out_file.write("Alignment Order:\n")
        out_file.write(center_star.title)
        for dna_seq in dna_sequences_list:
            out_file.write(dna_seq.title)
        for i in range(int(len(dna_sequences_list[0].sequence) / LINE_LENGTH) + 1):
            for seq in alignment:
                pos = LINE_LENGTH * i
                line = ""
                for base in seq[pos:pos+LINE_LENGTH]:
                    line += base
                out_file.write(line)
            out_file.write("\n")


def parse_list(string):
    out = []
    for c in string:
        try:
            out.append(int(c))
        except ValueError:
            pass
    return out


def find_center_star(dna_sequences_list, mat, indel_cost=1):
    min_cost = -1
    for i in range(len(dna_sequences_list)):
        cost = 0
        dna_sequence = dna_sequences_list[i]
        for j in range(len(dna_sequences_list)):
            if i != j:
                cost += edit_distance(dna_sequence.sequence, dna_sequences_list[j].sequence, mat, indel_cost)
        if min_cost == -1:
            min_cost = cost
            center_star = dna_sequence
        elif cost < min_cost:
            min_cost = cost
            center_star = dna_sequence
    return center_star


def fasta_parser(file_name):
    sep = ""
    dna_sequences_list = []
    dna_seq = DNASequence("", "")
    with open(file_name) as f:
        for line in f:
            if (not line.startswith(">")) and (line != "") and (line[:1] not in BASE_LIST):
                return []
            if line.startswith(">"):
                if dna_seq.sequence != "":
                    new_dna_seq = DNASequence(dna_seq.title, dna_seq.sequence)
                    dna_sequences_list.append(new_dna_seq)
                dna_seq.title = (sep + line.strip() + ' ').lstrip(">")
                dna_seq.sequence = ""
            else:
                dna_seq.sequence += line.strip()
    dna_sequences_list.append(dna_seq)
    return dna_sequences_list


# calculates the edit distance for two DNA sequences and a
# score matrix, if no cost matrix is given, unit costs will be used
def edit_distance(seq1, seq2, mat, indel_cost=1):
    # initialize the first row of matrix
    row1 = [i for i in range(len(seq1) + 1)]
    # only two rows are saved of the matrix to save computational space
    for i in range(1, len(seq2) + 1):
        # initialize fist entry for recursive calculation
        row2 = [0 for i in range(len(seq1) + 1)]
        row2[0] = i
        for j in range(1, len(seq1) + 1):
            # implements the recursive definition, view sequence analysis script
            row2[j] = min(row1[j - 1] + get_cost(seq1[j - 1], seq2[i - 1], mat), row1[j] + indel_cost, row2[j - 1] + indel_cost)
        row1 = row2
    return row1[len(row1) - 1]


# calculates score for 2 bases and a score matrix
def get_cost(base1, base2, mat):
    index1 = BASE_LIST.index(base1)
    index2 = BASE_LIST.index(base2)
    return mat[index1][index2]


def calculate_alignment(center_star, sequences, mat, indel_cost):
    alignment = []
    for seq in sequences:
        next_pair = pair_next_sequence(center_star.sequence, seq.sequence, mat, indel_cost)
        alignment = append_next_sequence(alignment, next_pair)
    return alignment


def append_next_sequence(alignment, next_pair):
    if not alignment:
        return make_alignment(next_pair)
    else:
        new_pair = make_alignment(next_pair)
        index = 0
        old_center = alignment[0]
        new_center = new_pair[0]
        new_sequence = new_pair[1]
        while True:
            if len(old_center) <= index:
                for i in range((len(new_center) - index) + 1):
                    new_sequence.append("-")
                break
            if len(new_center) <= index:
                for i in range((len(old_center) - index)+1):
                    append_new_gap(alignment)
                break
            if old_center == new_center:
                break
            if old_center[index] == new_center[index]:
                pass
            elif old_center[index] == '-':
                new_sequence.insert(index, '-')
                new_center.insert(index, '-')
            elif new_center[index] == '-':
                insert_new_gap(alignment, index)
            index += 1
    alignment.append(new_sequence)
    return alignment


def append_new_gap(alignment):
    for seq in alignment:
        seq.append("-")


def insert_new_gap(alignment, index):
    for seq in alignment:
        seq.insert(index, "-")


def make_alignment(paired_sequence):
    fst_row = []
    snd_row = []
    for i in paired_sequence:
        fst_row.append(i[0])
        snd_row.append(i[1])
    alignment = [fst_row, snd_row]
    return alignment


def pair_next_sequence(center_star, seq, mat, indel_cost):
    new_alignment = []
    index = 0
    new_alignment.append((center_star, seq))
    while 1:
        try:
            center_star_snippet = new_alignment[index][0]
        except IndexError:
            break
        seq_snippet = new_alignment[index][1]
        insert = []
        if len(seq_snippet) == 1 and len(center_star_snippet) == 1:
            del new_alignment[index]
            insert.append(center_star_snippet)
            insert.append(seq_snippet)
            new_alignment.insert(index, insert)
            index += 1
        elif len(seq_snippet) == 0 and len(center_star_snippet) == 1:
            del new_alignment[index]
            insert.append(center_star_snippet)
            insert.append('-')
            new_alignment.insert(index, insert)
            index += 1
        elif len(seq_snippet) == 1 and len(center_star_snippet) == 0:
            del new_alignment[index]
            insert.append('-')
            insert.append(seq_snippet)
            new_alignment.insert(index, insert)
            index += 1
        elif len(seq_snippet) == 0 and len(center_star_snippet) == 0:
            break
        else:
            j_coordinate = (int((len(seq_snippet)) / 2) + ((len(seq_snippet)) % 2))
            hirschberg_row_one = calculate_hirschberg_row(center_star_snippet, seq_snippet, mat, indel_cost, j_coordinate)
            edit_dist = edit_distance(center_star_snippet, seq_snippet, mat)
            center_star_snippet_rev = center_star_snippet[::-1]
            seq_snippet_rev = seq_snippet[::-1]
            hirschberg_row_two = calculate_hirschberg_row(center_star_snippet_rev, seq_snippet_rev, mat, indel_cost, j_coordinate)
            i_coordinate = calculate_i_coordinate(hirschberg_row_one, hirschberg_row_two, edit_dist)
            del new_alignment[index]
            new_alignment.insert(index, (center_star_snippet[:i_coordinate], seq_snippet[:j_coordinate]))
            new_alignment.insert(index + 1, (center_star_snippet[i_coordinate:], seq_snippet[j_coordinate:]))
    return new_alignment


def calculate_i_coordinate(hirschberg_row_one, hirschberg_row_two, edit_dist):
    for i in range(len(hirschberg_row_one)):
        if (hirschberg_row_one[i] + hirschberg_row_two[i]) - edit_dist == 0:
            return i
    return 0


def calculate_hirschberg_row(seq1, seq2, mat, indel_cost, limit):
    row1 = [i for i in range(len(seq1) + 1)]
    for i in range(1, limit+1):
        row2 = [0 for i in range(len(seq1) + 1)]
        row2[0] = i
        for j in range(1, len(seq1) + 1):
            # implements the recursive definition, view sequence analysis script
            row2[j] = min(row1[j - 1] + get_cost(seq1[j - 1], seq2[i - 1], mat),
                          row1[j] + indel_cost, row2[j - 1] + indel_cost)
        row1 = row2
    return row1


if __name__ == '__main__':
    main()
