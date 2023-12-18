import os
from collections import defaultdict
import csv


def read_pathway(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line_split = line.split('\t')
        line_split[2] = line_split[2].replace('\n', '')
        if organism == "Mouse":
            organism = 'Mouse - orthology from human'
            line_split[0] = line_split[0].capitalize()
            line_split[2] = line_split[2].capitalize()
        interactions.append((line_split[0], line_split[2], 'PathwayCommons', line_split[1],
                             '-', '-', '-', '-', '-', '-', '-', '-', '-'))
    return interactions


def read_encode(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line_split = line.split('\t')
        line_0_split = line_split[0].split()
        if organism == "Mouse":
            if line_0_split[-1] == "mm9":
                for i in range(2, len(line_split)):
                    if line_split[i] not in ["", "\n"]:
                        interactions.append((line_0_split[0].capitalize(), line_split[i].capitalize(),
                                            'Enrichr: Encode', '-'
                                            'Mouse', '-', '-', '-', '-', '-', '-', " ".join(line_0_split[1:-1])))
        if organism == "Human":
            if line_0_split[-1] == "hg19":
                for i in range(2, len(line_split)):
                    if line_split[i] not in ["", "\n"]:
                        interactions.append((line_0_split[0], line_split[i], 'Enrichr: Encode',
                                            'Human', '-', '-', '-', '-', '-', '-', '-', " ".join(line_0_split[1:-1])))
    return interactions


def read_chea(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line_split = line.split('\t')
        line_0_split = line_split[0].split()
        if "ChIP-ChIP" in line_0_split:
            chip_mode = "ChIP-ChIP"
        elif "ChIP-Seq" in line_0_split or "Chip-Seq" in line_0_split:
            chip_mode = "ChIP-Seq"
        else:
            chip_mode = "-"
        conditions = " ".join(line_0_split[2:]).replace("Mouse", "").replace("Human", "").replace("ChIP-Seq", "")\
            .replace("ChIP-ChIP", "").replace("Chip-seq", "")
        conditions_final = ' '.join(conditions.split())
        if organism == "Mouse":
            if "Mouse" in line_0_split:
                for i in range(2, len(line_split)):
                    if line_split[i] not in ["", "\n"]:
                        interactions.append((line_0_split[0].capitalize(), line_split[i].capitalize(), 'Enrichr: Chea',
                                            'Mouse', '-', '-', '-', '-', '-', '-', '-', conditions_final, chip_mode))
        if organism == "Human":
            if "Human" in line_0_split:
                for i in range(2, len(line_split)):
                    if line_split[i] not in ["", "\n"]:
                        interactions.append((line_0_split[0], line_split[i], 'Enrichr: Chea',
                                            'Human', '-',  '-', '-', '-', '-', '-', '-', conditions_final, chip_mode))
    return interactions


def read_grndb(file_path, organism, data_type):
    interactions = []
    file = open(file_path, 'r')
    if organism == "Human":
        conditions = file_path.split('/')[1].replace('-regulons.txt', '')
    else:
        conditions = file_path.split('/')[1].replace('whole_', '').replace('-regulons.txt', '')
    for line in file:
        line = line.split('\t')
        interactions.append((line[0], line[1], 'Grndb', organism, '-',
                             data_type, '-', conditions + ": " + line[2],
                             conditions + ": " + line[3], conditions + ": " + line[4], conditions + ": " +
                             line[5].replace("\n", ""), conditions, '-'))
    return interactions


def read_trrust(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split('\t')
        interactions.append((line[0], line[1], 'Trrust', organism, '-', '-', line[2], '-', '-', '-', '-', '-', '-'))
    return interactions


def read_regnet(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split('\t')
        if organism == "Mouse":
            line[0] = line[0].capitalize()
            line[2] = line[2].capitalize()
        if organism == "Human":
            line[0] = line[0].upper()
            line[2] = line[2].upper()
        interactions.append((line[0], line[2], 'RegNetwork', organism, '-', '-', '-', '-', '-', '-', '-', '-', '-'))
    return interactions


def read_all_in_folder_grndb(folder, organism, data_type):
    interactions = []
    files = os.listdir(folder)
    for file in files:
        interactions.extend(read_grndb(folder + '/' + file, organism, data_type))
    return interactions


def merge_lists(lists):
    merged_list = []
    for small_list in lists:
        merged_list.extend(small_list)
    return merged_list


def merge_all_columns_lists_csv(list_of_lists, filename):
    final_file = open(filename, 'w', newline='')
    csv_writer = csv.writer(final_file, delimiter=',')

    merged_list = [item for small_list in list_of_lists for item in small_list]
    all_columns_dict = defaultdict(lambda: [[] for _ in range(12)])
    merged_list.sort(key=lambda x: (x[0], x[1]))

    for row in merged_list:
        key = tuple(row[:2])
        values = row[2:]

        for i, val in enumerate(values):
            if val != '-':
                all_columns_dict[key][i].append(val)

    csv_writer.writerow(["Regulator", "Regulated gene", "Sources", "Organism", "Type of interaction"
                         "Data type", "Mode of regulation", "Best motif",
                         "NES", "Genie3Weight", "Confidence", "Location", "ChIP mode", "Orthology"])

    for key, lists_values in all_columns_dict.items():
        values_with_dashes = [','.join(map(str, set(l))) if l else '-' for l in lists_values]
        csv_writer.writerow([key[0], key[1]] + values_with_dashes)

    final_file.close()
    return final_file


if __name__ == '__main__':
    encode_mouse = read_encode('ENCODE_TF_ChIP-seq_2015.txt', "Mouse")
    chea_mouse = read_chea('ChEA_2022.txt', "Mouse")
    pathway_mouse = read_pathway('PathwayCommons12.All.hgnc.sif', 'Mouse')
    trrust_mouse = read_trrust('trrust_rawdata.mouse.tsv', 'Mouse')
    regnet_mouse = read_regnet('mouse.source', 'Mouse')
    grndb_mouse = read_all_in_folder_grndb('grndb_mouse', 'Mouse', 'Single cell')
    encode_human = read_encode('ENCODE_TF_ChIP-seq_2015.txt', "Human")
    chea_human = read_chea('ChEA_2022.txt', "Human")
    pathway_human = read_pathway('PathwayCommons12.All.hgnc.sif', 'Human')
    print(pathway_human)
    trrust_human = read_trrust('trrust_rawdata.human.tsv', 'Human')
    regnet_human = read_regnet('human.source', 'Human')
    grndb_human_singlecell = read_all_in_folder_grndb('grndb_human_single_cell', 'Human', 'Single cell')
    grndb_human_bulk = read_all_in_folder_grndb('grndb_human_bulk', 'Human', 'Bulk')
    interactions_encode_final = merge_all_columns_lists_csv([encode_mouse, chea_mouse, pathway_mouse,
                                                             trrust_mouse, regnet_human],'interactions_mouse_try.csv')
