import os
from collections import defaultdict
import csv


# methods to go through normalization database

def count_for_synonyms(synonyms):
    synonyms_count = {}
    for key, values in synonyms.items():
        if key not in synonyms_count:
            synonyms_count[key] = 1
        else:
            synonyms_count[key] += 1
        for value in values:
            if value not in synonyms_count:
                synonyms_count[value] = 1
            else:
                synonyms_count[value] += 1
    sorted_dict = dict(sorted(synonyms_count.items(), key=lambda item: item[1]))
    return sorted_dict


# def merge_normalization_table(file_human, file_mouse, final_file):
#     file_human = open(file_human, 'r')
#     file_mouse = open(file_mouse, 'r')
#     final_file = open(final_file, 'w')
#     reader_human = csv.reader(file_human, delimiter='\t')
#     reader_mouse = csv.reader(file_mouse, delimiter='\t')
#     final_file.write('TaxID' + ';' + 'GeneID' + ';' + 'Symbol' + ';' + 'Synonyms' + ';' +
#                      'Description' + ';' + 'Type of gene' + ';'
#                      + 'Symbol from nomenclature authority' + ';' + 'Full name from nomenclature authority'
#                      + ';' + 'Organism' + '\n')
#     first_line = True
#     for line in reader_human:
#         if first_line:
#             first_line = False
#             continue
#         for i in range(12):
#             if i in (0, 1, 2, 4, 8, 9, 10, 11):
#                 final_file.write(line[i] + ';')
#         final_file.write('Homo sapiens' + '\n')
#     first_line = True
#     for line in reader_mouse:
#         if first_line:
#             first_line = False
#             continue
#         for i in range(12):
#             if i in (0, 1, 2, 4, 8, 9, 10, 11):
#                 final_file.write(line[i] + ';')
#         final_file.write('Mus musculus' + '\n')
#     return final_file


def all_main_names(file_path):
    with open(file_path, 'r') as file_path:
        next(file_path)
        original_names = [line.split()[2] for line in file_path]
    return original_names


def read_synonyms(file):
    dictionary = {}
    with open(file, 'r') as file_human:
        for i, line in enumerate(file_human):
            if i == 0:
                continue
            columns = line.split()
            columns_2 = columns[2].capitalize()
            columns_4 = set(columns[4].split('|'))
            if columns_4 != {'-'}:
                dictionary[columns_2] = columns_4
            else:
                dictionary[columns_2] = []
    return dictionary


# methods to find tuples from databases


def read_pathway(file):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        interactions.append((line[0], line[2], 'PathwayCommons', '-', '-', '-', '-', '-', '-', '-'))
    return interactions


def read_enrichr(file):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        first_column = line[0].split('_')
        for i in range(1, len(line)):
            interactions.append((first_column[0], line[i], 'Enrichr', '-', '-', '-', '-', '-', '-', '-'))
    return interactions


def read_grndb(file, organism, data_type):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        if line[4] == 'High' or line[4] == 'Low':
            line.append(line[4])
            line[4] = '-'
        interactions.append((line[0], line[1], 'Grndb', organism, data_type, '-', line[2], line[3], line[4], line[5]))
    return interactions


def read_trrust(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        interactions.append((line[0], line[1], 'Trrust', organism, '-', line[2], '-', '-', '-', '-'))
    return interactions


def read_regnet(file, organism):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        interactions.append((line[0], line[2], 'RegNetwork', organism, '-', '-', '-', '-', '-', '-'))
    return interactions


def read_all_in_folder_grndb(folder, organism, data_type):
    interactions = []
    files = os.listdir(folder)
    for file in files:
        interactions.extend(read_grndb(folder+'/'+file, organism, data_type))
    return interactions

# method to find value in dictionary


def find_value_in_dict(dictionary, value):
    for key, value_list in dictionary.items():
        if value in value_list:
            return key, value_list
    return None


def find_most_known_name(synonyms, name, organism_names):
    main_names = all_main_names(organism_names)
    if name in main_names:
        return name
    for key, value_list in synonyms.items():
        if name in value_list:
            return key
    return name

# methods to find successors and predecessors


def find_successors_interactions(interactions):
    successors = defaultdict(list)
    interactions = open(interactions, 'r')
    for interaction in interactions:
        interaction = interaction.split(';')
        successors[interaction[0]].append((interaction[1], interaction[2]))
    return successors


def find_predecessors_interactions(interactions):
    predecessors = defaultdict(list)
    interactions = open(interactions, 'r')
    for interaction in interactions:
        interaction = interaction.split(';')
        predecessors[interaction[1]].append((interaction[0], interaction[2]))
    return predecessors


# methods to find all target genes for TF

def find_tgs_for_tf(successors, synonyms, organism_names):
    while True:
        gene = input("Input a transcription factor gene or type exit to exit the application: ")
        main_names = all_main_names(organism_names)
        if gene.lower() == "exit":
            break
        if gene.lower() not in [name.lower() for name in main_names]:
            known_name = find_most_known_name(synonyms, gene, organism_names)
            if known_name is not None:
                gene = known_name
        if any(gene.lower() == key.lower() for key in successors):
            print_info(gene, synonyms, main_names, successors, True, organism_names)
        else:
            print("This gene is not correct or is not among successors.")


# method to find all transcription factors for target gene

def find_tfs_for_tg(predecessors, synonyms, organism_names):
    while True:
        gene = input("Input a target gene or type exit to exit the application: ")
        main_names = all_main_names(organism_names)
        if gene.lower() == "exit":
            break
        if gene.lower() not in [name.lower() for name in main_names]:
            known_name = find_most_known_name(synonyms, gene, organism_names)
            if known_name is not None:
                gene = known_name
        if any(gene.lower() == key.lower() for key in predecessors):
            print_info(gene, synonyms, main_names, predecessors, False, organism_names)
        else:
            print("This gene is not correct or is not among predecessors.")


def print_info(gene, synonyms, main_names, successors_or_predecessors, succ_or_pred_bool, organism_names):
    if succ_or_pred_bool:
        print("gene " + gene + " regulates these genes")
    else:
        print("gene " + gene + " is regulated by these genes")
    for key, values in successors_or_predecessors.items():
        if key.lower() == gene.lower():
            for tf_or_tg in values:
                print(tf_or_tg, end='')
                synonyms_for_transcription_factor = find_value_in_dict(synonyms, tf_or_tg)
                if tf_or_tg[0] in main_names:
                    if tf_or_tg[0] in synonyms:
                        print(" Synonyms: ", end='')
                        for synonym in synonyms[tf_or_tg[0]]:
                            print(synonym, end=' ')
                    print()
                elif synonyms_for_transcription_factor is not None:
                    print(" Commonly known as: " + synonyms_for_transcription_factor[0], end='')
                    print(", Other synonyms: ", end='')
                    for synonym in synonyms_for_transcription_factor[1]:
                        print(synonym, end=' ')
                    print()
                else:
                    print()
    print('\n')


def find_info_about_tf_and_tg(interactions, synonyms, organism_names):
    while True:
        found_interaction = False
        transcription_factor = input("Input a transcription factor: ")
        target_gene = input("Input a target gene: ")
        transcription_factor = find_most_known_name(synonyms, transcription_factor, organism_names)
        target_gene = find_most_known_name(synonyms, target_gene, organism_names)
        table = make_table_from_file(interactions)

        for tf_synonym in [transcription_factor, *synonyms[transcription_factor]]:
            for tg_synonym in [target_gene, *synonyms[target_gene]]:
                searched_interaction = binary_search_table(tf_synonym, tg_synonym, table)
                if searched_interaction is not None:
                    print(searched_interaction)
                    found_interaction = True

        if not found_interaction:
            print("Interaction between these two genes was not found.")
        continue_inputting = input("If you want to exit, type E, else press enter: ")
        if continue_inputting == "E":
            break


def find_info_about_group(interactions_file, synonyms, organism_names):
    while True:
        found_interaction = False
        genes = input("Input genes separated by spaces: ")
        genes = genes.split(" ")
        table = make_table_from_file(interactions_file)
        normalized_genes = []
        for gene in genes:
            gene = find_most_known_name(synonyms, gene, organism_names)
            normalized_genes.append(gene)
        if len(normalized_genes) < 2:
            print("There must be at least two genes.")

        for tf in normalized_genes:
            for tg in normalized_genes:
                for tf_synonym in [tf, *synonyms[tf]]:
                    for tg_synonym in [tg, *synonyms[tg]]:
                        searched_interaction = binary_search_table(tf_synonym, tg_synonym, table)
                        if searched_interaction is not None:
                            print(searched_interaction)
                            found_interaction = True

        if not found_interaction:
            print("Interaction between these two genes was not found.")
        continue_inputting = input("If you want to exit, type E, else press enter: ")
        if continue_inputting == "E":
            break


def print_info_about_gene(synonyms, gene_table):
    while True:
        gene = input("Input a gene to find info about: ")
        appearances = 0
        for key, value in synonyms.items():
            if gene.lower() == key.lower() or gene.lower() in [name.lower() for name in value]:
                appearances += 1

        if appearances == 0:
            print("This gene was not found.")
        if appearances == 1:
            print_info_about_gene_final(gene, gene_table)
        if appearances == 2:
            print("More appearances.")
        continue_inputting = input("If you want to exit, type E, else press enter: ")
        if continue_inputting == "E":
            break


def print_info_about_gene_final(gene, gene_table):
    gene_table = open(gene_table, 'r')
    for line in gene_table:
        line = line.split(';')
        if line[2] == gene:
            print("TaxID: " + line[0] + '\n' +
                  "GeneID: " + line[1] + '\n' +
                  "Symbol: " + line[2] + '\n' +
                  "Synonyms: " + line[3] + '\n' +
                  "Description: " + line[4] + '\n' +
                  "Type of gene: " + line[5] + '\n' +
                  "Symbol from nomenclature authority: " + line[6] + '\n' +
                  "Full name nomenclature authority: " + line[7] + '\n' +
                  "Organism: " + line[8])
            break


def start_application(interactions_human, successors_human, predecessors_human, synonyms_human, gene_table_human,
                      interactions_mouse, successors_mouse, predecessors_mouse, synonyms_mouse, gene_table_mouse):
    while True:
        organism = input("Type \033[91mhuman\033[0m to find regulations about Homo Sapiens \n"
                         "\033[91mmouse\033[0m to find regulation about Mus musculus \n"
                         "\033[91mexit\033[0m to end the application: ")
        if organism.lower() == "human":
            start_application_with_organism(interactions_human, successors_human, predecessors_human,
                                            synonyms_human, gene_table_human)
        elif organism.lower() == "mouse":
            start_application_with_organism(interactions_mouse, successors_mouse, predecessors_mouse,
                                            synonyms_mouse, gene_table_mouse)
        elif organism.lower() == "exit":
            break
        else:
            print("Incorrect mode")


def start_application_with_organism(interactions, successors, predecessors, synonyms, gene_table):
    while True:
        mode = input("Type \033[91mTF\033[0m to choose TF and find TGs for it\n\033[91mTG\033[0m to choose TG and find"
                    " all TFs that regulate it \n\033[91minfo\033[0m to find info"
                    " about TF and TG\n\033[91mgroup\033[0m to find interactions in "
                    "group of genes\n\033[91mgene\033[0m to print info about gene"
                    "\n\033[91mexit\033[0m to end the application: ")
        if mode.upper() == "TF":
            find_tgs_for_tf(successors, synonyms, gene_table)
        elif mode.upper() == "TG":
            find_tfs_for_tg(predecessors, synonyms, gene_table)
        elif mode.lower() == "info":
            find_info_about_tf_and_tg(interactions, synonyms, gene_table)
        elif mode.lower() == "group":
            find_info_about_group(interactions, synonyms, gene_table)
        elif mode.lower() == "gene":
            print_info_about_gene(synonyms, gene_table)
        elif mode.lower() == "exit":
            break
        else:
            print("Incorrect mode.")
# method to merge interactions


def merge_lists(lists):
    merged_list = []
    for small_list in lists:
        for item in small_list:
            merged_list.append(item)
    return merged_list


# method to merge interaction from all databases that keeps sources

def merge_triples_lists(list_of_lists, filename):
    final_file = open(filename, 'w')
    lists = merge_lists(list_of_lists)
    triple_dict = defaultdict(list)
    for triple in lists:
        key = (triple[0], triple[1])
        if triple[2] not in triple_dict[key]:
            triple_dict[key].append(triple[2])

    for key, values in triple_dict.items():
        values_str = ','.join(values)
        final_file.write(key[0] + ';' + key[1] + ';' + values_str)
        final_file.write('\n')

    final_file.close()
    return final_file


def merge_ten_columns_lists(list_of_lists, filename):
    final_file = open(filename, 'w')
    lists = merge_lists(list_of_lists)
    ten_columns_dict = defaultdict(lambda: [set() for _ in range(8)])

    for row in lists:
        key = tuple(row[:2])  # Assuming the first two columns are the unique identifiers
        values = row[2:]      # The rest of the columns are treated as values for the tuple

        for i, val in enumerate(values):
            if val != '-':
                ten_columns_dict[key][i].add(val)

    sorted_items = sorted(ten_columns_dict.items(), key=lambda x: (x[0][0], x[0][1]))

    final_file.write("Regulator" + ";" + "Regulated gene" + ";" + "Sources" + ";" + "Organism" + ";"
                     "Data type" + ";" + "Mode of regulation" + ";" + "Best motif" + ";" +
                     "NES" + ";" + "Genie3Weight" + ";" + "Confidence" + "\n")

    for key, sets_list in sorted_items:
        values_with_dashes = [list(s) if s else ['-'] for s in sets_list]
        values_str = ';'.join(','.join(map(str, s)) for s in values_with_dashes)
        final_file.write(key[0] + ';' + key[1] + ';' + values_str)
        final_file.write('\n')

    final_file.close()
    return final_file


def make_table_from_file(table_file):
    table_file = open(table_file, 'r')
    table = []
    for line in table_file:
        line = line.split(';')
        table.append(line)
    return table


# binary searcher for all lines with given tf

def binary_search_table(tf, tg, table):
    low, high = 0, len(table) - 1

    while low <= high:
        mid = (low + high) // 2
        line = table[mid]
        if line[0] == tf and line[1] == tg:
            return table[mid]
        elif line[0] < tf or (line[0] == tf and line[1] < tg):
            low = mid + 1
        else:
            high = mid - 1

    return None


if __name__ == '__main__':
    #all_genes = merge_normalization_table('Homo_sapiens.gene_info', 'Mus_musculus.gene_info', 'Normalization_merged.csv')
    # all_genes = 'Normalization_merged.csv'
    # synonyms_1 = read_synonyms(all_genes)
    synonyms_1 = read_synonyms('Homo_sapiens.gene_info')
    synonyms_2 = read_synonyms('Mus_musculus.gene_info')
    all_genes_human = 'Homo_sapiens.gene_info'
    all_genes_mouse = 'Mus_musculus.gene_info'
    orig_names_human = all_main_names(all_genes_human)
    orig_names_mouse = all_main_names(all_genes_mouse)
    #print(count_for_synonyms(synonyms_1))
    # interactions_human_bulk = read_all_in_folder_grndb('grndb_human_bulk', 'Homo Sapiens', 'Bulk')
    # interactions_human_sc = read_all_in_folder_grndb('grndb_human_single_cell', 'Homo Sapiens', 'Single cell')
    # interactions_mouse = read_all_in_folder_grndb('grndb_mouse', 'Mus musculus', 'Single cell')
    # interactions_mouse_merged = merge_ten_columns_lists([interactions_mouse], 'interactions_mouse.csv')
    # final_interactions = merge_eight_columns_lists([interactions_3], 'interactions4.csv')
    # interactions_1 = read_pathway('homo-sapiens-9606.sif', synonyms_1)
    # interactions_2 = read_grndb('AML_TCGA-regulons.txt', synonyms_1)
    # interactions_3 = read_trrust('trrust_rawdata.human.tsv', synonyms_1)
    # interactions_4 = read_erichr('ARCHS4_Coexpression.gmt', synonyms_1)
    # final_interactions = merge_eight_columns_lists([interactions_1, interactions_2, interactions_3, interactions_4],
    #                                        'interactions3.csv')
    final_interactions = 'interactions3.csv'
    interactions_mouse = 'interactions_mouse.csv'
    successors_1 = find_successors_interactions(final_interactions)
    predecessors_1 = find_predecessors_interactions(final_interactions)
    successors_2 = find_successors_interactions(interactions_mouse)
    predecessors_2 = find_predecessors_interactions(interactions_mouse)
    start_application(final_interactions, successors_1, predecessors_1, synonyms_1, all_genes_human,
                      interactions_mouse, successors_2, predecessors_2, synonyms_2, all_genes_mouse)
