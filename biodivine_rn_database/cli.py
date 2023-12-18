from collections import defaultdict
import csv
from typing import List, Dict, Union, Optional, Tuple
import pickle

# method to create dictionary for finding synonyms


def read_synonyms(file: str) -> Dict[str, List[str]]:
    dictionary = {}
    with open(file, 'r') as file_human:
        for i, line in enumerate(file_human):
            if i == 0:
                continue
            columns = line.split('\t')
            columns_2 = columns[2]
            columns_4 = set(columns[4].split('|'))
            if columns_4 != {'-'}:
                dictionary[columns_2] = list(columns_4)
            else:
                dictionary[columns_2] = []
    return dictionary


# method to find value in dictionary


def find_value_in_dict(dictionary: Dict[str, List[str]], value: str) -> Optional[Tuple[str, List[str]]]:
    for key, value_list in dictionary.items():
        if value in value_list:
            return key, value_list
    return None


# method to find most known game for a gene


def find_most_known_name(synonyms: Dict[str, List[str]], name: str, organism: str) -> str:
    if name.lower() in [k.lower() for k in synonyms.keys()]:
        return name
    for key, value_list in synonyms.items():
        if name.lower() in [v.lower() for v in value_list]:
            if organism == "Human":
                return key.upper()
            if organism == "Mouse":
                return key.capitalize()
    if organism == "Human":
        return name.upper()
    return name.capitalize()


# methods to find successors and predecessors


def find_successors_interactions(interactions_file: str) -> Dict[str, List[Tuple[str, str]]]:
    successors = defaultdict(list)
    with open(interactions_file, 'r') as interactions:
        csv_reader = csv.reader(interactions)
        for interaction in csv_reader:
            successors[interaction[0]].append((interaction[1], interaction[2]))
    return successors


def find_predecessors_interactions(interactions_file: str) -> Dict[str, List[Tuple[str, str]]]:
    predecessors = defaultdict(list)
    with open(interactions_file, 'r') as interactions:
        csv_reader = csv.reader(interactions)
        for interaction in csv_reader:
            predecessors[interaction[1]].append((interaction[0], interaction[2]))
    return predecessors


# methods to find all target genes for TF


def find_tgs_for_tf(successors: Dict[str, List[Tuple[str, str]]],
                    synonyms: Dict[str, List[str]], organism: str) -> None:
    while True:
        gene = input("Input a transcription factor gene or type exit to exit this mode: ")
        if gene.lower() == "exit":
            break
        known_name = find_most_known_name(synonyms, gene, organism)
        if known_name is not None:
            gene = known_name
        if any(gene.lower() == key.lower() for key in successors):
            print_info(gene, synonyms, successors, True)
        else:
            print("This gene is not correct or is not among successors.")


# method to find all transcription factors for target gene

def find_tfs_for_tg(predecessors: Dict[str, List[Tuple[str, str]]],
                    synonyms: Dict[str, List[str]], organism: str) -> None:
    while True:
        gene = input("Input a target gene or type exit to exit this mode: ")
        if gene.lower() == "exit":
            break
        known_name = find_most_known_name(synonyms, gene, organism)
        if known_name is not None:
            gene = known_name
        if any(gene.lower() == key.lower() for key in predecessors):
            print_info(gene, synonyms, predecessors, False)
        else:
            print("This gene is not correct or is not among predecessors.")


def print_info(gene: str, synonyms: Dict[str, List[str]],
               successors_or_predecessors: Dict[str, List[Tuple[str, str]]], succ_or_pred_bool: bool) -> None:
    if succ_or_pred_bool:
        print("gene " + gene + " regulates these genes")
    else:
        print("gene " + gene + " is regulated by these genes")
    for key, values in successors_or_predecessors.items():
        if key.lower() == gene.lower():
            for tf_or_tg in values:
                synonyms_for_tf_or_tg = find_value_in_dict(synonyms, tf_or_tg[0])
                if tf_or_tg[0] not in synonyms.keys() and synonyms_for_tf_or_tg is not None:
                    print(tf_or_tg, end='')
                    print(", Commonly known as: " + synonyms_for_tf_or_tg[0], end='')
                    print(", Other synonyms: ", end='')
                    for synonym in synonyms_for_tf_or_tg[1]:
                        print(synonym, end=' ')
                else:
                    print(tf_or_tg, end='')
                print()


def find_info_about_tf_and_tg(table: List[List[str]], synonyms: Dict[str, List[str]], organism: str) -> None:
    while True:
        found_interaction = False
        transcription_factor = input("Input a transcription factor or exit to exit this mode: ")
        if transcription_factor == "exit":
            break
        target_gene = input("Input a target gene or exit to exit this mode: ")
        if target_gene == "exit":
            break

        transcription_factor = find_most_known_name(synonyms, transcription_factor, organism)
        target_gene = find_most_known_name(synonyms, target_gene, organism)

        if transcription_factor is not None and target_gene is not None:
            print(transcription_factor, target_gene)
            searched_interaction = binary_search_table(transcription_factor, target_gene, table, organism)
            if searched_interaction is not None:
                print(','.join(table[0]))
                print(searched_interaction)
                found_interaction = True

        if not found_interaction:
            print("Interaction between these two genes was not found.")


def find_info_about_group(table: List[List[str]], synonyms: Dict[str, List[str]], organism: str) -> None:
    while True:
        found_interaction = False
        genes = input("Input genes separated by spaces or exit to exit this mode: ")
        if genes == "exit":
            break
        genes_split = genes.split(" ")
        normalized_genes = []
        for gene in genes_split:
            gene = find_most_known_name(synonyms, gene, organism)
            normalized_genes.append(gene)
        if len(normalized_genes) < 2:
            print("There must be at least two genes.")
            continue

        for tf in normalized_genes:
            for tg in normalized_genes:
                if tf is not None and tg is not None:
                    searched_interaction = binary_search_table(tf, tg, table, organism)
                    if searched_interaction is not None:
                        if not found_interaction:
                            print(','.join(table[0]))
                        print(searched_interaction)
                        found_interaction = True

        if not found_interaction:
            print("Some of the genes are not correct or interaction between these genes was not found.")


def get_header_from_csv(csv_file: str) -> Optional[List[str]]:
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        header = next(csv_reader, None)
    return header


def find_main_names_for_synonyms(synonyms: Dict[str, List[str]]) -> None:
    while True:
        gene = input("Input symbol for a gene or exit to exit this mode: ")
        if gene.lower() == "exit":
            break
        possible_genes = []
        for key, values in synonyms.items():
            if gene.lower() == key.lower():
                possible_genes.append(key)
            if gene.lower() in [v.lower() for v in values]:
                possible_genes.append(key)
        print("These are main names of genes that can use this symbol " + str(possible_genes))


def print_info_about_gene(gene_table: str) -> None:
    gene_table_opened = open(gene_table, 'r')
    while True:
        found_interaction = False
        gene = input("Input a gene to find info about or exit to exit this mode: ")
        if gene.lower() == 'exit':
            break
        for line in gene_table_opened:
            line_split = line.split('\t')
            if line_split[2].lower() == gene.lower():
                found_interaction = True
                print("TaxID: " + line_split[0] + '\n' +
                      "GeneID: " + line_split[1] + '\n' +
                      "Symbol: " + line_split[2] + '\n' +
                      "Synonyms: " + line_split[4] + '\n' +
                      "Description: " + line_split[8] + '\n' +
                      "Type of gene: " + line_split[9] + '\n' +
                      "Symbol from nomenclature authority: " + line_split[10] + '\n')
                break
        if not found_interaction:
            print("Selected gene does not exist or is not among main names.")


def start_application(interactions_human: List[List[str]], successors_human: Dict[str, List[Tuple[str, str]]],
                      predecessors_human: Dict[str, List[Tuple[str, str]]], synonyms_human: Dict[str, List[str]],
                      gene_table_human: str, organism_human: str,
                      interactions_mouse: List[List[str]], successors_mouse: Dict[str, List[Tuple[str, str]]],
                      predecessors_mouse: Dict[str, List[Tuple[str, str]]], synonyms_mouse: Dict[str, List[str]],
                      gene_table_mouse: str, organism_mouse: str) -> None:
    while True:
        organism = input("Type \033[91mhuman\033[0m to find regulations about Homo Sapiens \n"
                         "\033[91mmouse\033[0m to find regulation about Mus musculus \n"
                         "\033[91mexit\033[0m to end the application: ")
        if organism.lower() == "human":
            start_application_with_organism(interactions_human, successors_human, predecessors_human,
                                            synonyms_human, gene_table_human, organism_human)
        elif organism.lower() == "mouse":
            start_application_with_organism(interactions_mouse, successors_mouse, predecessors_mouse,
                                            synonyms_mouse, gene_table_mouse, organism_mouse)
        elif organism.lower() == "exit":
            break
        else:
            print("Incorrect mode")


def start_application_with_organism(interactions: List[List[str]], successors: Dict[str, List[Tuple[str, str]]],
                                    predecessors: Dict[str, List[Tuple[str, str]]], synonyms: Dict[str, List[str]],
                                    gene_table: str, organism: str) -> None:
    while True:
        mode = input("Type \033[91mTF\033[0m to choose TF and find TGs for it\n\033[91mTG\033[0m to choose TG and find"
                     " all TFs that regulate it \n\033[91minfo\033[0m to find info"
                     " about TF and TG\n\033[91mgroup\033[0m to find interactions in "
                     "group of genes\n\033[91mgene\033[0m to print info about gene \n\033[91mnames\033[0m "
                     "to find all possible genes for given symbol"
                     "\n\033[91mexit\033[0m to choose organism: ")
        if mode.upper() == "TF":
            find_tgs_for_tf(successors, synonyms, organism)
        elif mode.upper() == "TG":
            find_tfs_for_tg(predecessors, synonyms, organism)
        elif mode.lower() == "info":
            find_info_about_tf_and_tg(interactions, synonyms, organism)
        elif mode.lower() == "group":
            find_info_about_group(interactions, synonyms, organism)
        elif mode.lower() == "gene":
            print_info_about_gene(gene_table)
        elif mode.lower() == "names":
            find_main_names_for_synonyms(synonyms)
        elif mode.lower() == "exit":
            break
        else:
            print("Incorrect mode.")


# method to merge interactions


def merge_lists(lists: List[List[str]]) -> List[str]:
    merged_list = []
    for small_list in lists:
        for item in small_list:
            merged_list.append(item)
    return merged_list


def make_table_from_file(csv_file: str) -> List[List[str]]:
    table = []
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file)
        for line in csv_reader:
            table.append(line)
    return table


# binary searcher for all lines with given tf

def binary_search_table(tf: str, tg: str, table: List[List[str]], organism: str) -> Optional[List[str]]:
    if organism == "Human":
        tf = tf.upper()
        tg = tg.upper()
    if organism == "Mouse":
        tf = tf.capitalize()
        tg = tg.capitalize()
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


def pickle_all_files():
    data_sources = [
        ('predecessors_human.pkl', find_predecessors_interactions('interactions_human.csv')),
        ('predecessors_mouse.pkl', find_predecessors_interactions('interactions_mouse.csv')),
        ('successors_human.pkl', find_successors_interactions('interactions_human.csv')),
        ('successors_mouse.pkl', find_successors_interactions('interactions_mouse.csv')),
        ('synonyms_human.pkl', read_synonyms('Homo_sapiens.gene_info')),
        ('synonyms_mouse.pkl', read_synonyms('Mus_musculus.gene_info')),
    ]

    for file_name, data in data_sources:
        with open(file_name, 'wb') as file:
            pickle.dump(data, file)


def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data


if __name__ == '__main__':
    synonyms_1 = load_pickle('synonyms_human.pkl')
    synonyms_2 = load_pickle('synonyms_mouse.pkl')
    successors_1 = load_pickle('successors_human.pkl')
    predecessors_1 = load_pickle('predecessors_human.pkl')
    successors_2 = load_pickle('successors_mouse.pkl')
    predecessors_2 = load_pickle('predecessors_mouse.pkl')
    all_genes_human = 'Homo_sapiens.gene_info'
    all_genes_mouse = 'Mus_musculus.gene_info'
    interactions_human = make_table_from_file('interactions_human.csv')
    interactions_mouse = make_table_from_file('interactions_mouse.csv')
    start_application(interactions_human, successors_1, predecessors_1, synonyms_1, all_genes_human, "Human",
                      interactions_mouse, successors_2, predecessors_2, synonyms_2, all_genes_mouse, "Mouse")
