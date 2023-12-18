from collections import defaultdict
import csv
from biodivine_rn_database.cli import find_most_known_name, make_table_from_file, binary_search_table, load_pickle
from typing import List, Tuple, Optional, Dict, Union


def return_tgs_for_tf(gene: str, is_human: bool, normalized: bool) -> List[Tuple[str, str]]:
    """
    Input parameters:
    gene: string, symbol of gene for which user wants to find all target_genes
    is_human: boolean, True if organism is human, False if organism is mouse
    normalized: True if gene symbols should be normalized and False if same as in original databases
    Output:
    Returns list of tuples of target gene and sources of information about this interaction or None if
    there are no successors or the input is incorrect.
    """
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
        successors = load_pickle('successors_human.pkl')
        gene = find_most_known_name(synonyms, gene, "Human")
        organism = "Human"
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
        successors = load_pickle('successors_mouse.pkl')
        gene = find_most_known_name(synonyms, gene, "Mouse")
        organism = "Mouse"
    if not normalized:
        if gene not in successors.keys():
            return []
        else:
            return successors[gene]
    return normalized_finder(gene, synonyms, organism, successors)


def return_tfs_for_tg(gene: str, is_human: bool, normalized: bool) -> List[Tuple[str, str]]:
    """
    Input parameters:
    gene: string, symbol of gene for which user wants to find all transcription factors
    is_human: boolean, True if organism is human, False if organism is mouse
    normalized: True if gene symbols should be normalized and False if same as in original databases
    Output:
    Returns tuples of transcription factor and sources of information about this interaction or None if there
    are no predecessors or the input is incorrect..
    """
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
        predecessors = load_pickle('predecessors_human.pkl')
        gene = find_most_known_name(synonyms, gene, "Human")
        organism = "Human"
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
        predecessors = load_pickle('predecessors_mouse.pkl')
        gene = find_most_known_name(synonyms, gene, "Mouse")
        organism = "Mouse"
    if not normalized:
        if gene not in predecessors.keys():
            return []
        else:
            return predecessors[gene]
    return normalized_finder(gene, synonyms, organism, predecessors)


def normalized_finder(gene: str, synonyms: Dict[str, List[str]], organism: str,
                      successors_or_predecessors: Dict[str, List[Tuple[str, str]]]) -> Optional[List[Tuple[str, str]]]:
    returned_genes = []
    for key, values in successors_or_predecessors.items():
        if key.lower() == gene.lower():
            for tf_or_tg in values:
                returned_genes.append((find_most_known_name(synonyms, tf_or_tg[0], organism), tf_or_tg[1]))
    return returned_genes


def return_info_about_two_genes(transcription_factor: str, target_gene: str, is_human: bool) -> Optional[List[str]]:
    """
    Input parameters:
    transcription_factor: string, transcription factor in the interaction
    target_gene: string, target gene in the interaction
    is_human: boolean, True if organism is human, False if organism is mouse
    Output:
    Returns list of information in format ["Regulator", "Regulated gene", "Sources", "Organism", "Type of interaction"
    "Data type", "Mode of regulation", "Best motif", "NES", "Genie3Weight", "Confidence", "Location",
     "ChIP mode", "Orthology"]
    about interaction between this transcription factor and target gene or None if there
    is no interactions or one of the arguments is incorrect.
    """
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
        organism = "Human"
        table = make_table_from_file('interactions_human.csv')
        tf = find_most_known_name(synonyms, transcription_factor, "Human")
        tg = find_most_known_name(synonyms, target_gene, "Human")
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
        organism = "Mouse"
        table = make_table_from_file('interactions_mouse.csv')
        tf = find_most_known_name(synonyms, transcription_factor, "Mouse")
        tg = find_most_known_name(synonyms, target_gene, "Mouse")
    if tf is not None and tg is not None:
        result = binary_search_table(tf, tg, table, organism)
        if result is not []:
            return result
    return None


def return_info_about_group(list_of_genes: List[str], is_human: bool, detailed_info: bool) -> \
        Optional[Union[List[Tuple[str, str, str]], List[List[str]]]]:
    """
    Input parameters:
    list_of_genes: List[string], list of symbols for genes in interactions
    is_human: boolean, True if organism is human, False if organism is mouse
    detailed_info: boolean, True if user detailed info about interactions should be returned and False if only
    transcription factor, target gene and sources for each interactions should be returned
    Output:
    If detailed_info is True, function returns list of information about interactions between genes in given group.
    If detailed_info is False, function returns list of tuples, which consist of transcription factor,
    in format ["Regulator", "Regulated gene", "Sources", "Organism", "Type of interaction"
    "Data type", "Mode of regulation", "Best motif", "NES", "Genie3Weight", "Confidence", "Location",
     "ChIP mode", "Orthology"]
    target gene among given genes and sources of these interactions.
    Function returns None if there are no interaction, just 0 or 1 genes in list_of_genes or
    incorrect genes in list_of_genes.
    """
    if len(list_of_genes) < 2:
        return None
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
        organism = "Human"
        table = make_table_from_file('interactions_human.csv')
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
        organism = "Mouse"
        table = make_table_from_file('interactions_mouse.csv')
    normalized_genes = []
    for gene in list_of_genes:
        gene = find_most_known_name(synonyms, gene, organism)
        normalized_genes.append(gene)
    all_pairs = []
    for tf in normalized_genes:
        for tg in normalized_genes:
            if tf is not None and tg is not None:
                searched_interaction = binary_search_table(tf, tg, table, organism)
                if searched_interaction is not None:
                    if detailed_info:
                        all_pairs.append(searched_interaction)
                    else:
                        all_pairs.append((searched_interaction[0], searched_interaction[1],
                                          searched_interaction[2]))

    return all_pairs


def return_info_about_gene(gene: str, is_human: bool) -> Optional[List[str]]:
    """
    Input parameters:
    gene: string, symbol of gene for which user wants to find info
    is_human: boolean, True if organism is human, False if organism is mouse
    Output:
    Returns list of information about gene in format [TaxID, GeneID, Symbol, Synonyms, Description, Type of gene,
    Symbol from nomenclature authority] or None if gene is not among
    most known symbols
    """
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
        gene_table = open('Homo_sapiens.gene_info', 'r')
        gene = find_most_known_name(synonyms, gene, "Human")
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
        gene_table = open('Mus_musculus.gene_info', 'r')
        gene = find_most_known_name(synonyms, gene, "Mouse")

    for line in gene_table:
        line_split = line.split()
        if line_split[2].lower() == gene.lower():
            return [line_split[0], line_split[1], line_split[2], line_split[4],
                    line_split[8], line_split[9], line_split[10]]
    return None


def return_possible_main_names_for_gene(gene: str, is_human: bool) -> List[str]:
    """
    Input parameters:
    gene: string, symbol of gene for which user wants to find all possible genes that can use this symbol
    is_human: boolean, True if organism is human, False if organism is mouse
    Output:
    Returns list of most known names of genes that can possibly use this symbol.
    """
    possible_genes = []
    if is_human:
        synonyms = load_pickle('synonyms_human.pkl')
    else:
        synonyms = load_pickle('synonyms_mouse.pkl')
    for key, values in synonyms.items():
        if gene.lower() == key.lower():
            possible_genes.append(key)
        if gene.lower() in [v.lower() for v in values]:
            possible_genes.append(key)
    return possible_genes


#if __name__ == '__main__':
    # print(return_tgs_for_tf('A2M', True, True))
    # print(return_tfs_for_tg('clint1', False, True))
    # print(return_info_about_two_genes('A2m', 'apod', False))
    # print(return_info_about_group(['ABCF1', 'A2m', 'apod'], False, True))
    # print(return_tfs_for_tg('a2m', True, True))
    # print(return_possible_main_names_for_gene('A2M', True))
    # print(return_tgs_for_tf('NAT10', True, True))
    # print(return_tgs_for_tf('NAT10', True, False))
