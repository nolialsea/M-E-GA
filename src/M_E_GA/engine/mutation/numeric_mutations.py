"""
numeric_mutations.py

Handles mutations specifically for numeric typed genes, applying
random changes (e.g., Gaussian, uniform, or any other approach).

Make sure to import and use as needed in the MutationManager.
"""

import random


def perform_numeric_gaussian_mutation(organism, index, generation, manager, std_dev=0.1):
    """
    Perform a Gaussian-based mutation on a numeric typed gene.

    :param organism: The encoded organism (list of hash keys).
    :param index: The position of the gene to mutate.
    :param generation: Current generation number (used for logging).
    :param manager: The MutationManager (so we can log, etc.).
    :param std_dev: Standard deviation for the Gaussian offset.
    :return: (organism, new_index, mutation_log)
    """
    ga = manager.ga
    hash_key = organism[index]

    # Retrieve the typed gene dict
    typed_gene = ga.encoding_manager.encodings.get(hash_key, None)
    if not typed_gene or typed_gene.get('__type__') != 'numeric':
        return organism, index + 1, None

    old_value = typed_gene['value']
    gene_range = typed_gene.get('range', (-999999, 999999))

    # Generate offset from normal distribution
    offset = random.gauss(0, std_dev)
    new_value = old_value + offset
    # Clip if range is given
    new_value = max(gene_range[0], min(gene_range[1], new_value))

    # Update the typed gene
    typed_gene['value'] = new_value

    mutation_log = {
        'type': 'numeric_gaussian_mutation',
        'generation': generation,
        'index': index,
        'old_value': old_value,
        'new_value': new_value,
        'std_dev': std_dev
    }
    manager.log_mutation_if_needed(mutation_log)
    return organism, index + 1, mutation_log
