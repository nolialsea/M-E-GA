"""
typed_mutations.py

Handles mutations for typed genes (numeric, numeric_vector, etc.).

We support different mutation modes:
  - Gaussian (existing approach)
  - Uniform
  - Cauchy

Additionally, we handle single-value numeric genes as well as multi-dimensional
vectors. A typed gene is typically a dict with at least:
  {
    '__type__': 'numeric' or 'numeric_vector',
    'value': <float> or 'values': <list of floats>,
    'range': (min, max) or list of (min, max) pairs,
    'mutation_mode': 'gaussian'|'uniform'|'cauchy' (optional)
  }
"""

import math
import random


def perform_typed_mutation(organism, index, generation, manager, default_std_dev=0.1):
    """
    Perform a typed mutation on a gene at the specified index in the organism.

    This function dispatches to the correct sub-operation (gaussian, uniform, cauchy),
    and handles single-value numeric genes or multi-dimensional numeric_vector genes.

    :param organism: The encoded organism (list of hash keys).
    :param index: The position of the typed gene to mutate.
    :param generation: Current generation number (used for logging).
    :param manager: The MutationManager (so we can access GA config, logger, etc.).
    :param default_std_dev: Standard deviation for Gaussian offset if no other info is provided.
    :return: (organism, new_index, mutation_log).
    """
    ga = manager.ga
    hash_key = organism[index]

    typed_gene = ga.encoding_manager.encodings.get(hash_key, None)
    if not typed_gene or not isinstance(typed_gene, dict) or '__type__' not in typed_gene:
        return organism, index + 1, None

    gene_type = typed_gene['__type__']
    if gene_type not in ['numeric', 'numeric_vector']:
        return organism, index + 1, None

    mutation_mode = typed_gene.get('mutation_mode', 'gaussian').lower()

    old_data_repr = _get_data_repr(typed_gene)

    # Decide final_mode
    if mutation_mode == 'gaussian':
        _apply_gaussian_mutation(typed_gene, default_std_dev)
        final_mode = 'gaussian'
    elif mutation_mode == 'uniform':
        _apply_uniform_mutation(typed_gene)
        final_mode = 'uniform'
    elif mutation_mode == 'cauchy':
        _apply_cauchy_mutation(typed_gene)
        final_mode = 'cauchy'
    else:
        # fallback to gaussian
        _apply_gaussian_mutation(typed_gene, default_std_dev)
        final_mode = 'gaussian'

    new_data_repr = _get_data_repr(typed_gene)

    mutation_log = {
        'type': f'typed_{final_mode}_mutation',
        'generation': generation,
        'index': index,
        'old_data': old_data_repr,
        'new_data': new_data_repr
    }
    manager.log_mutation_if_needed(mutation_log)
    return organism, index + 1, mutation_log


def _apply_gaussian_mutation(typed_gene, std_dev=0.1):
    """
    Apply Gaussian mutation to the typed gene (single or multi-dimensional).

    :param typed_gene: The dict describing the typed gene, with 'value' or 'values'.
    :param std_dev: Standard deviation for Gaussian offset.
    """
    # If we have 'values' => multi-dimensional
    if 'values' in typed_gene:
        # Possibly different ranges per dimension or single range
        ranges = typed_gene.get('range', None)
        for i in range(len(typed_gene['values'])):
            old_val = typed_gene['values'][i]
            offset = random.gauss(0, std_dev)
            new_val = old_val + offset
            new_val = _clip_value(new_val, ranges, i)
            typed_gene['values'][i] = new_val
    else:
        # single-value numeric
        old_val = typed_gene['value']
        offset = random.gauss(0, std_dev)
        new_val = old_val + offset
        new_val = _clip_value(new_val, typed_gene.get('range', None))
        typed_gene['value'] = new_val


def _apply_uniform_mutation(typed_gene):
    """
    Apply uniform mutation: set the gene's value(s) to a random sample in [min, max].
    """
    if 'values' in typed_gene:
        ranges = typed_gene.get('range', None)
        for i in range(len(typed_gene['values'])):
            r = _get_range_for_index(ranges, i)
            typed_gene['values'][i] = random.uniform(r[0], r[1])
    else:
        # Single value
        r = typed_gene.get('range', (-999999, 999999))
        typed_gene['value'] = random.uniform(r[0], r[1])


def _apply_cauchy_mutation(typed_gene):
    """
    Apply Cauchy-based mutation, using random.random() for scale.
    We do x_new = x_old + scale * tan(pi*(u-0.5)) approach.

    The scale can be a constant or we can do something fancier later.
    """
    scale = 0.2  # could be a gene param if needed
    if 'values' in typed_gene:
        ranges = typed_gene.get('range', None)
        for i, old_val in enumerate(typed_gene['values']):
            u = random.random()
            offset = scale * math.tan(math.pi * (u - 0.5))
            new_val = old_val + offset
            typed_gene['values'][i] = _clip_value(new_val, ranges, i)
    else:
        old_val = typed_gene['value']
        u = random.random()
        offset = scale * math.tan(math.pi * (u - 0.5))
        new_val = old_val + offset
        typed_gene['value'] = _clip_value(new_val, typed_gene.get('range', None))


def _clip_value(value, ranges, index=0):
    """
    Clip the value to [min, max] using either:
      - A single (min, max) if ranges is a 2-tuple
      - A list of 2-tuples for each dimension
      - If None, no clipping

    :param value: The new value
    :param ranges: The specified range(s)
    :param index: The dimension index if ranges is a list of tuples
    :return: The clipped value
    """
    if ranges is None:
        return value
    if isinstance(ranges[0], (int, float)) and len(ranges) == 2:
        # Single pair
        return max(ranges[0], min(ranges[1], value))
    else:
        # It's presumably a list of (min, max) pairs
        pair = _get_range_for_index(ranges, index)
        return max(pair[0], min(pair[1], value))


def _get_range_for_index(ranges, idx):
    """
    Helper to safely retrieve the range for a given index if 'ranges' is a list, or
    fallback to a default if it's weird.
    """
    if ranges is None:
        return (-999999, 999999)
    if isinstance(ranges[0], (int, float)) and len(ranges) == 2:
        # single pair for the entire vector
        return (ranges[0], ranges[1])
    if idx < len(ranges):
        return ranges[idx]
    return (-999999, 999999)


def _get_data_repr(typed_gene):
    """
    Return a short representation of the typed gene data for logging.

    :param typed_gene: The gene dict
    :return: A float or list of floats representing the data
    """
    if 'values' in typed_gene:
        return typed_gene['values'][:]
    else:
        return typed_gene['value']
