"""
typed_mutations.py

Handles mutations for typed genes (numeric, numeric_vector, etc.) using a
copy-on-mutation approach to avoid blowing up memory. We only clone the gene
when we actually mutate it, letting multiple individuals share typed genes
if none of them mutate those genes.

We support:
  - Gaussian
  - Uniform
  - Cauchy

Additionally handle single-value or multi-dimensional typed genes.

When we detect a typed gene is about to mutate, we:
  1) Clone the gene data (deep copy)
  2) Insert the clone as a new entry in the encodings, grabbing a new hash key
  3) Mutate the cloned version
  4) Replace the organism's codon with the new hash key
"""

import copy
import math
import random


def perform_typed_mutation(organism, index, generation, manager, default_std_dev=0.1):
    """
    Perform a typed mutation on the gene at the specified index in the organism.

    Using copy-on-mutation, we do:
      1) Clone the typed gene data so we don't mutate the shared reference
      2) Insert the cloned data as a new gene in encodings (new hash)
      3) Apply the chosen numeric mutation method

    :param organism: The encoded organism (list of hash keys).
    :param index: The position of the typed gene to mutate.
    :param generation: Current generation number (used for logging).
    :param manager: The MutationManager (so we can access GA config, logger, etc.).
    :param default_std_dev: Standard deviation for Gaussian offset if no other info is provided.
    :return: (organism, new_index, mutation_log).
    """
    ga = manager.ga
    old_hash_key = organism[index]

    # 1) Retrieve the existing typed gene data
    old_typed_gene = ga.encoding_manager.encodings.get(old_hash_key, None)
    if (not old_typed_gene) or not isinstance(old_typed_gene, dict) or '__type__' not in old_typed_gene:
        # Not a recognized typed gene; skip
        return organism, index + 1, None

    # 2) Clone the typed gene data so we don't mutate the shared original
    cloned_gene = copy.deepcopy(old_typed_gene)
    # Insert as a new gene => new hash key
    new_hash_key = ga.encoding_manager.add_gene(cloned_gene)
    # Update organism's codon to the new hash key
    organism[index] = new_hash_key

    # 3) Decide which mutation mode to apply
    mutation_mode = cloned_gene.get('mutation_mode', 'gaussian').lower()

    old_data_repr = _get_data_repr(old_typed_gene)
    # We'll mutate the cloned_gene in place
    if mutation_mode == 'gaussian':
        _apply_gaussian_mutation(cloned_gene, default_std_dev)
        final_mode = 'gaussian'
    elif mutation_mode == 'uniform':
        _apply_uniform_mutation(cloned_gene)
        final_mode = 'uniform'
    elif mutation_mode == 'cauchy':
        _apply_cauchy_mutation(cloned_gene)
        final_mode = 'cauchy'
    else:
        # fallback to gaussian
        _apply_gaussian_mutation(cloned_gene, default_std_dev)
        final_mode = 'gaussian'

    new_data_repr = _get_data_repr(cloned_gene)

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
        r = typed_gene.get('range', (-999999, 999999))
        typed_gene['value'] = random.uniform(r[0], r[1])


def _apply_cauchy_mutation(typed_gene):
    """
    Apply Cauchy-based mutation, using random.random() for scale.
    We do: x_new = x_old + scale * tan(pi*(u-0.5)).

    :param typed_gene: The dict describing the typed gene.
    """
    scale = 0.2  # could be configurable
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
      - A list of (min, max) pairs for each dimension
      - If None, no clipping (defaults to -999999..999999)

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
    Helper to retrieve the range for a given index if 'ranges' is a list,
    or fallback to a big bounding box.

    :param ranges: The range data, either a 2-tuple or a list of 2-tuples
    :param idx: The dimension index
    :return: (min_val, max_val)
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
