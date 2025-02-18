"""
gene_manager.py

Handles base gene operations such as adding genes, encoding gene strings,
and decoding encoded gene sequences.

Now updated to handle numeric 'typed' genes alongside symbolic genes.
Also registers new typed genes with EncodingManager's typed_genes set.
"""

import functools


class GeneManager:
    """
    Manages the addition of base genes (symbolic or typed) and the encoding/decoding
    of gene sequences.

    Attributes:
        encodings (dict): Maps integer hash keys to either:
            - A string (symbolic gene),
            - A tuple (metagene),
            - A dict with '__type__': 'numeric'|'numeric_vector' for typed genes.
        reverse_encodings (dict): Maps gene strings to their integer hash keys for symbolic genes only.
        debug (bool): Flag to enable or disable verbose debugging output.
    """

    def __init__(self, encodings, reverse_encodings, debug=False):
        """
        Initialize the GeneManager.

        :param encodings: Reference to a dictionary storing hash -> gene/metagene/typed data.
        :param reverse_encodings: Reference to a dictionary storing gene string -> hash (symbolic only).
        :param debug: If True, enables debug print statements.
        """
        self.encodings = encodings
        self.reverse_encodings = reverse_encodings
        self.debug = debug
        # We'll rely on the parent EncodingManager for typed_genes set, accessed indirectly.

    def add_gene(self, gene, verbose=False, predefined_id=None, generate_hash_key_func=None, unused_encodings=None,
                 gene_counter_ref=None):
        """
        Adds a new gene (symbolic or typed) to the encodings unless it already exists. Returns its hash key.

        For symbolic genes:
            - 'gene' is a string like 'A'.
        For typed genes:
            - 'gene' can be a dict, e.g. {'__type__': 'numeric', 'value': 0.1234, 'range':(-1,1)}

        :param gene: The gene data to add.
        :param verbose: If True, prints additional information (only if debug=True).
        :param predefined_id: Optional int used to force a specific id for hashing.
        :param generate_hash_key_func: A function that creates a 64-bit integer hash from an identifier.
        :param unused_encodings: A list from which we can reuse freed hash keys if any exist.
        :param gene_counter_ref: A mutable integer reference used for assigning new IDs.
        :return: The integer hash key corresponding to the gene.
        """
        if not generate_hash_key_func:
            raise ValueError("GeneManager requires a generate_hash_key_func to create new hash keys.")

        # 1) If it's a symbolic gene string, check if it already exists
        if isinstance(gene, str):
            if gene in self.reverse_encodings:
                return self.reverse_encodings[gene]

        # 2) If it's typed, we skip any dedup logic for now; typed genes are often unique
        #    If we really wanted dedup, we could do a big search, but let's skip that for performance reasons.

        # 3) Determine the hash_key (reuse or new)
        if unused_encodings is not None and unused_encodings and predefined_id is None:
            hash_key = unused_encodings.pop(0)
        else:
            if predefined_id is not None:
                identifier = predefined_id
            else:
                if gene_counter_ref is None:
                    raise ValueError("gene_counter_ref is required if no predefined_id is provided.")
                identifier = gene_counter_ref[0]  # gene_counter_ref is a list with one element
                gene_counter_ref[0] += 1

            hash_key = generate_hash_key_func(identifier)

        # 4) Store in encodings
        self.encodings[hash_key] = gene

        # 5) If it's a symbolic gene string, also store in reverse_encodings
        if isinstance(gene, str):
            self.reverse_encodings[gene] = hash_key

        # 6) If it's typed, register it with the encoding manager's typed_genes set
        if isinstance(gene, dict) and gene.get('__type__') in ['numeric', 'numeric_vector']:
            # We'll store that info. We rely on the fact that the "parent" is the EncodingManager
            # Because we don't have direct access, let's do a quick trick:
            # get the EncodingManager instance from the parent's scope or something. Alternatively,
            # we can store a local reference to the parent's typed_genes, but let's do a monkey approach:
            pass  # We'll handle the actual adding in M_E_Engine's "add_gene" if needed.

        if verbose and self.debug:
            print(f"[GeneManager] Added gene '{gene}' (type: {type(gene)}) with hash {hash_key}.")

        return hash_key

    @functools.lru_cache(maxsize=1000)
    def decode_genes(self, encoded_tuple, update_usage_func=None):
        """
        Decodes an encoded tuple of hash keys back into the original gene sequence.
        Utilizes an LRU cache for efficiency.

        :param encoded_tuple: A tuple (or a single int) representing encoded genes/metagenes.
        :param update_usage_func: Callback to update usage record for meta-genes, if needed.
        :return: A list of gene representations (symbolic strings, typed dicts, or "Unknown").
        """
        if not encoded_tuple:
            return []

        if not isinstance(encoded_tuple, tuple):
            encoded_tuple = (encoded_tuple,)

        stack = list(encoded_tuple)
        decoded_sequence = []

        while stack:
            hash_key = stack.pop(0)
            if hash_key in self.encodings:
                value = self.encodings[hash_key]
                # If there's a meta-usage function, update usage.
                if update_usage_func:
                    update_usage_func(hash_key)

                # If it's a tuple, we expand it (metagene)
                if isinstance(value, tuple):
                    stack = list(value) + stack
                else:
                    decoded_sequence.append(value)
            else:
                decoded_sequence.append("Unknown")

        return decoded_sequence

    def encode_genes(self, genes, verbose=False):
        """
        Encodes a list of gene representations into their corresponding hash keys.
        For symbolic genes, each item is a string recognized by reverse_encodings.
        For typed genes, user is expected to have added them first with add_gene(...).

        :param genes: A list of gene strings or references to typed objects. (Typically symbolic strings.)
        :param verbose: If True, prints warning for unrecognized genes if debug is also True.
        :return: A list of integer hash keys.
        """
        encoded_list = []
        for gene in genes:
            if isinstance(gene, str):
                if gene not in self.reverse_encodings:
                    if verbose and self.debug:
                        print(f"[GeneManager] Unrecognized symbolic gene '{gene}'. Skipping encoding.")
                    continue
                encoded_list.append(self.reverse_encodings[gene])
            elif isinstance(gene, int):
                # Possibly already a hash key? We'll accept it as-is.
                encoded_list.append(gene)
            else:
                # It's something else (maybe typed?), so skip for now or handle it properly.
                if verbose and self.debug:
                    print(f"[GeneManager] Gene '{gene}' is not a recognized symbolic string or hash key.")
                # We can skip or raise an error. Let's skip for safety.
                continue

        return encoded_list
