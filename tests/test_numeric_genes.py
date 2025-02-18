"""
test_numeric_genes.py

Unit tests to ensure numeric (typed) genes are properly handled:
 - Adding numeric genes
 - Mutating them with Gaussian noise
 - Ensuring decode doesn't break
 - Checking that numeric genes can coexist with symbolic ones
"""

import random
import unittest

from src.M_E_GA import M_E_GA_Base


class TestNumericGenes(unittest.TestCase):
    def setUp(self):
        """
        Initialize a minimal GA instance with a tiny gene set (symbolic),
        and we'll add numeric genes on the fly to see if they integrate well.
        """
        random.seed(20250)
        self.ga = M_E_GA_Base(
            genes=['A', 'B'],  # just a couple symbolic genes
            fitness_function=lambda org, ga: 0,  # not relevant for these tests
            population_size=2,
            max_generations=1,
            logging=False,
            experiment_name="TestNumericGenes"
        )
        # Make sure the GA is initialized so we have a functional encoding manager
        self.ga.initialize_population()

    def test_add_numeric_gene(self):
        """
        Ensure we can add a numeric typed gene and retrieve its hash key.
        """
        numeric_gene = {
            '__type__': 'numeric',
            'value': 3.14,
            'range': (-5, 5)
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_gene)
        self.assertIn(hash_key, self.ga.encoding_manager.encodings,
                      "Hash key for numeric gene should exist in encodings.")
        stored = self.ga.encoding_manager.encodings[hash_key]
        self.assertEqual(stored['value'], 3.14, "Numeric gene's value should match what we added.")
        self.assertEqual(stored['range'], (-5, 5), "Numeric gene's range should be stored correctly.")

    def test_mutation_on_numeric_gene(self):
        """
        Add a numeric gene to an organism, run a mutation pass, and see if the
        value changes as expected (Gaussian offset).
        """
        # Add numeric gene
        numeric_gene = {
            '__type__': 'numeric',
            'value': 0.0,
            'range': (-1, 1)
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_gene)

        # Build a tiny organism: [hash_key]
        organism = [hash_key]

        # Force 100% normal mutation rate to ensure we mutate
        self.ga.mutation_prob = 1.0
        # Turn off everything else
        self.ga.delimited_mutation_prob = 0
        self.ga.delimiter_insert_prob = 0
        self.ga.open_mutation_prob = 0
        self.ga.metagene_mutation_prob = 0

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        # mutated organism should be same length => [hash_key], but let's see if the value changed
        self.assertEqual(len(mutated), 1, "Length should remain 1 after numeric mutation.")
        mutated_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[mutated_key]
        old_val = 0.0
        new_val = mutated_gene['value']
        self.assertNotEqual(new_val, old_val,
                            "The numeric gene value should have been mutated from 0.0 to something else.")
        # Also check if it's in the range
        self.assertTrue(-1 <= new_val <= 1,
                        f"New value {new_val} should be clipped within the range [-1,1].")

    def test_decode_numeric_gene(self):
        """
        Ensure decode returns a dict for typed numeric genes, not 'Unknown'.
        """
        numeric_gene = {
            '__type__': 'numeric',
            'value': 2.71828,
            'range': (0, 10)
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_gene)
        decoded = self.ga.encoding_manager.decode((hash_key,))
        self.assertEqual(len(decoded), 1, "Decoding a single-gene tuple should yield 1 item.")
        self.assertIsInstance(decoded[0], dict, "Decoded numeric gene should be a dict, not a string.")
        self.assertIn('__type__', decoded[0], "Decoded numeric gene dict should contain '__type__' key.")
        self.assertEqual(decoded[0]['value'], 2.71828,
                         "Decoded numeric gene 'value' should match the original.")

    def test_coexist_symbolic_and_numeric(self):
        """
        Check if an organism can contain both symbolic and numeric genes,
        and the mutation logic doesn't blow up.
        """
        # Add a numeric gene
        numeric_gene = {
            '__type__': 'numeric',
            'value': 0.5,
            'range': (0, 1)
        }
        numeric_key = self.ga.encoding_manager.add_gene(numeric_gene)

        # symbolic 'A' is already in reverse_encodings
        symbolic_key = self.ga.encoding_manager.reverse_encodings['A']

        organism = [symbolic_key, numeric_key]
        # Force a moderate mutation rate
        self.ga.mutation_prob = 1.0

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        self.assertTrue(len(mutated) >= 1,
                        "Mutated organism should have at least 1 gene (but could have more if insertion happened).")

        # Just ensure it doesn't crash, and see if numeric gene changed
        # We'll do a quick decode to see if there's a dict
        decoded = self.ga.encoding_manager.decode(tuple(mutated))
        has_numeric = any(isinstance(x, dict) and x.get('__type__') == 'numeric' for x in decoded)
        self.assertTrue(has_numeric, "Decoded organism should still contain at least one numeric gene dict.")


if __name__ == '__main__':
    unittest.main()
