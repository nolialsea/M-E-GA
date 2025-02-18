"""
test_typed_genes.py

Unit tests to ensure typed numeric genes are properly handled:
 - Single-value and multi-dimensional (vector) numeric genes
 - Various mutation modes: gaussian, uniform, cauchy
 - Ensuring decode doesn't break
 - Checking that typed genes can coexist with symbolic ones
 - Handling edge cases like per-dimension ranges, missing ranges, unknown modes, etc.
"""

import random
import unittest

from src.M_E_GA import M_E_GA_Base


class TestTypedGenes(unittest.TestCase):
    def setUp(self):
        """
        Initialize a minimal GA instance with a tiny gene set (symbolic),
        then we add numeric typed genes on the fly to see if they integrate well.
        """
        random.seed(20250)
        self.ga = M_E_GA_Base(
            genes=['A', 'B'],  # just a couple symbolic genes
            fitness_function=lambda org, ga: 0,  # not relevant for these tests
            population_size=2,
            max_generations=1,
            logging=False,
            experiment_name="TestTypedGenes"
        )
        # Make sure the GA is initialized so we have a functional encoding manager
        self.ga.initialize_population()

    def test_add_numeric_gene(self):
        """
        Ensure we can add a single-value numeric typed gene and retrieve its hash key.
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
        Add a numeric gene to an organism, run a mutation pass, and see if
        the value changes as expected.
        """
        numeric_gene = {
            '__type__': 'numeric',
            'value': 0.0,
            'range': (-1, 1),
            'mutation_mode': 'gaussian'
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_gene)

        organism = [hash_key]

        # Force 100% normal mutation rate to ensure we do mutate
        self.ga.mutation_prob = 1.0
        self.ga.delimited_mutation_prob = 0
        self.ga.delimiter_insert_prob = 0
        self.ga.open_mutation_prob = 0
        self.ga.metagene_mutation_prob = 0

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        self.assertEqual(len(mutated), 1, "Length should remain 1 after numeric mutation.")

        mutated_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[mutated_key]
        new_val = mutated_gene['value']
        self.assertNotEqual(new_val, 0.0,
                            "The numeric gene value should have been mutated from 0.0 to something else.")
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
        numeric_gene = {
            '__type__': 'numeric',
            'value': 0.5,
            'range': (0, 1)
        }
        numeric_key = self.ga.encoding_manager.add_gene(numeric_gene)
        symbolic_key = self.ga.encoding_manager.reverse_encodings['A']

        organism = [symbolic_key, numeric_key]
        self.ga.mutation_prob = 1.0

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        self.assertTrue(len(mutated) >= 1,
                        "Mutated organism should have at least 1 gene after mutation.")
        decoded = self.ga.encoding_manager.decode(tuple(mutated))
        has_numeric = any(isinstance(x, dict) and x.get('__type__') == 'numeric' for x in decoded)
        self.assertTrue(has_numeric, "Decoded organism should still contain at least one numeric gene dict.")

    def test_numeric_vector_uniform(self):
        """
        Test a multi-dimensional typed gene with uniform mutation.
        """
        numeric_vector_gene = {
            '__type__': 'numeric_vector',
            'values': [0.0, 0.0, 0.0],
            'range': (-2, 2),
            'mutation_mode': 'uniform'
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_vector_gene)
        organism = [hash_key]

        self.ga.mutation_prob = 1.0
        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        new_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[new_key]
        self.assertIn('values', mutated_gene, "After uniform mutation, 'values' should still exist.")
        for val in mutated_gene['values']:
            self.assertTrue(-2 <= val <= 2, "Each value after uniform mutation must be within [-2, 2].")

    def test_numeric_vector_cauchy(self):
        """
        Test cauchy mutation for a numeric_vector typed gene.
        """
        numeric_vector_gene = {
            '__type__': 'numeric_vector',
            'values': [1.0, -1.0],
            'range': (-10, 10),
            'mutation_mode': 'cauchy'
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_vector_gene)
        organism = [hash_key]

        self.ga.mutation_prob = 1.0
        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        new_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[new_key]
        self.assertIn('values', mutated_gene, "Should still have 'values' after cauchy mutation.")
        self.assertEqual(len(mutated_gene['values']), 2, "We still have 2D data.")
        # We won't do a deep check, but let's see if they're in range
        for val in mutated_gene['values']:
            self.assertTrue(-10 <= val <= 10, "Values should remain in the specified range.")

    # ------------------------------------------------------------------------
    # Additional coverage tests
    # ------------------------------------------------------------------------

    def test_numeric_vector_per_dimension_range(self):
        """
        Verify that if we have a per-dimension range, each dimension gets clipped correctly.
        For instance, range = [(-1,1), (5,6)] => the first dimension is in [-1,1], second in [5,6].
        """
        numeric_vector_gene = {
            '__type__': 'numeric_vector',
            'values': [0.5, 5.5],
            'range': [(-1, 1), (5, 6)],
            'mutation_mode': 'uniform'
        }
        hash_key = self.ga.encoding_manager.add_gene(numeric_vector_gene)
        organism = [hash_key]

        self.ga.mutation_prob = 1.0

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        mutated_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[mutated_key]

        self.assertIn('values', mutated_gene)
        self.assertEqual(len(mutated_gene['values']), 2)
        self.assertTrue(-1 <= mutated_gene['values'][0] <= 1,
                        "Dimension 0 should remain within [-1,1].")
        self.assertTrue(5 <= mutated_gene['values'][1] <= 6,
                        "Dimension 1 should remain within [5,6].")

    def test_unknown_mutation_mode_fallback(self):
        """
        If a typed gene has an unrecognized 'mutation_mode', we fall back to Gaussian.
        We'll confirm that the mutation event logs as 'typed_gaussian_mutation'.
        """
        # Make sure we definitely mutate:
        self.ga.mutation_prob = 1.0

        weird_gene = {
            '__type__': 'numeric',
            'value': 0.0,
            'range': (-1, 1),
            'mutation_mode': 'banana'
        }
        hash_key = self.ga.encoding_manager.add_gene(weird_gene)
        organism = [hash_key]

        mutated_org, logs = self.ga.mutation_manager.mutate_organism(organism, generation=0, log_enhanced=True)
        self.assertGreater(len(logs), 0, "We should have at least one mutation event in the logs.")

        # The event type should reflect the fallback to Gaussian
        event_type = logs[-1]['mutation_event']['type']
        self.assertIn('gaussian', event_type,
                      "Unknown mode fallback should produce a 'typed_gaussian_mutation' event type.")

    def test_missing_range_omitted(self):
        """
        If 'range' is missing, the typed gene is defaulted to no clipping, i.e. (-999999, 999999).
        We'll confirm it doesn't crash and that mutation is possible.
        """
        self.ga.mutation_prob = 1.0

        no_range_gene = {
            '__type__': 'numeric',
            'value': 99.0,
            'mutation_mode': 'uniform'
        }
        hash_key = self.ga.encoding_manager.add_gene(no_range_gene)
        organism = [hash_key]

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        mutated_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[mutated_key]

        self.assertNotEqual(mutated_gene['value'], 99.0,
                            "We used uniform without a range, value should definitely have changed from 99.")

    def test_boundary_values_clipping(self):
        """
        Ensure that large offsets near boundary get clipped properly.
        We'll set the initial value near the boundary and do a forced large offset.
        """
        self.ga.mutation_prob = 1.0

        boundary_gene = {
            '__type__': 'numeric',
            'value': 1.95,  # near the top of our -2..2
            'range': (-2, 2),
            'mutation_mode': 'gaussian'
        }
        hash_key = self.ga.encoding_manager.add_gene(boundary_gene)
        organism = [hash_key]

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        mutated_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[mutated_key]

        self.assertTrue(-2 <= mutated_gene['value'] <= 2,
                        f"Value after big gaussian offset should still be clipped within [-2,2].")

    def test_empty_vector_values(self):
        """
        Check if a 'numeric_vector' typed gene with an empty 'values' array
        can be mutated without crashing. We expect it to remain empty.
        """
        self.ga.mutation_prob = 1.0

        empty_vector_gene = {
            '__type__': 'numeric_vector',
            'values': [],
            'range': (-1, 1),
            'mutation_mode': 'gaussian'
        }
        hash_key = self.ga.encoding_manager.add_gene(empty_vector_gene)
        organism = [hash_key]

        mutated = self.ga.mutation_manager.mutate_organism(organism, generation=0)
        new_key = mutated[0]
        mutated_gene = self.ga.encoding_manager.encodings[new_key]
        self.assertIn('values', mutated_gene, "Still should have the 'values' field.")
        self.assertEqual(len(mutated_gene['values']), 0,
                         "Empty vector should remain empty after mutation (no data to mutate).")


if __name__ == '__main__':
    unittest.main()
