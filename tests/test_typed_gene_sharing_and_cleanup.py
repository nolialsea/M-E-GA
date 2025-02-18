"""
test_typed_gene_sharing_and_cleanup.py

Ensures that:
  1) Typed genes remain shared if no mutation occurs (same hash keys across individuals).
  2) Copy-on-mutation spawns new hash keys for mutated individuals.
  3) Unused typed genes are cleaned from encodings after each generation.
"""

import random
import unittest

from src.M_E_GA import M_E_GA_Base

def dummy_fitness(org, ga):
    """
    A trivial fitness function returning 1.0 for all individuals.
    Not relevant to typed gene tests.
    """
    return 1.0

class TestTypedGeneSharingAndCleanup(unittest.TestCase):

    def setUp(self):
        """
        Create a small GA instance that uses typed genes. We'll run short generations
        and forcibly manipulate certain parameters to observe typed gene behavior.
        """
        # We'll seed it for reproducibility.
        random.seed(9999)

        # We'll pass in minimal symbolic genes, rely on typed for the tests
        self.ga = M_E_GA_Base(
            genes=["SymA", "SymB"],
            fitness_function=dummy_fitness,
            population_size=4,     # default 4; we'll do forced mutation with a custom instance
            num_parents=2,
            max_generations=2,
            max_individual_length=2,
            mutation_prob=0.0,     # initially off (we toggle in specific tests)
            logging=False,
            experiment_name="TestTypedGeneSharingAndCleanup"
        )
        self.ga.initialize_population()

    def add_typed_gene_to_population(self, value, range_=(-5, 5), mutation_mode="gaussian"):
        """
        Registers a typed gene in the GA's encoding manager,
        returns its hash key so we can insert it into individuals.
        """
        typed_gene = {
            '__type__': 'numeric',
            'value': value,
            'range': range_,
            'mutation_mode': mutation_mode
        }
        return self.ga.encoding_manager.add_gene(typed_gene)

    def test_no_mutation_shares_typed_genes(self):
        """
        If typed mutation never triggers, every individual referencing the same typed gene
        should keep exactly one shared hash key. No duplication allowed.
        """
        typed_hash = self.add_typed_gene_to_population(value=1.0, range_=(-1,1), mutation_mode="gaussian")

        # Insert typed_hash in each individual's genome
        for org in self.ga.population:
            org.clear()
            org.append(typed_hash)

        # Evaluate & produce next generation with mutation_prob=0.0 => no typed mutation
        fitness_scores = self.ga.population_manager.evaluate_population_fitness(self.ga.population)
        self.ga.population = self.ga.population_manager.select_and_generate_new_population(
            self.ga.population, fitness_scores, generation=0
        )
        # Cleanup step
        self.ga.encoding_manager.cleanup_typed_genes(self.ga.population)

        # Check that they still share the same typed_hash
        hashed_vals = [org[0] for org in self.ga.population]
        unique_hashes = set(hashed_vals)

        self.assertEqual(len(unique_hashes), 1,
                         f"All individuals should still share the same typed gene hash. Got {unique_hashes}."
                         )
        self.assertIn(typed_hash, unique_hashes,
                      f"Expected typed gene hash {typed_hash} to remain in use."
                      )

    def test_forced_mutation_copies_typed_genes(self):
        """
        If typed genes are mutated, each mutated gene should get a new hash key.
        We'll forcibly push typed mutation probability to 1.0 in a separate GA
        with population=2 so we definitely see distinct new hashes.
        """
        # Create a minimal GA with pop=2, ensuring no repeated parent pair combos
        ga2 = M_E_GA_Base(
            genes=["SymA"],
            fitness_function=lambda org, g: 1.0,
            population_size=2,
            max_individual_length=1,   # each organism has exactly 1 gene
            mutation_prob=1.0,         # force typed mutation
            logging=False,
            experiment_name="ForcedMutationCopiesTest"
        )
        ga2.initialize_population()

        # Insert one typed gene, forcing each individual's single gene to be typed_hash
        typed_hash = ga2.encoding_manager.add_gene(
            {
                '__type__': 'numeric',
                'value': 0.0,
                'range': (-1, 1),
                'mutation_mode': 'gaussian'
            }
        )
        for ind in ga2.population:
            ind.clear()
            ind.append(typed_hash)

        # Evaluate & produce next generation => each child should mutate that typed gene
        fitness_scores = ga2.population_manager.evaluate_population_fitness(ga2.population)
        new_pop = ga2.population_manager.select_and_generate_new_population(
            ga2.population, fitness_scores, generation=0
        )

        # Cleanup typed genes
        ga2.encoding_manager.cleanup_typed_genes(new_pop)

        # Expect each new individual's typed hash to be unique
        hashed_vals = [org[0] for org in new_pop]
        unique_hashes = set(hashed_vals)

        self.assertEqual(len(new_pop), 2,
                         "We expect exactly 2 individuals in the new population."
                         )
        self.assertTrue(
            len(unique_hashes) == 2,
            "With forced typed mutation, each individual should get a distinct typed hash key."
        )
        self.assertNotIn(
            typed_hash, unique_hashes,
            "The old typed hash should no longer appear if everything mutated."
        )

    def test_cleanup_removes_unused_typed_genes(self):
        """
        Demonstrate that typed genes not referenced by the new population
        are removed from encodings. We'll:
          1) Insert typed_hash_1 in half the population,
          2) Insert typed_hash_2 in the other half,
          3) Force the new generation to keep only typed_hash_2 => typed_hash_1 becomes unused => cleanup removes it.
        """
        typed_hash_1 = self.add_typed_gene_to_population(value=5.5, range_=(0,10), mutation_mode="uniform")
        typed_hash_2 = self.add_typed_gene_to_population(value=-2.0, range_=(-5,5), mutation_mode="cauchy")

        # Force each half to have only the typed gene we want
        half_size = len(self.ga.population)//2
        for i in range(half_size):
            self.ga.population[i].clear()
            self.ga.population[i].append(typed_hash_1)

        for i in range(half_size, len(self.ga.population)):
            self.ga.population[i].clear()
            self.ga.population[i].append(typed_hash_2)

        # We'll override the fitness function so typed_hash_2 gets top fitness
        def custom_fitness(org, ga):
            return 10.0 if org and org[0] == typed_hash_2 else 0.1

        self.ga.fitness_function = custom_fitness

        # Evaluate & produce next population => presumably all typed_hash_2
        fitness_scores = self.ga.population_manager.evaluate_population_fitness(self.ga.population)
        self.ga.population = self.ga.population_manager.select_and_generate_new_population(
            self.ga.population, fitness_scores, generation=0
        )

        # Cleanup: typed_hash_1 is presumably no longer used
        self.ga.encoding_manager.cleanup_typed_genes(self.ga.population)

        # typed_hash_1 => gone from encodings
        self.assertNotIn(
            typed_hash_1,
            self.ga.encoding_manager.encodings,
            "Unused typed_hash_1 should have been removed by cleanup_typed_genes."
        )
        # typed_hash_2 => still present
        self.assertIn(
            typed_hash_2,
            self.ga.encoding_manager.encodings,
            "typed_hash_2 should remain in encodings, as it's used by the new population."
        )

        # Finally, ensure the new population references typed_hash_2 only
        for org in self.ga.population:
            self.assertEqual(
                org[0], typed_hash_2,
                "All new population members should reference typed_hash_2."
            )


if __name__ == '__main__':
    unittest.main()
