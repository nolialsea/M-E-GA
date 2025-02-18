# Typed Genes in MEGA

With MEGA, genes aren’t limited to simple symbolic strings like `"A"`, `"B"`, etc. We also support **typed** genes that store numeric data—handy for a broad range of tasks like optimizing neural network weights, numerical parameters, or multi-dimensional vectors.

## Overview

Typed genes can either be:
- **Single-Value Numeric Genes**:  
  ```python
  {
    "__type__": "numeric",
    "value": 0.5,
    "range": (-1.0, 1.0),
    "mutation_mode": "gaussian"  # or "uniform", "cauchy"
  }
  ```
- **Multi-Dimensional Numeric Vector Genes**:
  ```python
  {
    "__type__": "numeric_vector",
    "values": [0.1, -0.2, 0.75],
    "range": [(-1.0, 1.0), (-2.0, 2.0), (0.5, 1.0)],
    "mutation_mode": "cauchy"
  }
  ```

### Key Points

1. **Mutation Modes**:
  - **Gaussian**: Offsets the gene’s current value(s) by `random.gauss(0, std_dev)`. Default `std_dev` is 0.1 if none is specified.
  - **Uniform**: Randomly picks a new value within the specified `[min, max]` range for each dimension.
  - **Cauchy**: Uses a Cauchy distribution to perturb the gene’s value(s).
  - **Fallback**: If a typed gene uses an unrecognized `mutation_mode`, MEGA automatically uses Gaussian.

2. **Ranges**:
  - **Single Pair**:  `(min_value, max_value)` for all values in a vector.
  - **Per-Dimension**: A list of `(min, max)` pairs, each matching a dimension in `values`.
  - If no `range` is specified, the system defaults to a big bounding box `(-999999, 999999)`, which effectively means “almost anything goes.”

3. **Coexisting with Symbolic Genes**:
  - You can freely mix symbolic and typed genes in the same organism.
  - The mutation engine automatically detects typed genes and applies the numeric mutation logic. Symbolic ones get normal point/swap/insertion/deletion.

4. **Integration**:
  - **Adding a Typed Gene**: Use `ga.encoding_manager.add_gene(...)` with a dict describing your typed gene.
  - **Mutation**: The GA loop calls the `MutationManager`, which detects typed genes and mutates them accordingly.
  - **Decoding**: `ga.encoding_manager.decode(...)` returns a dict for typed genes instead of a string. That way, your fitness function can do numeric computations, or re-encode them if needed.

## Example

```python
# Suppose we have a GA instance
my_ga = M_E_GA_Base(
    genes=["SymA", "SymB"],          # Some symbolic genes
    fitness_function=your_fitness_fn,
    population_size=20,
    max_generations=50
)
my_ga.initialize_population()

# Define a single-value numeric gene
num_gene = {
    "__type__": "numeric",
    "value": 2.5,
    "range": (0, 5),
    "mutation_mode": "gaussian"
}
# Add it to the manager, get the hash key
num_gene_hash = my_ga.encoding_manager.add_gene(num_gene)

# Build an organism that has both symbolic and typed genes
symA_hash = my_ga.encoding_manager.reverse_encodings["SymA"]
organism = [symA_hash, num_gene_hash]

# Mutate
mutated_org = my_ga.mutation_manager.mutate_organism(organism, generation=0)
decoded = my_ga.encoding_manager.decode(tuple(mutated_org))
print("Decoded with typed genes:", decoded)
# -> You’ll see something like ["SymA", {"__type__": "numeric", "value": 2.47, "range": (0, 5), ...}]

# All standard GA flows apply the same way with typed genes in the population.
```

---

### Gotchas and Tips

- Always ensure your typed gene dictionaries have `{"__type__": "numeric"}` or `{"__type__": "numeric_vector"}` so the GA knows how to handle them.
- If you supply `range` as a list of tuples for a numeric_vector, make sure the list length matches the length of `values`.
- If you do something truly “out there” with typed genes, be prepared for it to mutate in unexpected ways. The code is flexible but also expects well-structured data.

That’s it! Typed genes in MEGA let you blend symbolic + numeric evolution in a single framework.
