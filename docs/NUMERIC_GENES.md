# Using Numeric (Typed) Genes

In addition to symbolic string genes, MEGA supports “typed” genes where each gene can store numeric data,
like floats for neural network weights or parameter optimization.

### Adding a Numeric Gene
To create a numeric gene, you can do:

```python
from src.M_E_GA import M_E_GA_Base

# Create GA instance with some symbolic genes
ga = M_E_GA_Base(genes=['SymA', 'SymB'], fitness_function=your_fitness_func)

# Initialize population so we have a working EncodingManager
ga.initialize_population()

# Now define a numeric gene with a value and range
my_numeric_gene = {
    '__type__': 'numeric',
    'value': 0.5,       # The current numeric value
    'range': (-2.0, 2.0)# (optional) range for clipping
}

# Add it to the EncodingManager
hash_key = ga.encoding_manager.add_gene(my_numeric_gene)
print("Created numeric gene with hash:", hash_key)

# We can create an organism that uses both symbolic and numeric genes:
symb_hash = ga.encoding_manager.reverse_encodings['SymA']
organism = [symb_hash, hash_key]

# Mutate organism
mutated_org = ga.mutation_manager.mutate_organism(organism, generation=0)

# Decoding yields either a string for symbolic genes or a dict for numeric:
decoded = ga.encoding_manager.decode(tuple(mutated_org))
print("Decoded organism:", decoded)
```

### Mutation Behavior

- If the gene is identified as `{'__type__': 'numeric', ...}`, MEGA will apply a Gaussian-based mutation
  (`perform_numeric_gaussian_mutation`) with some default `std_dev` (0.1).
- The new value will be clipped into the specified `range` if provided.
- You can tweak `mutation_prob` on the GA to control the chance of any mutation per gene.

### Coexistence with Symbolic Genes

It’s perfectly fine to have an organism that contains a mix of numeric typed genes and symbolic string genes.  
The normal symbolic mutation logic (point, swap, etc.) coexists with numeric-specific mutation.
