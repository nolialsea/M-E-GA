# ROADMAP

## 1. Typed Genes

- **Done**: We introduced typed genes (`numeric` and `numeric_vector`) with multiple mutation modes (gaussian, uniform, cauchy).
- **Done**: Implemented a **copy-on-mutation** strategy so that typed genes remain shared if unmutated, but instantly diverge once mutation hits.
- **Done**: Added a cleanup step after each generation to remove any typed genes no longer referenced by the population, preventing a memory blow-up.
- **Done**: Created thorough test coverage (`test_typed_genes.py` + `test_typed_gene_sharing_and_cleanup.py`) to ensure typed gene sharing, forced mutation, and cleanup are all working correctly.

### Possible Future Enhancements

- Experiment with partial/hierarchical typed gene sharing or specialized caching for large numeric-based populations.
- Consider advanced typed gene mutation strategies (e.g., per-dimension adaptive rates).

---

## 2. Population Manager Updates

- **Done**: Allowed `max_individual_length=1` by adjusting the random range selection in population initialization, so we can handle single-gene organisms without empty-range errors.

### Potential Follow-Ups

- Provide even more flexible ways to define population initialization bounds, especially if you need zero-length organisms or special-case constraints.
- Explore domain-specific initialization strategies (like near-limits for typed genes) to jumpstart certain solutions.

---

## 3. Additional Ideas

- Investigate synergy between typed genes and meta-genes: how capturing or opening nested segments interacts with numeric data.
- Evaluate performance profiling for large-scale runs with heavy typed gene usage. Possibly offload part of typed gene management to GPU if you’re doing huge numeric vectors.
- Keep refining logging, real-time event hooking, and analytics integration for deeper insights into typed gene evolutions.  

