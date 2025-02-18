# ROADMAP

## 1. Typed Genes (Partially Implemented)

- **Done**: We introduced multi-dimensional typed genes (`numeric_vector`), plus `uniform` and `cauchy` mutations.
- Next Steps:
    - Possibly handle typed meta-genes more explicitly.
    - Possibly store dimension-specific ranges more elegantly.
    - Investigate advanced typed gene mutation strategies or hybrid approaches.

## 2. Further Enhancements

- Possibly enhance searching/deduplication for typed genes if needed.
- Potentially separate typed gene logic into its own manager if it becomes more complex.
- Keep an eye on performance overhead with large numeric-based or vector-based populations.
