# ROADMAP

## 1. Typed Genes (Partially Implemented)
- We added basic support for storing single-value numeric typed genes in `GeneManager`.
- We introduced a simple Gaussian mutation approach to mutate numeric genes.
- Next Steps:
   - Implement multi-dimensional typed genes (vectors) for neural net weights or advanced data structures.
   - Add more numeric mutation types (Cauchy, uniform, etc.).
   - Possibly handle typed meta-genes more explicitly.

## 2. Further Enhancements
- Enhance searching/deduplication for typed genes if needed.
- Potentially separate typed gene logic into its own manager if it grows complicated.
- Keep an eye on performance overhead with large numeric-based populations.
