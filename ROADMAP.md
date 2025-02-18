# ROADMAP

## 1. Typed Genes
- **Done**: Introduced multi-dimensional typed genes, plus multiple mutation modes (gaussian, uniform, cauchy).
- **Done**: Implemented copy-on-mutation so typed genes remain shared unless a child mutates them, preventing massive duplication.
- **Done**: Cleanup method to remove typed genes unused by any current organism, avoiding memory leaks.

## 2. Further Enhancements
- Optionally refine typed gene sharing vs. forced copying.
- Possibly integrate advanced typed gene mutation strategies or a specialized typed gene manager.
- Keep an eye on large populations with heavy typed gene usage.
