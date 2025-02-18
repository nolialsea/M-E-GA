# Meta-Genes in MEGA

This document provides an overview of Meta-Genes, a critical feature in the Mutable Encoding Enabled Genetic Algorithm (MEGA).

## 1. What Are Meta-Genes?

Meta-Genes are special “genes” that store entire genetic segments inside them. Instead of representing a single symbolic or typed gene, a Meta-Gene is essentially a compressed tuple of genes. This allows segments of the genome to be treated as a single unit, letting them be passed around, opened, or mutated in chunk form.

- **Compressing** a segment is called _capturing_ it. 
- **Expanding** that Meta-Gene back is called _opening_ it.

When a Meta-Gene is present in an organism, the GA sees it as just another gene (albeit a special one). However, if you decide to “open” it (expand), it unfolds into its constituent genes (with optional delimiters around them).

## 2. Capturing (Compressing) Segments

MEGA can detect delimited blocks (`Start` ... `End`) and wrap them up into a single Meta-Gene. For instance:
```
[Start, A, B, End]
```
might turn into:
```
[MetaGene123]
```
internally. This single integer codon references the **tuple** `(A, B)` in the encodings dictionary.  

### Why Capture?

1. **Speed up** repeated patterns by storing them once.  
2. **Make nested structures** feasible, as you can have Meta-Genes inside other Meta-Genes (like Russian dolls).  
3. **Prevent an explosion** of repeated sequences. 

## 3. Opening (Expanding) Meta-Genes

When you open a Meta-Gene, you replace that single codon with the stored sequence. You can either:
- Wrap it in `Start`/`End` delimiters again (the default).
- Or expand it “raw” (without delimiters) if `no_delimit=True`.

Example:

```
Before: [MetaGene123]
Open:   [Start, A, B, End]
```

or, if `no_delimit=True`:
```
Open:   [A, B]
```

## 4. LRU Usage & Deletion

Meta-Genes are tracked using an LRU (Least Recently Used) system to avoid memory bloat. Each time you reference a Meta-Gene, it’s marked as “used.” If a Meta-Gene stays unused for too many generations, MEGA automatically “deletes” it:
- **Deletion** means the stored tuple is inlined back into any higher-level Meta-Genes that referenced it (so references stay consistent).
- Freed references can be reused by new Meta-Genes.

## 5. Why Meta-Genes?

By capturing repeated or meaningful patterns into single “units,” the GA can:
- Explore higher-level building blocks of solutions more effectively.
- Retain learned structures across generations or even runs (transfer learning).
- Support advanced nested representations for complex tasks.

In a nutshell, Meta-Genes break the mold of “one gene = one simple token” and let MEGA treat segments as composable building blocks. They’re one of the core reasons MEGA stands out from typical GAs.

## 6. Further Reading

- Check out `TYPED_GENES.md` for how typed genes complement the Meta-Gene approach.
- The `M_E_Engine.py` and `meta_gene_manager.py` files in `src/M_E_GA/engine` hold the main logic for capturing, opening, and garbage-collecting Meta-Genes.
- Don’t forget to run the test suite (`test_M_E_Engine.py`, `test_mutation_manager.py`, etc.) to see how Meta-Genes are tested.

Enjoy the meta-ness, folks!
