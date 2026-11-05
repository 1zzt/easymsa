# MSA-Kit benchmark variants

Three experimental DNA/RNA aligners derived from malabz/MSA-Kit.

- `main` / `codex/fast-dedup`: compact traceback and duplicate-sequence caching.
- `codex/kmer-center`: normalized 4-mer profile center selection.
- `codex/affine-star`: affine gap penalties.

Sources, upstream MIT license, offline dependencies, and independent tests are included.
The programs were built and tested on a CPU Linux server. Synthetic local measurements are not EasyMSA official benchmark scores.

Upstream: https://github.com/malabz/MSA-Kit/tree/bbd5767a250654cb6b0c494444405e74bd115f02
