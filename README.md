# MSA-Kit website test: fast_dedup

Derived from malabz/MSA-Kit commit bbd5767a250654cb6b0c494444405e74bd115f02.
This is an experimental DNA/RNA program for testing the EasyMSA journey.
It is not a claim of superiority over established aligners.

Only the alignment core and version label are changed. Upstream MIT and
third-party notices are retained. Build/run entrypoints are in benchmark/.

Variants: fast_dedup uses rolling score rows, byte traceback and duplicate
sequence caching with upstream scoring/ties; kmer_center additionally selects
the sequence nearest the mean normalized 4-mer profile (length penalty 0.05);
affine_star uses the first center, match +2, mismatch -1, gap open 3 and extend 1.
All retain N-neutral and T/U-equivalent scoring, input order, original symbols
through the upstream benchmark adapter, and a single CPU execution thread.

Run `bash benchmark/build.sh`, then
`bash benchmark/run.sh INPUT_FASTA OUTPUT_FASTA 1 dna` (or rna).
OUTPUT_FASTA must be new, and its parent directory must exist.

Core tests: `cmake -S . -B build -DBUILD_TESTING=ON`,
`cmake --build build --parallel 2`, `ctest --test-dir build -L contract --output-on-failure`.
The inherited reference-layout tests are applicable only to fast_dedup.
tests/test_variant.py provides independent affine-score and property tests. Inapplicable upstream reference tests are disabled for modes 2/3.

Ordinary CLI: `./build/msakit -i examples/tiny.fasta -o -`.
For website submission use benchmark/run.sh; it preserves complete titles and
IUPAC ambiguity codes. All three retain quadratic traceback memory in pair length;
large genomes can still exceed website time/memory limits.

See README.zh-CN.md and benchmark/README.zh-CN.md for the original framework guide.
