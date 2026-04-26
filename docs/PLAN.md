# Project Plan -- De Bruijn Graph Genome Assembly

## Goal
Build a from scratch library (`debruijn`) that implements a minimal but correct Velvet style short read genome assembler. The library covers k-mer counting, de Bruijn graph construction, tip removal, bubble popping, Hierholzer's Eulerian path algorithm, and maximal non branching path extraction for contig output. Validate on simulated reads from the lambda phage genome and benchmark against SPAdes on a small bacterial or phage dataset.

## Layers

### Layer 0 -- Foundations
Read and write FASTA and FASTQ. Encode DNA as integer arrays. Handle reverse complements so that both strands contribute to the graph.

### Layer 1 -- K-mer operations
Canonical k-mer encoding (lexicographically smaller of k-mer and its reverse complement). Exact k-mer counting via a Python dict with numpy bucketing. Approximate k-mer counting via count min sketch for memory constrained settings.

### Layer 2 -- De Bruijn graph
Node set indexed by distinct (k-1)-mers. Edge set indexed by distinct k-mers. Edge weights = observed read coverage. Expose adjacency structure with in degree and out degree lookups and a neighbors iterator.

### Layer 3 -- Graph cleanup
Tip removal: iteratively prune paths of length less than 2k that end at a degree 1 node. Bubble popping: identify pairs of parallel paths between two shared nodes; keep the higher coverage path and remove the other. Low coverage edge removal: delete edges with coverage below a threshold proportional to the mean.

### Layer 4 -- Eulerian traversal
Hierholzer's algorithm in linear time for finding an Eulerian circuit or path. Handles the case where the graph is a disjoint union of connected components by returning a set of circuits.

### Layer 5 -- Contig extraction
Maximal non branching path collapse: any internal node with in degree and out degree both equal to 1 is merged with its successor. The resulting edges are contigs. Output as FASTA. Compute N50, total assembly length, number of contigs.

### Layer 6 -- Visualization
Plot the simplified graph (nodes and edges) using matplotlib and networkx. Plot the k-mer coverage histogram (used to set the coverage cutoff). Plot N50 as a function of k to help choose the right k.

### Layer 7 -- External comparison
Optional bridge to SPAdes via subprocess. Run SPAdes on the same read set and compare the N50 and fraction of reference genome recovered side by side with our assembler.

## Data

Lambda phage genome (about 48 kb) and simulated short reads. The lambda genome is small enough to run quickly and has published assemblies to compare against. We also include a handful of synthetic test sequences with known structure to verify the algorithm on idealized input.

## Success criteria

On error free simulated reads of the lambda phage genome with k=21 and 30x coverage, the assembler should produce a single contig equal to the reference genome (up to starting point for a circular genome). With 1 percent simulated error rate, after tip removal and bubble popping, the largest contig should cover more than 90 percent of the reference. On a small real phage read set, SPAdes should beat our assembler (it is a production tool) but our N50 should be within 2x of SPAdes.

## Testing

Brute force tests on toy sequences where the graph structure is fully known. Agreement between the de Bruijn Eulerian assembly and the known reference on error free short inputs. Idempotence of tip removal and bubble popping (running cleanup twice gives the same graph). N50 increases monotonically with k below the repeat threshold.
