# genome-dereplicator

Genome Dereplicator (gender) is a single-file Python script to dereplicate large numbers of genome assemblies and/or metagenomically-assembled genomes. The tool identifies groups of highly similar genomes based on the average nucleotide identity (ANI) threshold and reduces each group to a single representative genome. In other words, `gender` takes as an input a set of genome assemblies and removes genomes for which there are sufficiently close relatives (above sequence identity threshold), resulting in a smaller set where the assemblies are more unique.

I wrote `gender` inspired by the [Assembly Dereplicator](https://github.com/rrwick/Assembly-Dereplicator) tool developed by [Ryan Wick](https://github.com/rrwick). Compared to Assembly Dereplicator, `gender` runs faster, has lower RAM usage, provides information on genome clusters, and takes into account all *k*-mers in dereplication.

## Requirements

* Python >= 3.8
* [Kmer-db](https://github.com/refresh-bio/kmer-db) >= 1.10.0
   > Check with `kmer-db -h`
* [SeqKit](https://pandas.pydata.org/) >= 2.2.0 
   > Check with: `seqkit -h`


## Installation
No installation is required. You can simply clone it and run it:

```
git clone https://github.com/aziele/genome-dereplicator
cd genome-dereplicator
./gender.py --help
```

## Quick Start
`gender.py` takes as input a directory containing your genome assemblies in FASTA format, an output directory, and minimum sequence identity.

```
./gender.py --ani 0.95 example/indir example/outdir
```

## Method

### Average Nucleotide Identity (ANI) 

ANI is the most widely used measure of global sequence similarity between two genomes. The measure has been the gold standard whole-genome method for prokaryotic and virus species identification. For example, a pair of genomes with ANI of greater than 95% is considered of the same species. `gender` estimates ANI based on similarity/distance measures calculated by Kmer-db. By default, the estimation is based on Mash distance as it is well correlated to approximately `1 - ANI`, such that a Mash distance of `0.05` corresponds to an ANI of `0.95` (95% sequence identity between genomes). Alternatively, `gender` can use the `min` and `max` similarity measures to estimate ANI, which are defined as the intersection of *k*-mers between two genomes divided by the length of the smallest (min) or largest (max) of the two genomes.

### Clustering

`gender` clusters genome assemblies and chooses a single representative genome per cluster. First, the script uses Kmer-db to calculate ANI between all-vs-all genomes. It then builds an undirected graph where the genomes are vertices, and any two genomes above ANI are connected by an edge. Next, it identifies genome clusters using single-linkage clustering (it enables a genome to join a cluster as soon as any genome if that cluster is within the specified ANI threshold). Finally, it chooses a single representative genome in each cluster based on the largest N50 (in this way complete genomes are preferred over draft genomes). 

### Iterative clustering (batches)
`gender` uses the Assembly Dereplicator's method that allows you to control how many genomes will be clustered at once. Briefly, the tool takes as input the path to a directory containing the genomes to be dereplicated, rearranges them randomly, and separates them into smaller batches (50,000 genomes per batch by default). The clustering of each batch yields clusters with single representatives, and all non-representative genomes are removed from the next rounds of clustering. This procedure is repeated until an iteration fails to remove at least one genome from the current batch. Finally, the tool conducts a final dereplication round with all remaining representative genomes. 

Of note, although using batches reduces the RAM usage, it may also slightly change the number of representative genomes (final clusters). The smaller the batch size, the greater the chance you will get slightly more representative genomes. Depending on your purpose of dereplication, this could be advantageous. As stated in the Assembly Dereplicator's documentation, batches can help a bit with a common criticism of single-linkage clustering: long thin clusters (where nearby genomes of the same cluster have high ANI, but genomes at opposite ends of a cluster may be much farther from each other below ANI threshold). Clustering subsets at a time increases the chance that large clusters will be broken apart.


## Usage

```
./gender.py [options] <in_dir> <out_dir>
```

Positional arguments:

- `in_dir` - input directory with genome FASTA files (gzipped or not),
- `out_dir` - output directory with representative genome FASTA files and clustering information

Options:

- `--ani` - average nucleotide identity (ANI) above which two genomes will cluster together. With the default `0.95` threshold, genomes closer than ~95% identity will cluster together. [Default: `0.95`]
- `--batch_size` - controls how many genomes will be clustered at once (see Assembly Dereplicator](https://github.com/rrwick/Assembly-Dereplicator#batches). [Default: `50000`] 
- `--k` - k-mer size [Default: `21`]
- `--k_fraction` - a fraction (from 0 to 1) of *k*-mers that are used to estimate ANI [Default: `1`].
- `--distance` - a distance method to estimate ANI [Default: `mash`]
- `--threads` - number of CPU threads [Default: number of cores]

## Tests
For testing purposes, `gender` comes with automated tests and a toy data set (10 genome files of 100 kbp each with different ANI - e.g., 100%, 99%, 98%, and 95%). If you want to check that everything works as intended, just run:

```
./test.py
```

## Benchmark
I measured the running time and peak RAM memory of `gender` and Assembly Dereplicator (AD) using 16 threads on a Linux PC with Intel(R) Xeon(R) W-2295 3.00GHz processor.

### Phage genomes (n = 27,608)

<table>
<thead>
	<tr>
		<th>Tool</th>
		<th>ANI</th>
		<th><i>k</i>-mer fraction</th>
		<th>Sketch size</th>
		<th>Clusters</th>
		<th>Time [HH:MM:SS]</th>
		<th>RAM peak</th>
	</tr>
</thead>
<tbody>
	<tr>
		<td>gender</td>
		<td>95%</td>
		<td>0.1</td>
		<td>~10,000</td>
		<td>8949</td>
		<td>00:00:49</td>
		<td>2 GB</td>
	</tr>
   <tr>
		<td>AD</td>
		<td>95%</td>
		<td>~0.1</td>
		<td>10,000</td>
		<td>8945</td>
		<td>02:03:10</td>
		<td>233 GB</td>
	</tr>
	<tr>
		<td>gender</td>
		<td>95%</td>
		<td>1</td>
		<td>~100,000</td>
		<td>8942</td>
		<td>00:01:33</td>
		<td>7 GB</td>
	</tr>
   <tr>
		<td>AD</td>
		<td>95%</td>
		<td>~1</td>
		<td>100,000</td>
		<td>8942</td>
		<td>13:14:43</td>
		<td>233 GB</td>
	</tr>
</tbody>
</table>

### *Escherichia coli* genomes (n = 10,000)

<table>
<thead>
	<tr>
		<th>Tool</th>
		<th>ANI</th>
		<th><i>k</i>-mer fraction</th>
		<th>Sketch size</th>
		<th>Clusters</th>
		<th>Time [HH:MM:SS]</th>
		<th>RAM peak</th>
	</tr>
</thead>
<tbody>
	<tr>
		<td>gender</td>
		<td>99%</td>
		<td>0.002</td>
		<td>~10,000</td>
		<td>213</td>
		<td>00:01:03</td>
		<td>3 GB</td>
	</tr>
   <tr>
		<td>AD</td>
		<td>99%</td>
		<td>~0.002</td>
		<td>10,000</td>
		<td>203</td>
		<td>00:25:53</td>
		<td>36 GB</td>
	</tr>
	<tr>
		<td>gender</td>
		<td>99%</td>
		<td>0.02</td>
		<td>~100,000</td>
		<td>215</td>
		<td>00:02:32</td>
		<td>3 GB</td>
	</tr>
   <tr>
		<td>AD</td>
		<td>95%</td>
		<td>~0.02</td>
		<td>100,000</td>
		<td>209</td>
		<td>3:08:53</td>
		<td>37 GB</td>
	</tr>
	<tr>
		<td>gender</td>
		<td>99%</td>
		<td>1</td>
		<td>~5,000,000</td>
		<td>211</td>
		<td>00:55:49</td>
		<td>9.7 GB</td>
	</tr>
</tbody>
</table>
