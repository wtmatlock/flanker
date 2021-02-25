# Flanker

[Flanker](https://github.com/wtmatlock/flanker) is a tool for studying the homology of gene-flanking sequences. It will annotate FASTA/multi-FASTA files for specified genes, then write the flanking sequences to new FASTA files. There is also an optional step to cluster the flanks by sequence identity.


## Citation

If you use Flanker in your work please cite our [pre-print](https://www.biorxiv.org/content/10.1101/2021.02.22.432255v1):

Matlock W, Lipworth S, Constantinides B, Peto TEA, Walker AS, Crook D, Hopkins S, Shaw LP, Stoesser N. Flanker: a tool for comparative genomics of gene flanking regions. BioRxiv. 2021. doi: https://doi.org/10.1101/2021.02.22.432255


## Installation

### Bioconda

```
conda create -n flanker -c bioconda flanker
conda activate flanker
```
N.B. Bioconda can be tempermental for new releases, so please try the below methods if this doesn't work.

### Conda + pip

```
conda create -n flanker -c bioconda python=3 abricate=1.0.1 mash
conda activate flanker
pip install git+https://github.com/wtmatlock/flanker
```

### pip

```
pip install git+https://github.com/wtmatlock/flanker  # Install Mash and Abricate seperately
```

External dependencies:
* [Abricate](https://github.com/tseemann/abricate)
* [Mash](https://github.com/marbl/Mash)



## Quickstart

First we download some hybrid assemblies of plasmid genomes from *David, Sophia, et al. "Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread among Klebsiella pneumoniae." Proceedings of the National Academy of Sciences 117.40 (2020): 25043-25054.*

There are 44 plasmid genomes of which 16 and 28 contain *bla*KPC-2 and *bla*KPC-3, respectively. You'll need to take all replicons of interest (in this case the 44 plasmids) and concatenate these into a multi-FASTA file:

```
cat *fsa > david_plasmids.fasta
```

You should then rename the FASTA headers so that they match the original files (if this is not already the case). We have provided a simple script to do this (in /scripts/ directory on github):

```
ls *fsa | sed 's/[.]fsa//' > input_files
python multi_fa_rename.py david_plasmids.fasta input_files david_plasmids_renamed.fasta
mv david_plasmids_renamed.fasta david_plasmids.fasta
```

Now you are ready to use Flanker. In this example we are going to compare the flanking sequences around *bla*KPC-2. We are going to extract windows from 0bp (```--window```) to 5000bp (```--wstop```) base pairs in 100bp chuncks (```--wstep```) to the upsteam (```--flank upstream```) of the gene. We will include the gene (```--include_gene```) and use the default NCBI database.

```
flanker --flank upstream --window 0 --wstop 5000 --wstep 100 --gene blaKPC-2 --fasta_file david_plasmids.fasta --include_gene
```

You should now see many fasta files in the working directory containing upstream flanking regions from 0 to 4900 bp.



## Usage

| Required arguments  | Description |
| --- | --- |
| ```--fasta_file``` | Input .fasta file |
| ```--gene``` **OR** ```--list_of_genes``` | Space-delimited list of genes in the command line **OR** newline-delimited .txt of genes |


| Optional arguments | Description | Default|
| --- | --- | --- |
| ```--help``` | Displays help information then closes Flanker | ```False``` |
| ```--flank``` | Choose which side(s) of the gene to extract (upstream/downstream/both)| ```both``` |
| ```--mode``` | One of "default" - normal mode, "mm" - multi-allelic cluster, or "sm" - salami-mode| ```default``` |
| ```--circ``` | Add if your sequence is circularised | ```False``` |
| ```--include_gene``` | Add if you want the gene included in the output .fasta | ```False``` |
| ```--closest_match``` | Add if you want to find the closest match to your query by string edit distance | ```False```| 
| ```--database``` | Specify the database Abricate will use to find the gene(s) | ```ncbi``` |
| ```--verbose``` | Increase verbosity: 0 := only warnings, 1 := info, 2 := everything. | ```0``` |

| Window options | Description | Default |
| --- | --- | --- |
| ```--window``` | Flank length on either side of gene | ```1000``` |
| ```--wstop``` **AND** ```--wstep``` | For iterating: terminal flank length **AND** step size, ```--window``` becomes initial flank length | ```None``` |

| Clustering options | Description | Default |
| --- | --- | --- |
| ```--cluster``` | Use clustering mode? | ```False``` |
| ```--outfile``` | Prefix for clustering output file | - |
| ```--threshold``` | Mash distance threshold for clustering | ```0.001``` |
| ```--kmer_length```| *k*-mer length for Mash|```21```|
| ```--sketch_size```| Sketch size for Mash|```1000```|
| ```--threads```| Threads to use for Mash|```1```|

**N.B.** Gene queries use exact matching by default, so e.g. querying only ```bla``` will return nothing. Also be mindful that non-default databases, such as Resfinder, add indexing after annotation names e.g. ```blaCTX-M-15``` becomes ```blaCTX-M-15_1```. Please check your Abricate output if you are unsure of the naming conventions. In closest match mode, the annotation with the smallest [Levenshtein distance](https://en.wikipedia.org/wiki/Levenshtein_distance) to your query will be used.

Additionally, there should be **no whitespace** in your .fasta header names. This can be resolved using our [helper script](https://github.com/wtmatlock/flanker/blob/main/scripts/multi_fa_rename.py).

## Clustering

Having extracted flanking sequences around a gene, you might then want to cluster them into groups which share high sequence identity. Flanker does this using [single-linkage clustering](https://en.wikipedia.org/wiki/Single-linkage_clustering) of Mash distances. The method is very similar to that used by Ryan Wick in his [Assembly-Dereplicator](https://github.com/rrwick/Assembly-Dereplicator) package (and indeed we re-use several of his functions).

The ```-cl``` flag turns on clustering mode. The most important parameter is ```-tr``` - after running Mash, Flanker filters the ```Mash dist``` output so that only pairs with a mash dist smaller than this threshold are used to build the graph. You can also optionally tune the kmer-length ```-k``` and sketch size ```-s``` which mash uses, see [Mash-Docs](https://mash.readthedocs.io/en/latest/).

A seperate clustering file is produced for each window examined. The output is a comma separated file with two columns: sequence and cluster group. These can easily be combined for further processing:

```
cat out* | sed '/assembly/d'  > all_out
sed -i '1 i\flank,cluster' all_out
```

You can take this output and create figures similar to those in our manuscript (see the Binder on the Flanker GitHub page) or use in custom downstream applications.



## Multi-allelic mode

If you feed flanker a list of genes (```--list_of_genes```) in default mode (```--mode default```), flanker considers each of these in turn. However, if you turn on multi-allelic mode (```--mode mm```), it considers all genes in the list for each window. This allows you to detect flanking regions which are similar between different alleles of genes (e.g. *bla*KPC-2/3 etc) and between completely different genes.



## Salami mode

Salami mode (```--mode sm```)considers each window (of length ```--wstep```) from 0 to ```--wstop``` as a seperate contiguous sequence; in default mode these are concatenated together. This is intended to allow detection of recombination/mobile genetic elements which are occur in diverse genetic contexts. N.b at present Salami mode always starts at the edge of the gene and does not include the gene, therefore setting ```--w``` and ```-inc``` currently do nothing in this mode. As with the default mode, you can optionally perform clustering on the output from Salami Mode by using the ```-cl``` flag.

For instance, here we extract 100bp windows from 0-5000bp downstream of the *bla*TEM-1B gene.

```
flanker --fasta_file example.fasta --gene blaTEM-1B --wstop 5000 --wstep 100 --flank downstream --mode sm  
```

## Contributing

If you would like to contribute to this project, please [open an issue]((https://github.com/wtmatlock/flanker/issues)) or contact the authors directly.

To install and test a development build:

```
git clone https://github.com/wtmatlock/flanker
cd flanker
conda create -n flanker -c bioconda python=3 abricate=1.0.1 mash pytest
conda activate flanker
pip install --editable ./
pytest
```

[![Documentation Status](https://readthedocs.org/projects/flanker/badge/?version=latest)](https://flanker.readthedocs.io/en/latest/?badge=latest)
