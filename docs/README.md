# Flanker

Flanker is a tool for studying the homology of gene-flanking sequences. It will annotate FASTA or Multi-FASTA files for specified genes, then write the flanking sequences to new FASTA files.

## Installation

### Conda

This is the recommended method and by far the easiest. It will be available soon.

### During development

```
conda create -n flanker python=3 abricate biopython mash  # needs bioconda channel
pip install --editable git+https://github.com/wtmatlock/flanker
pip install --editable /path/to/flanker/package  # where 'package' contains setup.py
```

### Run tests

```
pytest
```

### From GitHub

Ensure you have all dependencies installed:

**Python dependencies:**

* multiprocessing
* sys
* argparse
* pandas
* numpy
* subprocess
* biopython
* pathlib
* time
* logging
* tempfile
* os
* collections
* glob
* networkx

**External software:**

* abricate
* mash

Then simply clone the repository:

```
  git clone https://github.com/wtmatlock/flanker
```

and check everything is working:

```
  python flanker.py --help
```

## Quickstart

First we download some hybrid assemblies of plasmid genomes from *David, Sophia, et al. "Integrated chromosomal and plasmid sequence analyses reveal diverse modes of carbapenemase gene spread among Klebsiella pneumoniae." Proceedings of the National Academy of Sciences 117.40 (2020): 25043-25054.*

There are 44 plasmid genomes of which 16 contain blaKPC-2 and 28 blaKPC-3.

You need to take all replicons of interest (in this case the 44 plasmids) and concatenate these into a multifasta file.

```
  cat *fsa > david_plasmids.fasta
```

You should then rename the fasta headers so that they match the original files. We have provided a simple script to do this.

e.g.

```
  ls *fsa | sed 's/[.]fsa//' > input_files
  python multi_fa_rename.py david_plasmids.fasta input_files david_plasmids_renamed.fasta
  mv david_plasmids_renamed.fasta david_plasmids.fasta
```

Now you are ready to use Flanker. In this example we are going to compare the flanking sequences around blaKPC-2. We are going to extract windows from 0 (```-w```) to 5000 (```-wstop```) base pairs in 100bp chuncks (```-wstep```) to the left (```-f left```) (downstream) of the gene. We will include the gene (```-inc```) and use the default resfinder database.

```
  python flanker.py -f left -w 0 -wstop 5000 -wstep 100 -p 8 -v 1 -g blaKPC-2_1 -i david_plasmids.fasta -inc
```

You should now see many fasta files in the working directory containing left flanking regions from 0 to 4900 bp.

## Usage

| Required arguments  | Description |
| --- | --- |
| ```--fasta_file``` | Input .fasta file |
| ```--gene```| Space-delimited list of genes to annotate |

| Optional arguments | Description | Default|
| --- | --- | --- |
| ```--help``` | Displays help information then closes Flanker | ```False``` |
| ```--flank``` | Choose which side(s) of the gene to extract (upstream/downstream/both)| ```both``` |
| ```--mode``` | One of "default" - normal mode, "mm" - multi-allelic cluster, or "sm" - salami-mode| ```default``` |
| ```--circ``` | Add if your sequence is circularised | ```False``` |
| ```--include_gene``` | Add if you want the gene included in the output .fasta | ```False``` |
| ```--database``` | Specify the database Abricate will use to find the gene(s) | ```resfinder``` |
| ```--verbose``` | Increase verbosity: 0 := only warnings, 1 := info, 2 := everything. | ```0``` |

| Window options | Description | Default |
| --- | --- | --- |
| ```--window``` | Flank length on either side of gene | ```1000``` |
| ```--wstop``` **AND** ```--wstep``` | For iterating: terminal flank length **AND** step size, ```--window``` becomes initial flank length | ```None``` |

| Clustering options | Description | Default |
| --- | --- | --- |
| ```--cluster``` | Use clustering mode? | ```False``` |
| ```--outfile``` | Prefix for clustering output file | - |
| ```--threshold``` | Mash distance threshold for clustering | ```0.01``` |

## Clustering

Having extracted flanking sequences around a gene you might then want to cluster them into groups which share high sequence identity. Flanker does this using single linkage clustering of Mash distances. The method is very similar to that used by Ryan Wick in his [Assembly-Dereplicator](https://github.com/rrwick/Assembly-Dereplicator) package (and indeed we borrow several of his functions).


| Option    | Explanation   |
|-------- | -------- |
| ```-tr``` | mash clustering threshold, flanks are grouped into the same cluster if they are <= this |
| ```-p``` | threads for running mash (will make little difference unless very large input) |
| ```-o``` | prefix for the output clustering files |

<br/><br/>

A seperate clustering file is produced for each window examined. The output is a comma separated file with two colomns: isolate and clustering group.

These can easily be combined for further processing e.g.

```
  cat out* | sed '/assembly/d'  > all_out
  sed -i '1 i\flank,cluster' all_out
```

## Multi-allelic mode

If you feed flanker a list of genes (```-lg```) in default mode (```-m Default```), flanker considers each of these in turn. If you turn on multi-allelic mode (-m ```mm```) however, it considers all genes in the list for each window. This allows you to detect flanking regions which are similar between different alleles of genes (e.g. blaCPC-2/3 etc) and between completely different genes. 
