# Flanker

Flanker is a tool for studying the homology of gene-flanking sequences. It will annotate FASTA or Multi-FASTA files for specified genes, then write the flanking sequences to new FASTA files.

# Usage

| Required arguments  | Description |
| --- | --- |
| ```--fasta_file``` | Input .fasta file |
| ```--gene``` **OR** ```--list_of_genes``` | Single gene to be found **OR** line-seperated list of genes|

| Optional arguments | Description | Default|
| --- | --- | --- |
| ```--help``` | Displays help information then closes Flanker | ```False``` |
| ```--flank``` | Choose which side(s) of the gene to extract (left/right/both)| ```both``` |
| ```--mode``` | One of "default" - normal mode, "mm" - multi-allelic cluster, or "sm" - salami-mode| ```default``` |
| ```--circ``` | Add if your sequence is circularised | ```False``` |
| ```--include_gene``` | Add if you want the gene included in the output .fasta | ```False``` |
| ```--database``` | Specify the database Abricate will use to find the gene(s) | ```resfinder``` |
| ```--verbose``` | Increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity. | ```0``` |

| Window options | Description | Default |
| --- | --- | --- |
| ```--window``` | Flank length on either side of gene | ```1000``` |
| ```--wstop``` **AND** ```--wstep``` | For iterating: terminal flank length **AND** step size, ```--window``` becomes initial flank length | ```None``` |

| Clustering options | Description | Default |
| --- | --- | --- |
| ```--cluster``` | Use clustering mode? | ```False``` |
| ```--outfile``` | Prefix for clustering output file | - |
| ```--threshold``` | Mash distance threshold for clustering | 0.01 |

# Dependencies

The following Python packages are required for ```flanker.py```:

- [argparse](https://docs.python.org/3/library/argparse.html)
- [biopython](https://biopython.org)
- [numpy](https://numpy.org)
- [pandas](https://pandas.pydata.org)
- [pathlib](https://docs.python.org/3/library/pathlib.html)
- [subprocess](https://docs.python.org/3/library/subprocess.html)
- [sys](https://docs.python.org/3/library/sys.html)

- [mash](https://github.com/marbl/Mash) - if you want to do clustering
- [abricate](https://github.com/tseemann/abricate)


