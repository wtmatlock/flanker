# About Flanker

Flanker is a tool for studying the homology of gene-flanking sequences. It will annotate FASTA or Multi-FASTA files for specified genes, then write the flanking sequences to new FASTA files.

It allows a relatively agnostic approach to the study of gene flanking sequences because there is no reliance on underlying reference databases/sequences.

## Usage

| Required arguments  | Description |
| --- | --- |
| ```--fasta_file``` | Input .fasta file |
| ```--gene```| Space-delimited list of genes to annotate |

<br/><br/>


| Optional arguments | Description | Default|
| --- | --- | --- |
| ```--help``` | Displays help information then closes Flanker | ```False``` |
| ```--flank``` | Choose which side(s) of the gene to extract (left/right/both)| ```both``` |
| ```--mode``` | One of "default" - normal mode, "mm" - multi-allelic cluster, or "sm" - salami-mode| ```default``` |
| ```--circ``` | Add if your sequence is circularised | ```False``` |
| ```--include_gene``` | Add if you want the gene included in the output .fasta | ```False``` |
| ```--database``` | Specify the database Abricate will use to find the gene(s). There are several shipped with abricate but it is also possible (and easy) to configure your own | ```resfinder``` |
| ```--verbose``` | Increase verbosity: 0 = only warnings, 1 = info, 2 = debug. No number means info. Default is no verbosity. | ```0``` |

<br/><br/>


| Window options | Description | Default |
| --- | --- | --- |
| ```--window``` | Flank length on either side of gene | ```1000``` |
| ```--wstop``` **AND** ```--wstep``` | For iterating: terminal flank length **AND** step size, ```--window``` becomes initial flank length | ```None``` |

<br/><br/>


| Clustering options | Description | Default |
| --- | --- | --- |
| ```--cluster``` | Use clustering mode? | ```False``` |
| ```--outfile``` | Prefix for clustering output file | - |
| ```--threshold``` | Mash distance threshold for clustering | 0.01 |
| ```--threads``` | Threads to use for mash (makes little difference unless dealing with huge volumes of data) | 1 |
