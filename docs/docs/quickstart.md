#Quickstart

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
