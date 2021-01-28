# Salami mode

Caveat: this is still under development and should be treated as experimental, use at your own risk. Please report any bugs.

Salami mode considers each window (of length ```-wstep```) from ```-w``` to ```-wstop``` as a seperate entity; in default mode these are concatenated together. This is intended to allow detection of recombination/mobile genetic elements which are occur in diverse genetic contexts.

Example:
```
  python flanker.py -i example.fasta  -g blaTEM-1B_1 -w 0 -w 4900 -f left  -m SM  
```

Here we extract 100bp windows from 0-4900 bp to the left of the blaTEM-1B gene.

![Schematic](/images/flanker_schematic.jpg)
*Figure 1 - Illustrating the difference between default and salami mode*
