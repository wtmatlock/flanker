# Clustering

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
