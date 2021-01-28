# Multi-allelic mode

If you feed flanker a list of genes (```-lg```) in default mode (```-m Default```), flanker considers each of these in turn. If you turn on multi-allelic mode (-m ```mm```) however, it considers all genes in the list for each window. This allows you to detect flanking regions which are similar between different alleles of genes (e.g. blaCPC-2/3 etc) and between completely different genes. 
