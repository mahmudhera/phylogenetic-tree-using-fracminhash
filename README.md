# Phylogenetic tree using FracMinHash

In this experiment, we will compute pairwise distance using FMH, and then use the distance matrix to compute Phylogenetic tree.

Steps:

1. Create folder for genomes
1. Create a file : list of genomes
1. Write code that will
      Open all genomes
      Then, construct fmh sketches for a bunch of different seeds (1 - 100)
1. After these, write code that will read genome list, seeds, and output pairwise distances using fmh
