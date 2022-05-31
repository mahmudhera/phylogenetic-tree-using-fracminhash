# Phylogenetic tree using FracMinHash

In this experiment, we will compute pairwise distance using FMH, and then use the distance matrix to compute Phylogenetic tree.

Steps:

1. Create folder for genomes
1. Create a file : list of genomes. Line1: name, line2: folder path
1. Write code that will
      Open all genomes
      Then, construct fmh sketches for a bunch of different seeds (1 - 100)
1. After these, write code that will read genome list, seeds, and output pairwise distances using fmh


Tree building method that takes into consideration uncertainties
Look into bootstrapping techniques


Robensen-Folds distance


Add more bacterial genomes
Start with a known Phylo tree based off of 16S
GTDB te start korbo??

copy, WABI format

check bootstrapping values. sample dist matrix a lot of times. each time see if tree branches change or not


PCOA -- compare the two dist matrices -- MDS -- multi dimensional scaling

take the organisms -- do MSA -- get pairwide evo dist
