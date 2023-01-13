# AlgBio
Algorithmic bioinformatic tasks for MIPT coursework

## Task 1 - De Brujn Assembly

Based on the given .fastq files we have to give an output of assembled contigs.\
Folder for data and outputs - `./data` and `./contigs` respectively.\
Following packages are required:
* [networkx](https://networkx.org/)
* [pyfastx](https://github.com/lmdu/pyfastx)
* tqdm

## Task 2 - Mass spectrometry

Based on the given unknown spectrum we have to determin the proteins' amino acid chain:\
![](https://github.com/khomi-a/AlgBio/blob/main/2.%20Mass%20spec%20Tyrocidine/spectrum.png)\
Following packages are required:
* numpy

## Task 3 - Synteny search

We have to:
* find synteny blocks for X-chromosomes of mouse and human.
* find 2-break distances between circualized chromosomes.
* determine *maxDistance* and *minSize* influence on these distances.

![](https://github.com/khomi-a/AlgBio/blob/main/3.%20Synteny/anchors.jpg)\
Following packages are required:
* [Bio](https://pypi.org/project/bio/)