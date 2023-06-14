## kdpp
k-DPP for subset selection using VLMC distances

This repository contains an example on how to use DPPy and Tensorflow to sample exactly from k-DPP and evaluate log probability.
We use Nextflow to implement a pipeline to compute VLMC pairwise distances between sequencing reads.
To use the Nextflow pipeline, you have to install the VLMC library [dvstar](https://github.com/Schlieplab/dvstar).

Input: Pairwise distances between VLMC

Output: Samples from k-DPP

#### Requirement
* Nextflow
* BBTools such as bbduk.sh and reformat.sh
* [dvstar](https://github.com/Schlieplab/dvstar) to generate pairwise distances between VLMC 
* [DPPy](https://github.com/guilgautier/DPPy)
* [Tensorflow](https://www.tensorflow.org/install)


#### TODO
Improve conversion from distance to similarity and implement Gibbs sampling.

#### Acknowledgement

Alexander Schliep, Joel Gustafsson for help with [dvstar](https://github.com/Schlieplab/dvstar)
