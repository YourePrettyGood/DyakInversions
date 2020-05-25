# DyakInversions
Scripts used for analyses of inversions in Drosophila yakuba

## Diversity and divergence analyses:

The `pi` directory contains results from the diversity (Tajima's pi) and divergence (Nei's Dxy and the Hudson, Slatkin, and Maddison F_ST estimator) analyses.

The analyses within Dyak took a relatively standard approach, averaging only sites where less than 50% of individuals had missing data. Windowed F_ST was calculated as a ratio of averages (i.e. calculate the windowed averages of expected heterozygosity and Dxy first, then take the ratio of these averages as the final F_ST for the window. This (and the use of the Hudson, Slatkin, and Maddison (1992) Genetics F_ST estimator) is based on recommendations from Bhatia et al. (2013) Genome Research.

Analyses involving Dsan required projection of alleles from the *D. santomea* STO CAGO 1482 assembly coordinate space into *D. yakuba* NY73PB assembly coordinate space, since all Dyak pseudoreferences are in NY73PB coordinate space. These projections were performed using `projectPseudorefOnMAF.pl` using a MAF of 1-to-1 regions derived from pairwise alignment of the two assemblies using [LAST](http://last.cbrc.jp).

## Approximate Bayesian Computation (ABC) for inversion age estimation:

The ABC inversion age estimates can be found in the `pi/Dsan_projections/Dyak_ABC` directory. I haven't included the summary statistic values from the simulations from the prior distributions, as these files exceed the storage limit of Github repositories. However, the ABCreg draws from the posterior and the observed values of the summary stats are provided, along with the R script and resulting plot.

**TODO: I need to add the Python scripts for prior simulation to this repository.**

## STRUCTURE-based inversion karyotyping:

The inversion karyotyping analyses are in the `MavericK` directory.

Inversion karyotyping is performed across a group of samples by identifying biallelic SNPs across the samples, splitting into windows, and estimating the number of "ancestral populations" (K) and their contributions to the extant samples (the Q matrix of Pritchard, Stephens, and Donnelly (2000) Genetics). I simultaneously estimate K and the Q matrix for each window using [MavericK](http://bobverity.com/home/maverick/what-is-maverick/) version 1.0.4.

The Q matrices are aligned across windows based on a focal set of samples with known karyotypes. The algorithm for this step is very naive, and doesn't really enable detection of novel inversions, so any suggestions to improve it would be very welcome.

I've tried quite a few renditions of this pipeline (e.g. originally with STRUCTURE, then Structurama, now MavericK, and more recent attempts to get ALStructure working with it). There are other published pipelines performing similar but not identical tasks. For example, the [lostruct](https://github.com/petrelharp/local_pca) package from Peter Ralph (which doesn't seem to be able to karyotype easily, but definitely detects the boundaries of inversions well), and [asaph](https://github.com/rnowling/asaph) from Nowling, Manke, and Emrich (2019) BioRxiv (again, works for detection well, but I'm not sure if it works for karyotyping).

The pipeline is pretty darn slow even for tens of Drosophila samples (G=150 Mbp). ALStructure was very appealing in terms of performance and the potential for more accurate inference of Q matrices, but the K estimator used (derived from Leek (2011) BioMetrics, ish) seems to consistently overestimate K for 100 kb windows. ALStructure also errors out when the K estimate is less than 2 (i.e. the majority of windows), but that is easy enough to handle.

Again, any suggestions for alternative approaches would be greatly appreciated! These Dyak inversions have very obvious signals (even just from heterozygosity and divergence from reference), so they're quite useful for diagnostics and troubleshooting.

## Differential expression between inversions and ASE analysis:

RNA-seq data for differential expression analyses are derived from Rogers et al. (2014) G3 and Rogers, Shao, and Thornton (2017) PLoS Genetics. These reads were selectively-aligned (Srivastava et al. (2019) BioRxiv) using [salmon](https://github.com/COMBINE-lab/salmon) commit df3e6ab against the extracted transcriptome of *D. yakuba* NY73PB (plus the genome as decoy, as per the selective-alignment algorithm) with 100 bootstraps per library.

Transcript quantification results were loaded into R version 3.5.2 using tximport version 1.15.5, and differential expression was assessed using DESeq2 version 1.22.2, with log fold changes shrunken using the empirical Bayes approach with t-distributed prior as per Zhu, Ibrahim, and Love (2019) Bioinformatics (apeglm version 1.4.2). Significance was assessed using the Wald test, and FDR controlled at alpha=0.05 using the Benjamini and Hochberg (1995) procedure.

Transcript quantification for the allele-specific expression (ASE) analysis was performed using the non-competitive Allele-Specific Expression [ncASE](https://github.com/YourePrettyGood/ncASE) pipeline from Gutin (2019) Princeton University Dissertations. Briefly, reads were mapped to the extracted *D. yakuba* NY73PB transcriptome and a *D. yakuba* Tai18E2 transcriptome in NY73PB coordinate space (i.e. extracted from a Tai18E2 pseudoreference). Variants were called on both transcriptomes using BCFtools, and allele-specific read counts were extracted from each VCF using a custom awk script `vcf_to_ncASE.awk`. In order to minimize the effects of reference bias, the ncASE workflow filters for biallelic SNPs where alleles match between variant calls on each transcriptome, and allele depth on one parental transcriptome is within a threshold deviance of allele depth on the other parental transcriptome. For our analyses, we chose a threshold of 10% deviance. Counts for each transcript are then summarized as the weighted sum of per-SNP allele depths, with the weight of SNP i defined as w_{i} = \(\sum_{j \in S} 1 - \frac{\|p_{i} - p_{j}\|}{L}\)^{-1}

where S is the set of SNPs in the transcript, p_{i} is the position of SNP i, and L is the average read length. The resulting transcript count matrices were rounded to the nearest integer, and imported into R using DESeq2's functions for count matrix import. As for the DE analyses, we shrunk log fold changes with apeglm, assessed significance with the Wald test, and controlled FDR at alpha=0.05 using the Benjamini and Hochberg (1995) procedure.

GO enrichment analyses were performed using GOseq (Young et al. (2014) Genome Biology) based on functional annotations generated with Trinotate.

**TODO: Add the R script used for these DE and ASE analyses to this repo.**
