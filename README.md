OutFLANK v 0.2
========


A procedure to find Fst outliers based on an inferred distribution of neutral Fst


Version 0.2 was released in November 2017. This version fixes a bug in the OutFLANK algorithm that failed to correctly remove low heterozygosity loci when calculating neutral mean FST and df.


First, read the description of how to use OutFLANK in "OutFLANK readme.pdf" available at github.com/whitlock/OutFLANK.


After you install the package and call the library, you can view the vignette.


The vignette illustrates best practices and data checks when using the package.


The vignette can be installed by typing:
`browseVignettes("OutFLANK")` in R


This will open a web browser, and click on the `html` link.


You can also view the vignette on the web at https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html

# Input and interpretation safeguards

Before running OutFLANK, verify the following points:

- OutFLANK is designed for diploid, genotype-derived allele-count input. Haploid data are not directly supported; do not silently duplicate or recode haploid calls as diploid genotypes.
- Population labels and sample counts must align exactly across the input matrix and population metadata. Check spelling, ordering, and the NumberOfSamples value.
- Missing genotypes must use one documented encoding consistently. Do not convert unknown or malformed calls into reference or alternate homozygotes.
- Preserve locus and sample identifiers during format conversion and inspect allele counts before analysis.
- Large SNP datasets can be memory-intensive; test a representative subset, monitor memory, and record filtering and thinning choices.
- Report OutFLANK, R, dependency versions, input encoding, thresholds, and all function arguments.

Interpret FST, corrected FST, heterozygosity, degrees of freedom, p-values, q-values, and outlier flags together; an outlier flag alone is not proof of local adaptation.
