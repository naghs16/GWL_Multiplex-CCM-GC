# GWL_Multiplex-CCM-GC
The source codes are provided in this repository used to preprocess the data and reproduce the results from the statistical approach in Naghipour et al. (2023).

The codes are performed in a R environment by several packages.

If you use the codes, please cite our recent work in your publications:

Naghipour L., Kamran K.V., Ghanbari A., Nourani V., Aalami M.T. and Zhang Y. (2023) Disentangling groundwater level fluctuations from the analysis of spurious connectivity in multiplex networksl. Under Review in Environmental Modelling and Software 

For any questions/suggestions about the codes/approach, please contact naghipour [dot] l [at] tabrizu [dot] ac [dot] ir.

# Computer Requirements
Intel® Pentium® CPU G2030 @ 3.00GHz 3GHz with 2 GByte RAM and 128 GByte Hardware (the analysis used around 5 GByte of the Hardware)

# Setup Instructions
If you already have a R environment, you could ignore this part and continue with explanations of the next part. Otherwise, you should install R based on your operating system (i.e. Windows etc). Microsoft R Open is the enhanced distribution of R from Microsoft Corporation, and this distribution is highly recommended by the scientists.  

## Prerequisites
In the folowing, all instructions are performed in R 4.0.2 with Ubuntu desktop operating system (as well as Windows 10). 

The source codes need several packages, including rEDM, igraph, brainGraph.

```{r}
# install the packages
install.packages(c('rEDM','igraph','brainGraph'))
```

Note that there is a regular update for the package rEDM, while the source code CCM.R is operated based on an old version of this package. Please follow the instruction released by the Sugihara Lab to be sure that an implementation of 'EDM' algorithms is based on the latest version. 

# Run the Codes
The modeling procedure consists of five steps, including I, II, III, IV, V.

I) At first, groundwater level measurements are preprocessing by running four scripts as,

```{r}
Rscript Data_Generate.R
Rscript Preprocessing.Rmd
Rscript generate_data_modeling.Rmd
Rscript Statistics.R
```

II) In the second step, the CCM and GC methods are used to infer causal relationships for all possible pairs of interactions among the nodes by running Rscript CCM.R and GC.R with one argument, namely layer.type (FOR, INV) as,

```{r}
Rscript CCM.R FOR
Rscript GC.R FOR
```

III) Then, the values of z-score are obtained according to the surrogates data-sets for the MON and MUX dynamics by Rscript ZScore_MON.R and ZScore_MUX.R, respectively. 

```{r}
Rscript ZScore_MON.R 
Rscript ZScore_MUX.R 
```

IV) In the next step, the networks are constructed from Rscript Graph.R using the adjacency matrices obtained by threshording the z-score between each pair of the nodes.

```{r}
Rscript Graph.R
```

V) Finally, the connections of the networks are characterized in local and global scales by Local-Measures.R and Global-Measures.R, respectively. Moreover, directionality is inferred according to Directionality.R as,

```{r}
Rscript Local-Measures.R 
Rscript Global-Measures.R 
Rscript Directionality.R
```
