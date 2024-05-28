# PanDA
PanDA is a joint discriminant analysis method aimed at fusing multi-omics datasets through finding a discriminant common latent space. PanDA captures cross-omics interaction and consistency and uses an uncorrelated constraint to ensure that the extracted latent components for each omics (omics-specific components) are not highly correlated. Because the components extracted using PanDA contain valuable discriminant information, we refer to them as discriminant components. These components can be used as inputs to several multi-omics analysis tools to enable efficient, improved downstream analysis. Here and in our paper, we demonstrated the advantages of PanDA over ten integrative multi-omics methods through four distinct downstream analyses: single-cell multi-omics data visualization, patient (or tumor) classification, biomarker identification, and clinical outcome prediction.

<p align="center">
  <img src="docs/panda results.png" width="900">
</p>

### Installation
PanDA can be installed by simply running the following code:
```
## Install development version
remotes::install_github("WuLabMDA/PANDA") 
```
or by running:
``` r
remotes::install_github("muhammadaminu47/PANDA")
```

## Input data

The input dataset for all the experiments reported in our paper can be found in the `./data/` folder of this repository. 
   
### Tutorials

1. [Basic PanDA example on simulated dataset](https://)
2. [PanDA identifies important markers related to breast cancer](https://)
3. [PanDA application to clinical outcome prediction](https://)

The benchmarking results on the TCGA multi-omics cancer datasets can be obtained by running the corresponding code in the `./code/` folder.

### Example work flow
An example of the `PanDA` work flow to get started. Here we used the processed TCGA multiomics dataset from the mixOmics package

```{r}
library(PANDA)
library(mixOmics) # import the mixOmics library
data(breast.TCGA) # extract the TCGA data

data = list()
data$mirna <- t(breast.TCGA$data.train$mirna)
data$mrna <- t(breast.TCGA$data.train$mrna)
data$protein <- t(breast.TCGA$data.train$protein)

Y <- breast.TCGA$data.train$subtype # use the subtype as the outcome variable
subtype <- factor(Y)
```

Extract discriminant latent components using PanDA.

```{r}
labels <- as.numeric(Y)
numComponents <- 10
PanDAModel <- PanDA(data,labels,numComponents)
```

Plot the discriminant latent representations for the different omics data

```{r}
mirnaComp <- as.data.frame(DOmicsModel[["DOmicsComponents"]][["mirnaComponents"]])
mrnaComp <- as.data.frame(DOmicsModel[["DOmicsComponents"]][["mrnaComponents"]])
proteinComp <- as.data.frame(DOmicsModel[["DOmicsComponents"]][["proteinComponents"]])

library(plotly)
cols = c('#BF382A', '#0C4B8E', "#fc8d59")
fig1 <- plot_ly(mrnaComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`, color = ~subtype, colors = cols)
fig1 <- fig1 %>% add_markers()
fig1 <- fig1 %>% layout(scene = list(xaxis = list(title = 'DC 1'),
                                     yaxis = list(title = 'DC 2'),
                                     zaxis = list(title = 'DC 3')))
fig1

```

<p align="center">
  <img src="docs/mrna breast comp.png" width="500">
</p>

```{r}
fig2 <- plot_ly(mirnaComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`, color = ~subtype, colors = cols)
fig2 <- fig2 %>% add_markers()
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = 'DC 1'),
                                     yaxis = list(title = 'DC 2'),
                                     zaxis = list(title = 'DC 3')))
fig2
```

<p align="center">
  <img src="docs/mirna breast comp.png" width="500">
</p>

```{r}
fig3 <- plot_ly(proteinComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`, color = ~subtype, colors = cols)
fig3 <- fig3 %>% add_markers()
fig3 <- fig3 %>% layout(scene = list(xaxis = list(title = 'DC 1'),
                                     yaxis = list(title = 'DC 2'),
                                     zaxis = list(title = 'DC 3')))
fig3
```

<p align="center">
  <img src="docs/protein breast comp.png" width="500">
</p>



### Support

For any question, request or bug report please create a new issue in this repository. 

### Contributions

We welcome contributions and suggestions from the community. If you have any idea, please submit it as an issue, which we will look into and ask for further explannation if necessary. 

### Reporting issues

PanDA is under continuous development. If you encounter an issue, please make sure you install all the required packages necessary to run the codes. If that does not solve your problem, please open a new issue detailing your encountered problem by providing a code and a demo example. We will try to look into your issue and hopefully provide you with a solution. Thanks.

##  Cite PanDA
To cite this work and access a comprehensive description of the PanDA method, please refer to:
