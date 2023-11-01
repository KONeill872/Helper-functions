# Differential expression analysis methods 

The script contains wrapper functions for: 
- DESeq2 (Wald test): `DEA_DESeq2(...)`
- Limma trend: `DEA_limma(...)`
- Limma voom: `DEA_limma(..., voom=T)`
    - With one random effect `DEA_limma(..., voom=T, random=T)`
- Limma voom with quality weights: `DEA_limma(..., voomWeights =T)`
    - With one random effect `DEA_limma(..., voomWeights=T, random=T)`
- edgeR (LTR): `DEA_edgeR(...)`
- Linear mixed model `DEA_glmm(..., lmm=T, glmm=F)`
- General linear mixed model: `DEA_glmm(...)`

### Inputs requirements

#### Counts data
All methods required DGEList for input. 

#### Design table 
Dataframe containing metadata for the experiment. Sample should match names in DGEList.

| sample   | batch | condition |
|----------|-------|-----------|
| sample A | 1     | T         |
| sample B | 1     | N         |
| sample C | 2     | T         |
| sample D | 2     | N         |
| sample E | 3     | T         |
| sample F | 3     | N         |

####  Model 
Example model for DESeq2, Limma and edgeR methods: `model <- ~ batch + condition`  
Example model for LMM or GLMM (requires random effect in the model): `model <- (1| batch) + condition`  

#### Can add gene name annotation from Ensembl annotation
Insert gene annotation to `geneNames`.

#### Output files can be used in `plotFunctions.R` for PCA, volcano, MA and heatmap. 
