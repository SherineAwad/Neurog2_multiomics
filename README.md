# Neurog2 Multiomics 

# re-align FASTQ files using cellranger-arc to mm10 

## Cell ranger web summary 

[TH1 web summary](TH1_mm10.html)


[TH2 web summary](TH2_mm10.html)


| Metric                              | Count  |
|-------------------------------------|--------|
| Number of ATAC cells                | 23,969 |
| Number of RNA barcodes              | 22,105 |

## Before filtering 


| cell_type | ATAC_only | ATAC+RNA |
|-----------|-----------|----------|
| **Total** | 2,305     | 21,664   |



| Condition | Control | KO    |
|-----------|---------|-------|
| **Total** | 11,409  |12,560 |


## QC before filtering 

![Pre filtering QC](mNeurog2_beforeFilterQC.png)

## After filtering 

We used the following filtering parameters:

| Parameter       | Description                           | Default Value |
|-----------------|---------------------------------------|---------------|
| `minTSS`        | Minimum TSS enrichment for ATAC cells | 10            |
| `minFrags`      | Minimum number of ATAC fragments      | 5,000         |
| `minGexUMI`     | Minimum number of RNA UMIs per cell   | 1,000         |
| `maxGexUMI`     | Maximum number of RNA UMIs per cell   | 15,000        |
| `minGexGenes`   | Minimum number of genes detected      | 500           |
| `maxGexGenes`   | Maximum number of genes detected      | 5,000         |


Number of cells per sample after filtering:

| Condition | Control | KO    |
|-----------|---------|-------|
|**Total**  |  7179   | 7205  |



![After filtering QC](mNeurog2_postFilterQC.png)



## Adding UMAP and clustering 

![Clsuters_Samples_Sample](mNeurog2_SamplesUMAP_bySample.png)

![Clusters_Samples_ClustersUMAP](mNeurog2_SamplesUMAP_byCluster.png)

![ClustersUMI](mNeurog2_perClustersnUMI.png)




