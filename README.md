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

![Pre filtering QC](Neurog2_beforeFilterQC.png)

## After filtering 

We used the following filtering parameters:

We used the following filtering parameters:

| Parameter       | Description                           | Default Value |
|-----------------|---------------------------------------|---------------|
| `minTSS`        | Minimum TSS enrichment for ATAC cells | 10            |
| `minFrags`      | Minimum number of ATAC fragments      | 5,000         |
| `minGexUMI`     | Minimum number of RNA UMIs per cell   | 1,000         |
| `maxGexUMI`     | Maximum number of RNA UMIs per cell   | 30,000        |
| `minGexGenes`   | Minimum number of genes detected      | 500           |
| `maxGexGenes`   | Maximum number of genes detected      | 7,000         |



### Number of cells per sample before filtering

| Sample  | Cells |
|---------|-------|
| Control | 11409 |
| KO      | 12560 |

### Number of cells per sample after filtering

| Sample  | Cells |
|---------|-------|
| Control | 7282  |
| KO      | 7560  |




![After filtering QC](Neurog2_postFilterQC.png)



## Adding UMAP and clustering 

![Clsuters_Samples_Sample](Neurog2_SamplesUMAP_bySample.png)

![Clusters_Samples_ClustersUMAP](Neurog2_SamplesUMAP_byCluster.png)

![ClustersUMI](Neurog2_perClustersnUMI.png)


# Marker Genes



<img src="Neurog2_UMAP_Zscore_Sox9.png" width="200"> <img src="Neurog2_UMAP_Zscore_Rho.png" width="200"> <img src="Neurog2_UMAP_Zscore_Gad1.png" width="200">

<img src="Neurog2_UMAP_Zscore_Emx1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Pax6.png" width="200"> <img src="Neurog2_UMAP_Zscore_Apoe.png" width="200">

<img src="Neurog2_UMAP_Zscore_Ascl1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Chat.png" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx2.png" width="200">

<img src="Neurog2_UMAP_Zscore_Otx2.png" width="200"> <img src="Neurog2_UMAP_Zscore_Neurog2.png" width="200"> <img src="Neurog2_UMAP_Zscore_Ccr2.png" width="200">

<img src="Neurog2_UMAP_Zscore_Prdx6.png" width="200"> <img src="Neurog2_UMAP_Zscore_Gfap.png" width="200"> <img src="Neurog2_UMAP_Zscore_Elavl4.png" width="200">

<img src="Neurog2_UMAP_Zscore_Lhx1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Notch1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Bsn.png" width="200">

<img src="Neurog2_UMAP_Zscore_Elavl3.png" width="200"> <img src="Neurog2_UMAP_Zscore_Tie1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Vim.png" width="200">

<img src="Neurog2_UMAP_Zscore_Slc17a7.png" width="200"> <img src="Neurog2_UMAP_Zscore_Sox11.png" width="200"> <img src="Neurog2_UMAP_Zscore_Acta2.png" width="200">

<img src="Neurog2_UMAP_Zscore_Rlbp1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Prdm1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Calb1.png" width="200">

<img src="Neurog2_UMAP_Zscore_Rpe65.png" width="200"> <img src="Neurog2_UMAP_Zscore_Malat1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx4.png" width="200">

<img src="Neurog2_UMAP_Zscore_Insm1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Kcnj8.png" width="200"> <img src="Neurog2_UMAP_Zscore_Slc1a3.png" width="200">

<img src="Neurog2_UMAP_Zscore_Nrl.png" width="200"> <img src="Neurog2_UMAP_Zscore_Arr3.png" width="200"> <img src="Neurog2_UMAP_Zscore_Foxn4.png" width="200">

<img src="Neurog2_UMAP_Zscore_Hes1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Isl1.png" width="200"> <img src="Neurog2_UMAP_Zscore_Calb2.png" width="200">

<img src="Neurog2_UMAP_Zscore_Pou4f2.png" width="200"> <img src="Neurog2_UMAP_Zscore_Hes5.png" width="200"> <img src="Neurog2_UMAP_Zscore_Cabp5.png" width="200">

<img src="Neurog2_UMAP_Zscore_Sebox.png" width="200"> <img src="Neurog2_UMAP_Zscore_Atoh7.png" width="200"> <img src="Neurog2_UMAP_Zscore_Abca8a.png" width="200">

<img src="Neurog2_UMAP_Zscore_Slc6a9.png" width="200"> <img src="Neurog2_UMAP_Zscore_Olig2.png" width="200"> <img src="Neurog2_UMAP_Zscore_Aqp4.png" width="200">

<img src="Neurog2_UMAP_Zscore_Rbfox3.png" width="200"> <img src="Neurog2_UMAP_Zscore_Tfap2a.png" width="200"> <img src="Neurog2_UMAP_Zscore_Glul.png" width="200">

<img src="Neurog2_UMAP_Zscore_Csf1r.png" width="200"> <img src="Neurog2_UMAP_Zscore_Pax2.png" width="200">



