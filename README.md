# Neurog2 Multiomics 


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

![Pre filtering QC](Neurog2_preFilterQC.png)

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

## Per Clusters QC 

![Per Clusters QC](Neurog2_perClusterQC.png)

# Marker Genes

<img src="Neurog2_UMAP_Zscore_Sox9.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Rho.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Gad1.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Emx1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Pax6.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Apoe.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Ascl1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Chat.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx2.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Otx2.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Neurog2.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Ccr2.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Prdx6.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Gfap.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Elavl4.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Lhx1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Notch1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Bsn.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Elavl3.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Tie1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Vim.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Slc17a7.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Sox11.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Acta2.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Rlbp1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Prdm1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Calb1.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Rpe65.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Malat1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx4.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Insm1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Kcnj8.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Slc1a3.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Nrl.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Arr3.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Foxn4.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Hes1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Isl1.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Calb2.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Pou4f2.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Hes5.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Cabp5.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Sebox.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Atoh7.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Abca8a.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Slc6a9.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Olig2.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Aqp4.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Rbfox3.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Tfap2a.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Glul.png?v=2" width="200">

<img src="Neurog2_UMAP_Zscore_Csf1r.png?v=2" width="200"> <img src="Neurog2_UMAP_Zscore_Pax2.png?v=2" width="200">


# Refining 

### After removing clusters 1, 2, 3, 4, 5, 6, and 21

<img src="Neurog2_refined_ClustersUMAP.png?v=2" width =600> 

<img src="Neurog2_CellTypeUMAP_Filtered.png?v=2" width =600>


 
