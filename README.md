# Neurog2 Multiomics 


## QC before filtering 



Sample distribution changes:
   Sample Pre_Filter Post_Filter Retention_Rate
1 Control      11409        7283          63.84
2      OE      12560        7561          60.20




## After filtering 

We used the following filtering parameters:

| Parameter       | Description                           | Default Value |
|-----------------|---------------------------------------|---------------|
| `minTSS`        | Minimum TSS enrichment for ATAC cells | 10            |
| `minFrags`      | Minimum number of ATAC fragments      | 5,000         |
| `minGexUMI`     | Minimum number of RNA UMIs per cell   | 1,000         |
| `maxGexUMI`     | Maximum number of RNA UMIs per cell   | 30,000        |
| `minGexGenes`   | Minimum number of genes detected      | 500           |
| `maxGexGenes`   | Maximum number of genes detected      | 7,000         |




### Number of cells per sample after filtering

| Sample  | Cells |
|---------|-------|
| Control |  7238 |
| OE      |  7561 |



![After filtering QC](Neurog2_postFilterQC.png?v=5)

## Adding UMAP and clustering 

![Clsuters_Samples_Sample](Neurog2_SamplesUMAP_bySample.png?v=5)

![Clusters_Samples_ClustersUMAP](Neurog2_SamplesUMAP_byCluster.png?v=5)

![ClustersUMI](Neurog2_perClustersnUMI.png?v=5)



# Marker Genes

<img src="Neurog2_UMAP_Zscore_Hes5.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Vim.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Abca8a.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Slc1a3.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Prdx6.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Pax6.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Aqp4.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Notch1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Hes1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Apoe.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Gfap.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Rlbp1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Lhx4.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Slc6a9.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Slc17a7.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Bsn.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Tfap2a.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Neurog2.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Emx1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Atoh7.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Sox11.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Insm1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Otx2.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Olig2.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Prdm1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Chat.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Foxn4.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Ascl1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Elavl4.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Isl1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Cabp5.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Elavl3.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Gad1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Sebox.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Calb1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Calb2.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Csf1r.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Rbfox3.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx2.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Glul.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Sox9.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Malat1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Arr3.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Nrl.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Rho.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Acta2.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Pou4f2.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Tie1.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Kcnj8.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Lhx1.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Rpe65.png?v=5" width="200">

<img src="Neurog2_UMAP_Zscore_Pax2.png?v=5" width="200"> <img src="Neurog2_UMAP_Zscore_Ccr2.png?v=5" width="200">




# ⚠️🚨 **ATTENTION! IMPORTANT NOTICE** WILL NEED TO RE-ANNOTATE 🚨⚠️

# Annotations 


![Annotations](Neurog2_annotations.png)


# DGE Heatmap 


![HEATMAP](Neurog2_heatmap.png?v=2)


# Peaks 

![PEAKS](peaks.png) 

# Motifs

![Motifs](motifs.png) 



 
