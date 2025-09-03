# Neurog2 Multiomics 

## Before Filtering QC Summary

**Cell types:**

| Cell type  | Number of cells |
|------------|----------------|
| ATAC_only  | 2              |
| ATAC+RNA   | 31             |

**Number of cells per sample before filtering:**

| Sample   | Number of cells |
|----------|----------------|
| Control  | 20             |
| KO       | 13             |

## Before Filtering QC Plots

The figure below shows QC metrics for the Neurog2 project before filtering:

![Neurog2 Before Filtering QC](Neurog2_beforeFilterQC.png)

## No. of cells 

```r
# 1️⃣ Clean ATAC barcodes from proj_ALL
atac_clean <- gsub(".*#", "", rownames(proj_ALL@cellColData))
cat("Number of ATAC cells:", length(atac_clean), "\n")

# 2️⃣ Clean RNA barcodes from the combined SE
rna_clean <- gsub("-1$", "", colnames(seRNAcombined))
cat("Number of RNA barcodes:", length(rna_clean), "\n")

# 3️⃣ Optional: count how many are in both
common_cells <- intersect(atac_clean, rna_clean)
cat("Number of barcodes present in both ATAC and RNA:", length(common_cells), "\n")


| Metric                               | Count   |
|--------------------------------------|--------|
| Number of ATAC cells                  | 33     |
| Number of RNA barcodes                | 22,232 |
| Number of barcodes present in both    | 0      |



## Cell Ranger html 

[TH1](TH1.html)
[TH2](TH2.html)

