# Paths
MM10 = /nfs/turbo/umms-thahoang/sherine/REFERENCES/refdata-cellranger-arc-mm10-2020-A-2.0.0
FASTQS = /nfs/turbo/umms-thahoang/sherine/Neurog2_multiomics/align
CELLRANGER = /nfs/turbo/umms-thahoang/sherine/tools/cellranger-arc-2.0.2
# Default target: build both samples
all: TH1 TH2

# Rules for each sample with increased cores and memory
TH1:
	$(CELLRANGER)/cellranger-arc count \
		--id=TH1 \
		--reference=$(MM10) \
		--fastqs=$(FASTQS) \
		--sample=13793-TH-1-GEX,13793-TH-1-ATAC \
		--localcores=16 \
		--localmem=48

TH2:
	$(CELLRANGER)/cellranger-arc count \
		--id=TH2 \
		--reference=$(MM10) \
		--fastqs=$(FASTQS) \
		--sample=13793-TH-2-GEX,13793-TH-2-ATAC \
		--localcores=16 \
		--localmem=48

# Optional: clean outputs
clean:
	rm -rf TH1 TH2

