ORDER OF RUNNING CODE:

(For data downloading script, please see `/scratch/st-kdkortha-1/nshen7/vmrseq/vmrseq-experiments/data/raw_counts/raw_counts_Liu2021/README.txt`)

0. `metadata_eda`: EDA on metadata

1. `raw_data_subsetting`: removes non-CpG sites in raw data to reduce occupied space in disk 

2. `data_processing` and `data_processing2` at the same time: integrate individual cell info into SummarizedExperiment object (according to cell type and sample)

3. `data_qc`: perform QC