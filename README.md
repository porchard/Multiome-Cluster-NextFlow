# Notes

snRNA and QC pipelines must be run first. Each library must be from a single species (i.e., no barnyard experiments).

# Steps

First, update nextflow.config.

## Per-library RNA decontamination (decontX) and clustering + joint clustering

This step generates the cleaned RNA count matrices to be used for clustering, and performs the clustering on a per-library and joint basis.

Generate the config:

```bin
python bin/make-clean-and-cluster-config.py /path/to/rnaseq/results /path/to/qc/results > config.json
```

The config should look like this:

```python
{
    "libraries": {
        "6616-NM-1-mm10": {
            "barcodes": "/lab/work/porchard/PL6616/work/rnaseq/results/starsolo/6616-NM-1-mm10/6616-NM-1-mm10.Solo.out/GeneFull_ExonOverIntron/raw/barcodes.tsv", # list of barcodes (corresponding to the market matrix file of RNA UMI counts). Should contain ALL barcodes, not just pass QC barcodes
            "features": "/lab/work/porchard/PL6616/work/rnaseq/results/starsolo/6616-NM-1-mm10/6616-NM-1-mm10.Solo.out/GeneFull_ExonOverIntron/raw/features.tsv", # list of genes (corresponding to the market matrix file of RNA UMI counts)
            "matrix": "/lab/work/porchard/PL6616/work/rnaseq/results/starsolo/6616-NM-1-mm10/6616-NM-1-mm10.Solo.out/GeneFull_ExonOverIntron/raw/matrix.mtx", # market matrix format file of RNA UMI counts, corresponding to the barcodes and features above.
            "pass_qc_barcodes": "/lab/work/porchard/PL6616/work/qc/results/atac-doublet-detection/6616-NM-1-mm10.rna.singlets.txt", # list of pass QC barcodes (should pass QC thresholds and be singlets)
            "thresholds": {
                "max_contamination": "0.2", # nuclei with estimated ambient contamination > this value (based on decontX output) are filtered out
                "soup_max_umi": "10" # barcodes with <= this number of total UMIs in the RNA UMI count matrix are taken to represent ambient background (passed to decontX as background matrix)
            }
        }
    }
}
```

Edit the config if desired, then run the pipeline:

```bin
nextflow run -resume -with-trace -params-file config.json -with-report -qs 300 --results /path/to/results /path/to/Multiome-Cluster-NextFlow/clean-and-cluster.nf
```

You may wish to run under nohup so that the pipeline continues to run in the background and does not terminate upon logging out of the server (`nohup nextflow run ... &`)

## Output
* `soup/*`: Market matrix files for each library, representing the RNA UMI soup matrix (background matrix for decontX)
* `pass-qc-nuclei-counts/*`: Market matrix files for each library, representing the RNA UMI matrix of only pass-QC barcodes
* `initial-cluster/*`: Output of Seurat clustering on RNA UMI matrix of pass-QC barcodes, prior to running decontX. This clustering is used as the clustering for decontX. **You should check that this looks reasonable as it may substantially impact the quality of decontX results**
* `decontX/*`: Output of decontX, including decontaminated matrices.
* `decontX-filtered/*`: DecontX decontaminated matrices, after removing nuclei with high estimated contamination levels (i.e. applying the max_contamination threshold from the config)
* `merged-counts/*`: Merged market matrix RNA UMI matrix (merged across all libraries; i.e. merging the matrices in `decontX-filtered`). A prefix is added to each barcode to distinguish between identical barcodes from different libraries.
* `cluster/no-contamination-filter/*`: Output of Seurat clustering, using decontX decontaminated matrices but prior to removing nuclei with high estimated contamination levels (i.e., using matrices in `decontX/`)
* `cluster/contamination-filter/*`: Output of Seurat clustering, using decontX decontaminated matrices after removing nuclei with high estimated contamination levels (i.e., using matrices in `decontX-filtered/`)
* `cluster/joint/*`: Output of joint Seurat clustering across all libraries, using decontX decontaminated matrices after removing nuclei with high estimated contamination levels (i.e., using matrices in `merged-counts/`)