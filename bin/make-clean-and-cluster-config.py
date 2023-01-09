#!/usr/bin/env python
# coding: utf-8

import sys
import glob
import json
import os

RNA_RESULTS = sys.argv[1]
QC_RESULTS = sys.argv[2]

STARSOLO_DIRS = glob.glob(os.path.join(RNA_RESULTS, 'starsolo', '*'))
PASS_QC_BARCODES = glob.glob(os.path.join(QC_RESULTS, 'atac-doublet-detection', '*.rna.singlets.txt'))

# infer library names
rna_matrix = {os.path.basename(f): os.path.join(f, os.path.basename(f) + '.Solo.out', 'GeneFull_ExonOverIntron', 'raw', 'matrix.mtx') for f in STARSOLO_DIRS}
rna_features = {os.path.basename(f): os.path.join(f, os.path.basename(f) + '.Solo.out', 'GeneFull_ExonOverIntron', 'raw', 'features.tsv') for f in STARSOLO_DIRS}
rna_barcodes = {os.path.basename(f): os.path.join(f, os.path.basename(f) + '.Solo.out', 'GeneFull_ExonOverIntron', 'raw', 'barcodes.tsv') for f in STARSOLO_DIRS}
pass_qc_barcodes = {os.path.basename(f).split('.')[0]: f for f in PASS_QC_BARCODES}

# ensure that all keys are same
assert(set(rna_matrix.keys()) == set(rna_features.keys()))
assert(set(rna_matrix.keys()) == set(rna_barcodes.keys()))
assert(set(rna_matrix.keys()) == set(pass_qc_barcodes.keys()))

libraries = {
    library: {
        'matrix': rna_matrix[library],
        'features': rna_features[library],
        'barcodes': rna_barcodes[library],
        'pass_qc_barcodes': pass_qc_barcodes[library],
        'thresholds': {
            'soup_max_umi': '10',
            'max_contamination': '0.2'
        }
    } 
    for library in rna_matrix.keys()
}


CONFIG = {'libraries': libraries}
print(json.dumps(CONFIG, sort_keys = True, indent = 4))