#!/usr/bin/env nextflow

nextflow.enable.dsl=2
markers = ['Myh7', 'Myh1', 'Myh2', 'Vwf', 'Pax7', 'Ptprc', 'Myh11', 'Myh3', 'Myh4', 'Fbn1', 'Pdgfra', 'Acta2', 'Chrna1']


process subset_nuclei {

    publishDir "${params.results}/pass-qc-nuclei-counts"
    memory '50 GB'
    tag "${library}"
    container 'docker://porchard/mm:20230104'

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv'), path('keep-barcodes.txt')

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    mm subset --matrix matrix.mtx --features features.tsv --barcodes barcodes.tsv --keep-barcodes keep-barcodes.txt --prefix ${library}.
    """

}


process get_soup {

    publishDir "${params.results}/soup"
    memory '50 GB'
    tag "${library}"
    container 'docker://porchard/mm:20230104'

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv')

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    mm countbarcodes matrix.mtx features.tsv barcodes.tsv > barcode-counts.txt
    cat barcode-counts.txt | awk '\$2<=${params.libraries[library]['thresholds']['soup_max_umi']}' | cut -f1 > soup-barcodes.txt
    mm subset --matrix matrix.mtx --features features.tsv --barcodes barcodes.tsv --keep-barcodes soup-barcodes.txt --prefix ${library}.
    """

}


process initial_cluster {

    publishDir "${params.results}/initial-cluster"
    memory '15 GB'
    tag "${library}"
    container 'library://porchard/default/r-general:20220112'

    input:
    tuple val(library), path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    tuple val(library), path("${library}.clusters.txt"), emit: clusters

    """
    seurat.R --markers ${markers.join(',')} --prefix ${library}. --sctransform --pcs 15 --resolution 0.2 --matrix $matrix --features $features --barcodes $barcodes
    """

}


process decontX {

    publishDir "${params.results}/decontX"
    memory '25 GB'
    tag "${library}"
    container 'docker://porchard/decontx:20230104'

    input:
    tuple val(library), path('nuclei.matrix.mtx'), path('nuclei.features.tsv'), path('nuclei.barcodes.tsv'), path('soup.matrix.mtx'), path('soup.features.tsv'), path('soup.barcodes.tsv'), path(clusters)

    output:
    tuple val(library), path("${library}.decontaminated-counts-rounded.matrix.mtx"), path("${library}.decontaminated-counts-rounded.features.tsv"), path("${library}.decontaminated-counts-rounded.barcodes.tsv"), emit: decontaminated
    tuple val(library), path("${library}.contamination.txt"), emit: contamination
    path("*.png")

    """
    decontX.R --matrix nuclei.matrix.mtx --features nuclei.features.tsv --barcodes nuclei.barcodes.tsv --soup_matrix soup.matrix.mtx --soup_features soup.features.tsv --soup_barcodes soup.barcodes.tsv --prefix ${library}. --clusters $clusters
    """

}


process filter_decontX {

    publishDir "${params.results}/decontX-filtered"
    memory '50 GB'
    tag "${library}"
    container 'docker://porchard/mm:20230104'

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv'), path(contamination)

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    cat $contamination | awk '\$2<=${params.libraries[library]['thresholds']['max_contamination']}' | cut -f1 > keep-barcodes.txt
    mm subset --matrix matrix.mtx --features features.tsv --barcodes barcodes.tsv --keep-barcodes keep-barcodes.txt --prefix ${library}.
    """

}


process cluster_no_contamination_filter {

    publishDir "${params.results}/cluster/no-contamination-filter"
    memory '15 GB'
    tag "${library}"
    container 'library://porchard/default/r-general:20220112'

    input:
    tuple val(library), path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    tuple val(library), path("${library}.clusters.txt"), emit: clusters

    """
    seurat.R --markers ${markers.join(',')} --prefix ${library}. --pcs 15 --resolution 0.2 --matrix $matrix --features $features --barcodes $barcodes
    """

}


process cluster_contamination_filter {

    publishDir "${params.results}/cluster/contamination-filter"
    memory '15 GB'
    tag "${library}"
    container 'library://porchard/default/r-general:20220112'

    input:
    tuple val(library), path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    tuple val(library), path("${library}.clusters.txt"), emit: clusters

    """
    seurat.R --markers ${markers.join(',')} --prefix ${library}. --pcs 15 --resolution 0.2 --matrix $matrix --features $features --barcodes $barcodes
    """

}


process add_prefix_to_barcodes {

    memory '50 GB'
    tag "${library}"

    input:
    tuple val(library), path('matrix.mtx'), path('features.tsv'), path('barcodes.tsv')

    output:
    tuple val(library), path("${library}.matrix.mtx"), path("${library}.features.tsv"), path("${library}.barcodes.tsv")

    """
    ln -s matrix.mtx ${library}.matrix.mtx
    ln -s features.tsv ${library}.features.tsv
    cat barcodes.tsv | perl -pe 's/^/${library}-/' > ${library}.barcodes.tsv
    """

}

process merge_matrices {

    publishDir "${params.results}/merged-counts"
    memory '50 GB'
    container 'docker://porchard/mm:20230104'

    input:
    tuple val(x), path(matrices), path(features), path(barcodes)

    output:
    tuple path("merged.matrix.mtx"), path("merged.features.tsv"), path("merged.barcodes.tsv")

    """
    mm merge --matrices ${matrices.join(' ')} --features ${features.join(' ')} --barcodes ${barcodes.join(' ')} --prefix merged.
    """

}

process cluster_joint {

    publishDir "${params.results}/cluster/joint"
    memory '50 GB'
    container 'library://porchard/default/r-general:20220112'

    input:
    tuple path(matrix), path(features), path(barcodes)

    output:
    path("*.png")
    path("*.txt")
    path("merged.clusters.txt")

    """
    seurat.R --markers ${markers.join(',')} --prefix merged. --pcs 15 --resolution 0.2 --matrix $matrix --features $features --barcodes $barcodes
    """

}


workflow {

    libraries = params.libraries.keySet()

    library_mtx_features_barcodes_passqcbarcodes = Channel.from(libraries.collect({it -> [it, file(params.libraries[it].matrix), file(params.libraries[it].features), file(params.libraries[it].barcodes), file(params.libraries[it].pass_qc_barcodes)]}))
    nuclei = subset_nuclei(library_mtx_features_barcodes_passqcbarcodes)
    first_clusters = initial_cluster(nuclei)
    soup = get_soup(library_mtx_features_barcodes_passqcbarcodes.map({it -> it[0..(it.size() - 2)]}))
    d = nuclei.combine(soup, by: 0).combine(first_clusters.clusters, by: 0) | decontX
    cluster_no_contamination_filter(d.decontaminated)
    filtered = d.decontaminated.combine(d.contamination, by: 0) | filter_decontX
    cluster_contamination_filter(filtered)
    add_prefix_to_barcodes(filtered).map({it -> ['x'] + it[1..3]}).groupTuple() | merge_matrices | cluster_joint

}
