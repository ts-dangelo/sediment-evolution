nextflow.enable.dsl=2

params.faa_dir   = "/home/tdangelo/epc/prodigal_mmseqs/aero_misag_gene_faas_sc/"
params.aln_dir   = "/home/tdangelo/epc/prodigal_mmseqs/aero_misag_gene_aln_sc"
params.ffn_dir   = "/home/tdangelo/epc/prodigal_mmseqs/aero_misag_gene_ffn_sc"
params.codon_dir = "/home/tdangelo/epc/prodigal_mmseqs/aero_misag_gene_codon_sc"
params.tree_dir  = "/home/tdangelo/epc/prodigal_mmseqs/aero_misag_gene_aln_sc_iqtree"
params.hyphy_out = "/home/tdangelo/epc/hyphy"

params.devel = false

workflow {

    file(params.aln_dir).mkdirs()
    file(params.tree_dir).mkdirs()
    file(params.codon_dir).mkdirs()
    file("${params.hyphy_out}/fubar").mkdirs()
    file("${params.hyphy_out}/fel").mkdirs()
    file("${params.hyphy_out}/absrel").mkdirs()
    file("${params.hyphy_out}/busted").mkdirs()

    // Source channel (optional devel)
    def faa_files = params.devel ?
        Channel.fromPath("${params.faa_dir}/*.faa").take(3) :
        Channel.fromPath("${params.faa_dir}/*.faa")

    // Run mafft
    def mafft_ch = faa_files
        .map { file -> tuple(file.baseName, file) }
        | mafft_align

    // Run IQ-TREE
    def tree_ch = mafft_ch | iqtree_aa

    // Run pal2nal
    def codon_ch = mafft_ch | pal2nal

    // Join trees + codons for HyPhy
    def hyphy_input = tree_ch.join(codon_ch)

    // Run HyPhy
    hyphy_input | hyphy_fubar
    hyphy_input | hyphy_fel
    hyphy_input | hyphy_absrel
    hyphy_input | hyphy_busted
}


process mafft_align {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    publishDir params.aln_dir, mode: 'copy'

    input:
    tuple val(cluster), path(faa)

    output:
    tuple val(cluster), path("${cluster}.aln")

    script:
    """
    mafft --auto ${faa} > ${cluster}.aln
    """
}


process iqtree_aa {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    errorStrategy 'ignore'

    publishDir params.tree_dir, mode: 'copy'

    input:
    tuple val(cluster), path(aln)

    output:
    tuple val(cluster),
          path("${cluster}.aln.contree")

    script:
    """
    iqtree \
      -s ${aln} \
      -st AA \
      -m TEST \
      -bb 1000 \
      -nt AUTO \
      > ${cluster}.iqtree.log 2>&1
    """
}



process pal2nal {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    publishDir params.codon_dir, mode: 'copy'

    input:
    tuple val(cluster), path(aln)

    output:
    tuple val(cluster),
          path("${cluster}.codon.fna")

    script:
    """
    ffn="${params.ffn_dir}/${cluster}.ffn"

    pal2nal.pl ${aln} \$ffn -output fasta -nomismatch \
        > ${cluster}.codon.fna

    """
}


process hyphy_fubar {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    errorStrategy 'ignore'

    publishDir "${params.hyphy_out}/fubar", mode: 'copy'

    cpus 8
    memory '4 GB'

    input:
    tuple val(cluster),
          path(tree),
          path(codon)

    output:
    path("${cluster}.fubar.json")
    path("${cluster}.fubar.stdout")

    script:
    """
    hyphy fubar \
        --alignment ${codon} \
        --tree ${tree} \
        --output ${cluster}.fubar.json \
        > ${cluster}.fubar.stdout
    """
}


process hyphy_fel {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    errorStrategy 'ignore'

    publishDir "${params.hyphy_out}/fel", mode: 'copy'

    cpus 8
    memory '4 GB'

    input:
    tuple val(cluster),
          path(tree),
          path(codon)

    output:
    path("${cluster}.fel.json")
    path("${cluster}.fel.stdout")

    script:
    """
    hyphy fel \
        --alignment ${codon} \
        --tree ${tree} \
        --output ${cluster}.fel.json \
        > ${cluster}.fel.stdout
    """
}


process hyphy_absrel {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    errorStrategy 'ignore'

    publishDir "${params.hyphy_out}/absrel", mode: 'copy'

    cpus 8
    memory '8 GB'

    input:
    tuple val(cluster),
          path(tree),
          path(codon)

    output:
    path("${cluster}.absrel.json")
    path("${cluster}.absrel.stdout")

    script:
    """
    hyphy absrel \
        --alignment ${codon} \
        --tree ${tree} \
        --output ${cluster}.absrel.json \
        > ${cluster}.absrel.stdout
    """
}

process hyphy_busted {

    label 'hyphy'
    container '/mnt/s1/labs/lindsay/epc/nextflow/hyphy.sif'

    tag "$cluster"

    errorStrategy 'ignore'

    publishDir "${params.hyphy_out}/busted", mode: 'copy'

    cpus 8
    memory '8 GB'

    input:
    tuple val(cluster),
          path(tree),
          path(codon)

    output:
    path("${cluster}.busted.json")
    path("${cluster}.busted.stdout")

    script:
    """
    hyphy busted \
        --alignment ${codon} \
        --tree ${tree} \
        --output ${cluster}.busted.json \
        > ${cluster}.busted.stdout
    """
}
