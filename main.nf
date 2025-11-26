#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UHVDB/nucleotidecluster
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/UHVDB/nucleotidecluster
----------------------------------------------------------------------------------------
*/


// Define processes
process SEQKIT_SPLIT2 {
    label 'process_high'

    input:
    tuple val(meta), path(fna)

    output:
    tuple val(meta), path("split_fnas/*")   , emit: fnas
    tuple val(meta), path(".command.log")   , emit: log
    tuple val(meta), path(".command.sh")    , emit: script

    script:
    """
    seqkit split2 \\
        ${fna} \\
        --threads ${task.cpus} \\
        --by-size ${params.chunk_size} \\
        --out-dir split_fnas \\
        --extension ".gz"
    """
}

process KMERDB_BUILD {
    label 'process_super_high'

    input:
    tuple val(meta), path(fna)

    output:
    tuple val(meta), path("ref.kdb")        , emit: kmerdb
    tuple val(meta), path(".command.log")   , emit: log
    tuple val(meta), path(".command.sh")    , emit: script

    script:
    """
    echo "${fna}" > ref_kdb.txt

    # build kmer-db database for reference
    kmer-db \\
        build \\
        -k 25 \\
        -f 0.2 \\
        -t ${task.cpus} \\
        -multisample-fasta \\
        ref_kdb.txt \\
        ref.kdb
    """
}

process ALIGN {
    label 'process_super_high'

    input:
    tuple val(meta) , path(chunk_fna)
    tuple val(meta2), path(kmerdb), path(full_fna)

    output:
    tuple val(meta), path("${meta.id}.vclust_ani.parquet.zst")      , emit: parquet
    tuple val(meta), path("${meta.id}.vclust_ani.ids.parquet.zst")  , emit: ids
    tuple val(meta), path(".command.log")                           , emit: log
    tuple val(meta), path(".command.sh")                            , emit: script

    script:
    """
    echo "${chunk_fna}" > query_kdb.txt

    # build kmer-db database for query
    kmer-db \\
        new2all \\
        -sparse \\
        -min num-kmers:20 \\
        -min ani-shorter:${params.min_ani} \\
        -t ${task.cpus} \\
        -multisample-fasta \\
        ${kmerdb} \\
        query_kdb.txt \\
        ${meta.id}.new2all.csv

    # convert kmer-db output to distances
    kmer-db \\
        distance \\
        ani-shorter \\
        -sparse \\
        -min ${params.min_ani} \\
        -t ${task.cpus} \\
        ${meta.id}.new2all.csv \\
        ${meta.id}.dist.csv

    # convert distances to lz-ani format
    convert_to_lzani.py \\
        -i ${meta.id}.dist.csv \\
        -o ${meta.id}.lzani_prefilter.txt

    # align sequences with LZ-ANI
    lz-ani \\
        all2all \\
        --in-fasta ${full_fna} \\
        -o ${meta.id}.vclust_ani.tsv \\
        --out-format query,reference,tani,gani,ani,qcov,rcov \\
        -t ${task.cpus} \\
        --multisample-fasta true \\
        --out-type tsv \\
        --out-filter ani ${params.min_ani} \\
        --flt-kmerdb ${meta.id}.lzani_prefilter.txt ${params.min_qcov}

    tsv_to_parquet.py --input ${meta.id}.vclust_ani.tsv --output ${meta.id}.vclust_ani.parquet.zst
    tsv_to_parquet.py --input ${meta.id}.vclust_ani.ids.tsv --output ${meta.id}.vclust_ani.ids.parquet.zst
    rm ${meta.id}.new2all.csv ${meta.id}.dist.csv ${meta.id}.lzani_prefilter.txt
    """
}

process COMBINEANIS {
    label 'process_medium'
    publishDir path: file("${params.output}").getParent(), mode:'copy'

    input:
    tuple val(meta), path(parquets)

    output:
    tuple val(meta), path("${output}.ani.tsv.zst")  , emit: tsv
    tuple val(meta), path(".command.log")           , emit: log
    tuple val(meta), path(".command.sh")            , emit: script

    script:
    output = file("${params.output}").getName()
    """
    combine_anis.py \\
        -i "./*.parquet.zst" \\
        -o ${output}.ani.tsv

    zstd ${output}.ani.tsv
    """
}

process MCL {
    label 'process_super_high'
    publishDir path: file("${params.output}").getParent(), mode:'copy'

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${output}.mcl.zst")  , emit: mcl
    tuple val(meta), path(".command.log")       , emit: log
    tuple val(meta), path(".command.sh")        , emit: script

    script:
    output = file("${params.output}").getName()
    """
    # cut columns for MCL
    csvtk cut \\
        ${tsv} \\
        --tabs \\
        --fields query,reference,gani \\
        --delete-header \\
        --num-cpus ${task.cpus} \\
        --out-tabs \\
        --out-file ${output}.mcl_input.tsv

    # run mcl
    mcl \\
        ${output}.mcl_input.tsv \\
        --abc \\
        -sort revsize \\
        -te ${task.cpus} \\
        -o ${output}.mcl

    zstd ${output}.mcl
    """
}

process CLUSTY {
    label 'process_super_high'
    publishDir path: file("${params.output}").getParent(), mode:'copy'

    input:
    tuple val(meta) , path(tsv)
    tuple val(meta2), path(parquet)

    output:
    tuple val(meta), path("${output}.clusty.tsv.zst")   , emit: tsv
    tuple val(meta), path(".command.log")               , emit: log
    tuple val(meta), path(".command.sh")                , emit: script

    script:
    output = file("${params.output}").getName()
    """
    # convert IDs parquet to TSV
    parquet_to_tsv.py \\
        --input ${parquet} \\
        --output ${output}.ids.tsv \\
        --header True

    unzstd ${tsv} -f -o ${meta.id}.clusty_input.tsv

    # run clusty
    clusty \\
        --objects-file ${output}.ids.tsv \\
        --algo cd-hit \\
        --id-cols query reference \\
        --distance-col ani \\
        --similarity \\
        --min ani ${params.min_ani} \\
        --min qcov ${params.min_qcov} \\
        --out-representatives \\
        ${meta.id}.clusty_input.tsv \\
        ${output}.clusty.tsv

    zstd ${output}.clusty.tsv
    rm -rf ${output}.ids.tsv ${meta.id}.clusty_input.tsv
    """
}

// Run entry workflow
workflow {

    main:
    // Check if output file already exists
    def output_file = file("${params.output}")
    def input_fna = channel.fromPath(params.input_fna).map { fna ->
            [ [ id: "${fna.getBaseName()}" ], fna ]
        }

    if (!output_file.exists()) {

        // Split input FNA into chunks
        SEQKIT_SPLIT2(
            channel.fromPath(params.input_fna).map { fna ->
                [ [ id: "${fna.getBaseName()}" ], fna ]
            }
        )

        ch_split_fnas = SEQKIT_SPLIT2.out.fnas
            .map { _meta, fnas -> fnas }
            .flatten()
            .map { fna ->
                [ [ id: fna.getBaseName() ], fna ]
            }

        // Create kmer-db database from input FNA
        KMERDB_BUILD(input_fna)

        // Align split FNAs against kmer-db database
        ALIGN(
            ch_split_fnas,
            KMERDB_BUILD.out.kmerdb.join(input_fna).collect()
        )

        // Combine ANI scores
        COMBINEANIS(
            ALIGN.out.parquet.map { _meta, parquet -> [ [ id:'combined'], parquet ] }.groupTuple(sort:'deep')
        )

        if ( (params.cluster == 'mcl') || (params.cluster == 'MCL') ) {
            // Cluster combined ANI scores with MCL
            MCL(
                COMBINEANIS.out.tsv
            )
        } else if ( (params.cluster == 'clusty') || (params.cluster == 'Clusty') ) {
            CLUSTY(
                COMBINEANIS.out.tsv,
                ALIGN.out.ids.first()
            )
        }

    } else {
        println "[UHVDB/nucleotidecluster]: Output file [${params.output}] already exists!"
    }

    // Delete intermediate and Nextflow-specific files
    def remove_tmp = params.remove_tmp
    workflow.onComplete {
        if (output_file.exists()) {
            def work_dir = new File("./work/")
            def nextflow_dir = new File("./.nextflow/")
            def launch_dir = new File(".")

            work_dir.deleteDir()
            nextflow_dir.deleteDir()
            launch_dir.eachFileRecurse { file ->
                if (file.name ==~ /\.nextflow\.log.*/) {
                    file.delete()
                }
            }

            if (remove_tmp) {
                def tmp_dir = new File("./tmp/")
                tmp_dir.deleteDir()
            }
        }
    }
}
