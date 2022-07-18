nextflow.enable.dsl=2

workflow {
    Channel
        .fromSRA('PRJNA475130')
        .view()
        .set{reads}

	
    process getreads {
      input: 
        set sample_id, file(reads_file) from reads

      output: 
        file("fastqc_${sample_id}_logs") into fastqc_ch

      script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads_file}
        """
}
}