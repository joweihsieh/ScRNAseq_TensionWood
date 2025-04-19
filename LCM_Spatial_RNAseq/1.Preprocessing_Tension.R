
library(ShortRead)


sample_names <- c("Tension_Bio1_L1", "Tension_Bio1_L2", "Tension_Bio1_L3", 
    "Tension_Bio2_L1", "Tension_Bio2_L2", "Tension_Bio2_L3",
    "Tension_Bio3_L1", "Tension_Bio3_L2", "Tension_Bio3_L3")

#### Separation into multiple samples

whitelist <- readLines("whitelist.txt")

fastq_file <- "LCMSpTrTW_barcoded2_trimmed_multi.fastq.gz"
fastq_data <- readFastq(fastq_file)


for (i in 1:length(whitelist)) {
  current_string <- whitelist[i]
  matching_read_indices <- grep(current_string, id(fastq_data)) 
  matching_reads <- fastq_data[matching_read_indices]
  
  if (length(matching_reads) > 0) {
      print(sample_names[i])
      output_file <- paste0(sample_names[i], ".fastq.gz")
      
      writeFastq(matching_reads, output_file, compress = TRUE)
      
      message(paste("Saved matching reads for", current_string, "to", output_file))
    }
  }
}


#### mapping

for (i in 1:length(sample_names)) {
  fastq_file <- paste0(sample_names[i], ".fastq.gz")
  
  if (!file.exists(fastq_file)) {
    warning(paste("File", fastq_file, "does not exist. Skipping this sample."))
    next
  }
  
  output_prefix <- paste0("mapped_", sample_names[i], "_")  
  star_command <- paste(
    "STAR --runThreadN 10",
    "--genomeDir /home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/genome/STAR_GenomeIndex_150",
    "--readFilesIn", fastq_file,
    "--readFilesCommand zcat",
    "--outSAMtype BAM SortedByCoordinate",
    "--outFileNamePrefix", output_prefix
  )
  
  tryCatch({
    system(star_command)
    message(paste("Mapping for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("Mapping for", sample_names[i], "failed. Error message:", e$message))
  })
}

############################## conda activate subread for processing mapped bam files 
library(ShortRead)


sample_names <- c("Tension_Bio1_L1", "Tension_Bio1_L2", "Tension_Bio1_L3", 
    "Tension_Bio2_L1", "Tension_Bio2_L2", "Tension_Bio2_L3",
    "Tension_Bio3_L1", "Tension_Bio3_L2", "Tension_Bio3_L3")



# FeatureCounts
for (i in 1:length(sample_names)) {
  output_prefix <- paste0("mapped_", sample_names[i], "_")
  bam_file <- paste0(output_prefix, "Aligned.sortedByCoord.out.bam")
  feature_counts_bam <- paste0(bam_file, ".featureCounts.bam")
  
  if (!file.exists(bam_file)) {
    warning(paste("BAM file for", sample_names[i], "does not exist. Skipping this sample."))
    next
  }
  
  # add gene tag into BAM
  featureCounts_command <- paste(
    "featureCounts -a /home/woodydrylab/FileShare/LCM_spatial_RNAseq/TNNR24PRR02378_2xRNAseq_lib/01.Raw_read/genome/Ptrichocarpa_533_v4.1.gene.gtf",
    "-o", paste0("mapped_", sample_names[i], "_gene_counts.txt"),
    "-R BAM",
    bam_file
  )
  
  tryCatch({
    system(featureCounts_command)
    message(paste("Feature counting for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("Feature counting for", sample_names[i], "failed. Error message:", e$message))
    next
  })
  
  if (!file.exists(feature_counts_bam)) {
    warning(paste("FeatureCounts BAM file for", sample_names[i], "not found. Skipping sorting."))
    next
  }
  
  # sorting
  sorted_bam_file <- paste0(output_prefix, "Aligned.sortedByCoord.out.featureCounts.sorted.bam")
  sorting_command <- paste(
    "samtools sort -o", sorted_bam_file, feature_counts_bam
  )
  
  tryCatch({
    system(sorting_command)
    message(paste("Sorting for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("Sorting for", sample_names[i], "failed. Error message:", e$message))
    next
  })
  
  # indexing
  tryCatch({
    system(paste("samtools index", sorted_bam_file))
    message(paste("Indexing for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("Indexing for", sample_names[i], "failed. Error message:", e$message))
  })
}




# UMI-Tools Grouping
for (i in 1:length(sample_names)) {
  output_prefix <- paste0("mapped_", sample_names[i], "_")
  bam_file <- paste0(output_prefix, "Aligned.sortedByCoord.out.featureCounts.sorted.bam")
  
  if (!file.exists(bam_file)) {
    warning(paste("BAM file for", sample_names[i], "does not exist. Skipping UMI grouping for this sample."))
    next
  }
  
  grouped_tsv_output <- paste0(output_prefix, "grouped_output_featureCounts.tsv")
  grouped_bam_output <- paste0(output_prefix, "Aligned.sortedByCoord.out.featureCounts.gp.bam")
  log_output <- paste0(output_prefix, "grouping_featureCounts.log")
  
  tryCatch({
    system(paste("samtools index", bam_file))
    message(paste("Indexing for", bam_file, "completed successfully."))
  }, error = function(e) {
    warning(paste("Indexing for", bam_file, "failed. Error message:", e$message))
    next
  })
  
  umi_tools_command <- paste(
    "umi_tools group",
    "--extract-umi-method=read_id",
    "-I", bam_file,
    "--group-out", grouped_tsv_output,
    "--umi-separator=_",
    "--output-bam",
    "--stdout", grouped_bam_output,
    "--log", log_output
  )

  tryCatch({
    system(umi_tools_command)
    message(paste("UMI grouping for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("UMI grouping for", sample_names[i], "failed. Error message:", e$message))
    next
  })
}


# UMI-Tools Deduplication
for (i in 1:length(sample_names)) {
  output_prefix <- paste0("mapped_", sample_names[i], "_")
  grouped_bam_input <- paste0(output_prefix, "Aligned.sortedByCoord.out.featureCounts.gp.bam")
  deduplicated_bam_output <- paste0(sample_names[i], "_deduplicated.bam")
  dedup_log_output <- paste0(sample_names[i], "_deduplication.log")
  
  if (!file.exists(grouped_bam_input)) {
    warning(paste("Grouped BAM file for", sample_names[i], "not found. Skipping deduplication for this sample."))
    next
  }

  tryCatch({
    system(paste("samtools index", grouped_bam_input))
    message(paste("Indexing for", grouped_bam_input, "completed successfully."))
  }, error = function(e) {
    warning(paste("Indexing for", grouped_bam_input, "failed. Error message:", e$message))
    next
  })

  dedup_command <- paste(
    "umi_tools dedup",
    "--stdin", grouped_bam_input,
    "--stdout", deduplicated_bam_output,
    "--log", dedup_log_output,
    "--extract-umi-method=read_id"
  )

  tryCatch({
    system(dedup_command)
    message(paste("Deduplication for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("Deduplication for", sample_names[i], "failed. Error message:", e$message))
  })
}




# UMI-Tools Count
for (i in 1:length(sample_names)) {
  output_prefix <- paste0(sample_names[i], "_")
  deduped_bam_input <- paste0(output_prefix, "deduplicated.bam")
  umi_counts_output <- paste0(output_prefix, "counts.tsv")
  umi_counts_log_output <- paste0(output_prefix, "counts.log")

  if (!file.exists(deduped_bam_input)) {
    warning(paste("deduped BAM file for", sample_names[i], "not found. Skipping UMI counting for this sample."))
    next
  }

  tryCatch({
    system(paste("samtools index", deduped_bam_input))
    message(paste("Indexing for", deduped_bam_input, "completed successfully."))
  }, error = function(e) {
    warning(paste("Indexing for", deduped_bam_input, "failed. Error message:", e$message))
    next
  })

  count_command <- paste(
    "umi_tools count",
    "--extract-umi-method=read_id",
    "--gene-tag=XT",
    "--per-gene",
    "--stdin", deduped_bam_input,
    "--log", umi_counts_log_output,
    "--stdout", umi_counts_output
  )

  tryCatch({
    system(count_command)
    message(paste("UMI counting for", sample_names[i], "completed successfully."))
  }, error = function(e) {
    warning(paste("UMI counting for", sample_names[i], "failed. Error message:", e$message))
  })
}

