#### QIAGEN QIAseq UPX 3' Transcriptome kit
#https://www.biostars.org/p/372926/
#https://www.biostars.org/p/372926/
#### The UMI and cell index sequences are extracted from the R2 reads. A specific read structure is expected from the R2 reads. 
#because we use a custom sequencing primer that binds directly before the cell index. The expected R2 read structure from a MiSeq® or HiSeq® instrument is:
#<AAGCAGTGGTATCAACGCAGAGTAC><cell_index><UMI><ACG><poly-T>
### https://umi-tools.readthedocs.io/en/latest/reference/extract.html

### extract - normal
umi_tools extract --extract-method=regex \
	--stdin="LCMSpTrNW_S1_L001_R2_001.fastq.gz" \
    --read2-in="LCMSpTrNW_S1_L001_R1_001.fastq.gz" \
	--bc-pattern="(?P<discard_1>AAGCAGTGGTATCAACGCAGAGTAC{s<=2})(?P<cell_index>.{10})(?P<umi_1>.{12}).*" \
    --read2-stdout \
    --whitelist="whitelist.txt" \
    --filter-cell-barcode \
    --log="umi_extract_log2.txt" \
| gzip > LCMSpTrNW_barcoded2.fastq.gz




cutadapt -u 2 -a file:customadapt.fasta --times 3 -m 25 -o LCMSpTrNW_barcoded2_trimmed_multi.fastq.gz LCMSpTrNW_barcoded2.fastq.gz > trim_with_adapter_multi.log 2>&1
fastqc LCMSpTrNW_barcoded2_trimmed_multi.fastq.gz

### extract - tension

umi_tools extract --extract-method=regex \
	--stdin="LCMSpTrTW_S2_L001_R2_001.fastq.gz" \
    --read2-in="LCMSpTrTW_S2_L001_R1_001.fastq.gz" \
	--bc-pattern="(?P<discard_1>AAGCAGTGGTATCAACGCAGAGTAC{s<=2})(?P<cell_index>.{10})(?P<umi_1>.{12}).*" \
    --read2-stdout \
    --whitelist="whitelist.txt" \
    --filter-cell-barcode \
    --log="umi_extract_log.txt" \
| gzip > LCMSpTrTW_barcoded2.fastq.gz


cutadapt -u 2 -a file:customadapt.fasta --times 3 -m 25 -o LCMSpTrTW_barcoded2_trimmed_multi.fastq.gz LCMSpTrTW_barcoded2.fastq.gz > trim_with_adapter_multi.log 2>&1
fastqc LCMSpTrTW_barcoded2_trimmed_multi.fastq.gz
