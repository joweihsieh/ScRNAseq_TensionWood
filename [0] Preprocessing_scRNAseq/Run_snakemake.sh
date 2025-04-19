###############
# The snakemake system was previously established here https://github.com/joweihsieh/Woodformation1136_SingleCell
############### two-sample in a pair

### between Bio1 (Tung) and othrers

# Bio1 (Tung) vs Normal Bio5
nohup snakemake --use-conda -c 48 PtrVert1

# Bio1 (Tung) vs Normal Bio2
nohup snakemake --use-conda -c 48 PtrVert2

# Bio1 (Tung) vs Normal Bio3
nohup snakemake --use-conda -c 48 PtrVert3

# Bio1 (Tung) vs Normal Bio4
nohup snakemake --use-conda -c 48 PtrVert4



# Bio1 (Tung) vs Tension Bio1
nohup snakemake --use-conda -c 48 PtrTens1

# Bio1 (Tung) vs Tension Bio2
nohup snakemake --use-conda -c 48 PtrTens2

# Bio1 (Tung) vs Tension Bio3
nohup snakemake --use-conda -c 48 PtrTens3

# Bio1 (Tung) vs Tension Bio4
nohup snakemake --use-conda -c 48 PtrTens4

# Bio1 (Tung) vs Opposite Bio1
nohup snakemake --use-conda -c 48 PtrOp1

# Bio1 (Tung) vs Opposite Bio2
nohup snakemake --use-conda -c 48 PtrOp2

# Bio1 (Tung) vs Opposite Bio3
nohup snakemake --use-conda -c 48 PtrOp3

# Bio1 (Tung) vs Opposite Bio4
nohup snakemake --use-conda -c 48 PtrOp4


### within Opposite

# Opposite Bio1 vs Opposite Bio2
nohup snakemake --use-conda -c 6 PtrOp12_rep

# Opposite Bio1 vs Opposite Bio3
nohup snakemake --use-conda -c 6 PtrOp13_rep

# Opposite Bio1 vs Opposite Bio4
nohup snakemake --use-conda -c 6 PtrOp14_rep

# Opposite Bio2 vs Opposite Bio3
nohup snakemake --use-conda -c 6 PtrOp23_rep

# Opposite Bio2 vs Opposite Bio4
nohup snakemake --use-conda -c 6 PtrOp24_rep

# Opposite Bio3 vs Opposite Bio4
nohup snakemake --use-conda -c 6 PtrOp34_rep

### within Normal

# Normal Bio5 vs Normal Bio2
nohup snakemake --use-conda -c 6 PtrVert12_rep

# Normal Bio5 vs Normal Bio3
nohup snakemake --use-conda -c 6 PtrVert13_rep

# Normal Bio5 vs Normal Bio4
nohup snakemake --use-conda -c 6 PtrVert14_rep

# Normal Bio2 vs Normal Bio3
nohup snakemake --use-conda -c 6 PtrVert23_rep

# Normal Bio2 vs Normal Bio4
nohup snakemake --use-conda -c 6 PtrVert24_rep

# Normal Bio3 vs Normal Bio4
nohup snakemake --use-conda -c 6 PtrVert34_rep


### within Tension

# Tension Bio1 vs Tension Bio2
nohup snakemake --use-conda -c 6 PtrTens12_rep

# Tension Bio1 vs Tension Bio3
nohup snakemake --use-conda -c 6 PtrTens13_rep

# Tension Bio1 vs Tension Bio4
nohup snakemake --use-conda -c 6 PtrTens14_rep

# Tension Bio2 vs Tension Bio3
nohup snakemake --use-conda -c 6 PtrTens23_rep

# Tension Bio2 vs Tension Bio4
nohup snakemake --use-conda -c 6 PtrTens24_rep

# Tension Bio3 vs Tension Bio4
nohup snakemake --use-conda -c 6 PtrTens34_rep

############### all 13 samples

nohup snakemake --use-conda -c 48 PtrTensionWoodALL_13samples
