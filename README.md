# Insights into *Xanthomonas arboricola* pv. *juglandis* in California: population genomics, copper resistance, and type III effectors

## Read Processing
FastP v.0.23.4 (Chen et al., 2018) was used for read processing.
``` bash
fastp -i $READ1 -I $READ2 \
  -o $FILTERED/$SID.F.1.fq.gz -O $FILTERED/$SID.F.2.fq.gz \
  -q 30 \
  --detect_adapter_for_pe \
  -h $REPORTS/$SID.html
```
## Genome Assembly and Annotation
We followed these steps for de novo genome assembly and annotation:
1. De novo assembly using Platanus_B v1.3.2 (Kajitani et al. 2020).
2. Generate assembly stats via QUAST v5.1.0rc1 (Gurevich et al. 2013).
3. Scaffold contigs using Ragout v.2.3 (Kolmogorov et al. 2014).
4. Polish scaffolds with Pilon v.1.24(Walker et al. 2014).
5. Assess completeness and contamination with CheckM v.1.2.1(Parks et al. 2015) and BUSCO v6.0.0 (Tegenfeldt et al. 2025).
6. Annotate with NCBI PGAP v.2024-07-18.build7555(Tatusova et al. 2016). 
7. Assemble plasmids using MOB-suite v3.1.9 (Robertson and Nash 2018).

### Platanus_B
``` bash
platanus_b assemble \
-o "$STRAIN" \
-f "$READ1" "$READ2" \
-t 16 -m 64
```

### QUAST 
``` bash
quast.py "${STRAIN}_contig.fa" -o "$STATS/Contigs/$STRAIN"
```

### Ragout
```bash
ragout -o "$STRAIN" -s sibelia --refine -t 16 "${STRAIN}.rcp"
```

### Pilon
```bash
# BWA-mem2 v.2.2.1
bwa-mem2 index -p "$INDEX/${STRAIN}" "${STRAIN}_scaffolds.fasta"
bwa-mem2 mem -t 16 "$INDEX/${STRAIN}" \
"$READ1" "$READ2" | samtools view -b | samtools sort -@ 16 -o "$BAMS/${STRAIN}.bam"

# Samtools v.1.19.2
samtools index "$BAMS/${STRAIN}.bam"

# Pilon 
java -jar -Xmx64G $PILONJAR --genome "${STRAIN}_scaffolds.fasta" \
--bam "$BAMS/${STRAIN}.bam" \
--output "${STRAIN}_Pilon" \
--outdir "$RESULTS/$STRAIN"
```

### CheckM
```bash
# CheckM v.1.2.1
checkm taxonomy_wf -t 16 genus Xanthomonas -x fasta "$CONTIGS" "$OUTDIR/CheckM_Tax_Contigs"
```

### PGAP
```bash
# PGAP v.2024-07-18.build7555
"$PGAPPATH/pgap.py" -o "$STRAIN" --no-self-update -n --ignore-all-errors --taxcheck "$YAMLDIR/$STRAIN.generic.yaml"
```

### MOB-Recon
```bash
# MOB_recon v3.1.9
mob_recon -n 8 -s "$STRAIN" -i "$GENOMEDIR/$STRAIN.fna" -o "$DIRECTORY"
```


## Pangenome 

### PPanGGoLiN 

```bash
# PPanGGoLiN v.2.2.0
# run all analyses
ppanggolin all --anno "$POP.tsv" -c 16
# generate the MSA of the accessory genome
ppanggolin msa -p pangenome.h5 --partition accessory -o CA_XAJ_115_Accessory_MSA
```

### Phylogenetic tree of accessory genome

```bash
# IQtree v.2.2.2.6 
# ModelFinder Accessory 
iqtree2 -s "$MSADIR/accessory_genome_alignment.aln" -mem 32G -T AUTO -m MFP --prefix "$ACC/MFout"
# Phylogenetic reconstruction
iqtree2 -s "$MSADIR/accessory_genome_alignment.aln" -mem 32G -T AUTO -m JTT+F+I+R2 -B 1000 --prefix "$ACC/IQout"
```

# Core-Genome Phylogenetic Analyses
We followed these steps for the recombination-corrected core-genome phylogenetic and population structure analyses:

1. Core-genome alignment using Parsnp v.2.0.6
2. Detect and remove recombination using ClonalFrameML v.1.13
3. Phylogenetic reconstruction using IQtree v.2.2.2.6
4. Population structure analysis using RhierBAPS v.1.1.4  

## Parsnp
```bash
parsnp -p 8 -u --vcf -d $GENOMES -r $XAJREF -o $DIRECTORY
```

## ClonalFrameML
```bash
ClonalFrameML "$INDIR/parsnp.tree" "$INDIR/parsnp.xmfa" "$OUTDIR/CA_XAJ_115_v2_CFML" -emsim 100 -output_filtered true -xmfa_file true -show_progress true
```

## IQtree

```bash
# ModelFinder Accessory 
iqtree2 -s "$ALIGNMENT" -mem 32G -T AUTO -m MFP --prefix "$MFDIRECTORY/MFout"

# Substitution models used
# CA_XAJ: HKY+F+I 
# CA_XAJ and 24 NCBI RefSeq: TVM+F+G4
# CA_XAJ, 24 NCBI RefSeq, and Xhg: TVM+F+G4 

# Phylogenetic reconstruction
# with outgroup
iqtree2 -s "$ALIGNMENT" -mem 32G -T AUTO -m "$MODEL" -B 1000 -o "$OUTGROUP" --prefix "$IQDIR/IQout"
# without outgroup
iqtree2 -s "$ALIGNMENT" -mem 32G -T AUTO -m "$MODEL" -B 1000 --prefix "$IQDIR/IQout"
```

## RhierBAPS
``` r
# RhierBAPS v.1.1.4 
library(rhierbaps)

# Resource: https://github.com/gtonkinhill/rhierbaps

# Load fasta and generate matrix
Fasta = "CA_XAJ_CFML.filtered.fasta"
snp.matrix = load_fasta(Fasta)

# Running hierBAPS
hb.results = hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE, n.cores = 8)

# save the partition file
write.csv(hb.results$partition.df, file = file.path("Outputs/CA_XAJ_115_hierbaps_partition_npop20.csv"), col.names = TRUE, row.names = FALSE)
```

# TXSS Effectors

## Type III effectors

### tBLASTn
```bash
# ncbi-blast v.2.16.0+ 
makeblastdb -in $GENOMEDIR/$STRAIN.fna -input_type fasta -dbtype nucl

# alignment output
tblastn -query $T3SSE \
-db $GENOMEDIR/$STRAIN.fna \
-out $BLASTDIR/$STRAIN/$STRAIN.output \
-evalue 1e-10 \
-num_threads 8 \
-max_target_seqs 5 \
-max_hsps 5

# tab delimited output
tblastn -query $T3SSE \
-db $GENOMEDIR/$STRAIN.fna \
-out $BLASTDIR/$STRAIN/$STRAIN.txt \
-evalue 1e-10 \
-num_threads 8 \
-max_target_seqs 5 \
-max_hsps 5 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
```

## DeepSecE 
``` bash
# DeepSecE v.0.1.2
python3 predict_genome.py \
--fasta_path "$PROTS/$STRAIN.faa" \
--model_location "$MODEL" \
--data_dir "$TXSS" \
--out_dir "$OUTDIR/$STRAIN" \
--save_attn \
```

# CAZymes

```bash 
# db_can3 run_dbcan v.4.1.4
run_dbcan "$PROT/$STRAIN.faa" \
protein \
-c "$GFFDIR/$STRAIN.gff" \
--out_dir "$CAZYDIR/$STRAIN" \
--out_pre "${STRAIN}." \
-dt 16 \
--db_dir "$DBDIR" 
```
