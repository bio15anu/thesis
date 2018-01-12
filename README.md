# Thesis Pipeline
All files pertaining to masters thesis are contained within this repository

See below for project pipeline commands


----------------------------
### RELEVANT DIRECTORIES ###

```
~/bin/ # software binaries
~/scripts # custom scripts are stored here for easy access
~/Documents/Databases # where blastp, pfam, merops etc. databases are stored
~/Documents/RawData # containing the original full transcriptome assembly and sequencing reads
~/Documents/RawData/Reads # sequencing reads that mapped to Emuscae draft genome
~/Documents/RawData/stats # Trinity_full.fasta assembly statistics
~/Documents/Data # containing all data and analyses derived during the present thesis
~/Documents/Data/Alignments # all data pertaining to protein alignments and related analyses
~/Documents/Data/Annotations # functional annotations of the final transcriptome (pertaining to Trinotate)
~/Documents/Data/Assembly # data pertaining to the final (filtered) transcriptome assembly
~/Documents/Data/Assembly/Matrices # abundance estimates on the final transcriptome assembly
~/Documents/Data/DiffExpr # all data pertaining to differential expression analyses
~/Documents/Data/DiffExpr/GSEA # all data pertaining to gene set enrichment analyses
```

---------------------------------
### PREPARING LOCAL DATABASES ###

#### Starting directory
    ~/Documents/Data/Annotations/>

#### build the SQL database for annotations (also retrieve current Pfam and uniprot_sprot.pep database)
    ~/bin/Trinotate-3.0.2/admin/Build_Trinotate_Boilerplate_SQLite_db.pl Emuscae
    mv Pfam-A.hmm.gz ~/Documents/Databases/pfam/Pfam-A-hmm.gz
    mv uniprot_sprot.pep ~/Documents/Databases/uniprot_sprot/uniprot_sprot.pep

#### decompress pfam and prepare (in pfam directory)
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm

#### build sprot databases (in uniprot_sprot directory)
    makeblastdb -in uniprot_sprot.pep -dbtype prot
#### parse uniprot_sprot for fungi and diptera
```
~/scripts/sprot_fetchID.py uniprot_sprot.dat Diptera.txt Diptera
~/scripts/sprot_fetchID.py uniprot_sprot.dat Fungi.txt Fungi
~/scripts/sprot_extractID.py Diptera.txt uniprot_sprot.pep Diptera.uniprot_sprot.pep
~/scripts/sprot_extractID.py Fungi.txt uniprot_sprot.pep Fungi.uniprot_sprot.pep
cat Diptera.uniprot_sprot.pep Fungi.uniprot_sprot.pep > combine.uniprot_sprot.pep
makeblastdb -in combine.uniprot_sprot.pep -dbtype prot
```

#### build MEROPS database (in merops directory)
    wget "ftp://ftp.ebi.ac.uk/pub/databases/merops/current_release/pepunit.lib"
    makeblastdb -in pepunit.lib -dbtype prot


-----------------------------------------------
### PREFILTERING THE TRANSCRIPTOME ASSEMBLY ###

#### Starting directory
    ~/Documents/RawData/>

#### generate RSEM counts for Trinity_full.fasta (I originally ran a test of this on one sample with the --prep_reference option included, but since that had been done already I did not include the option again in the command here)
```
cd Reads
ls | grep -v Archive | while read file; do id=`echo $file | sed 's/\_mapped//g'`; ~/bin/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts ../Trinity_full.fasta --seqType fq --left ${file}/${id}_mapped_JGIDraftgenome.1 --right ${file}/${id}_mapped_JGIDraftgenome.2 --est_method RSEM --aln_method bowtie --trinity_mode --output_dir ${file}/${id}/; done
```

#### copy all RSEM.isoforms.results and RSEM.genes.results files for each sample to single directory
    mkdir ~/Documents/RawData/Reads/isoforms
    mkdir ~/Documents/RawData/Reads/genes

#### generate concatenated matrix files for "isoforms" level
```
cd isoforms
~/bin/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM 8_6.isoforms.results 8_7.isoforms.results 9_1_1.isoforms.results 9_6.isoforms.results MdEm1.isoforms.results MdEm2.isoforms.results MdEm3.isoforms.results Em1glen.isoforms.results Em2glen.isoforms.results Em3glen.isoforms.results sp_5.isoforms.results sp_6_1.isoforms.results sp_6_2.isoforms.results
```

#### generate concatenated matrix files for "genes" level
```
cd ../genes
~/bin/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM 8_6.genes.results 8_7.genes.results 9_1_1.genes.results 9_6.genes.results MdEm1.genes.results MdEm2.genes.results MdEm3.genes.results Em1glen.genes.results Em2glen.genes.results Em3glen.genes.results sp_5.genes.results sp_6_1.genes.results sp_6_2.genes.results
```

#### TransDecoder Filtering (150 bp / 50 aa)
```
cd ../../
TransDecoder.LongOrfs -t Trinity_full.fasta -m50
hmmscan --cpu 8 --domtblout Trinity_full.fasta.transdecoder_dir/pfam.domtblout ~/Documents/Databases/pfam/Pfam-A.hmm Trinity_full.fasta.transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t Trinity_full.fasta --retain_pfam_hits Trinity_full.fasta.transdecoder_dir/pfam.domtblout --single_best_orf
```

#### filter Trinity_full.fasta to contain only sequences that encode proteins in Trinity_trans.fasta file
    ~/scripts/assembly_prefilter_with_Transpep.py Trinity_full.fasta.transdecoder_dir/Trinity_full.fasta.transdecoder.pep Trinity_full.fasta Trinity_trans.fasta

#### generate gene_trans_map for Trinity_trans.fasta
    ~/bin/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity_trans.fasta > Trinity_trans.fasta.gene_trans_map

#### generate Trinity_isoH.fasta using the full isoform TMM count matrix to filter the Trinity_trans.fasta file
    ~/bin/trinityrnaseq-Trinity-v2.4.0/util/filter_low_expr_transcripts.pl --matrix matrix.TMM.EXPR.matrix --transcripts Trinity_trans.fasta --highest_iso_only --gene_to_trans_map Trinity_trans.fasta.gene_trans_map > Trinity_isoH.fasta

#### filter the full count matrices to contain only most highly expressed isoforms that encode proteins
```
ls Reads/isoforms | grep matrix | grep -v info | while read file; do cat Reads/isoforms/${file} | head -1 > ${file}; cat Trinity_isoH.fasta | grep ">" | tr -d ">" | cut -d " " -f1 | while read line; do cat Reads/isoforms/${file} | grep $line ; done >> ${file}; done
```

#### use CD-HIT-EST to filter Trinity_isoH.fasta file at 97% identity to build final assembly Trinity_isoH_cdhit10.fasta
    cd-hit-est -i Trinity_isoH.fasta -o ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta -c 0.97 -n 10

#### filter Trinity_full.fasta.transdecoder.pep to contain only sequences from final assembly Trinity_isoH_cdhit10.fasta
```
cat Trinity_isoH_cdhit10.fasta | grep ">" | cut -d " " -f1 | tr -d ">" | while read line; do cat Trinity_full.fasta.transdecoder_dir/Trinity_full.fasta.transdecoder.pep | grep -A1 $line ; done > ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep
```

#### generate SuperTranscript from full assembly (for later antiSMASH analysis)
    ~/bin/trinityrnaseq-Trinity-v2.5.1/Analysis/SuperTranscripts/Trinity_gene_splice_modeler.py --trinity_fasta Trinity_full.fasta --out_prefix Trinity_full.SuperTranscript


------------------------------------------------------------------------------
### FINAL (PREFILTERED) TRANSCRIPTOME PREPARATION AND ABUNDANCE ESTIMATION ###

#### Starting directory
    ~/Documents/Data/Assembly/>

#### abundance estimates on the final transcriptome (remapping reads)
```
cd ~/Documents/RawData/Reads
ls | grep -v Archive | grep -v isoforms | grep -v genes | while read file; do id=`echo $file | sed 's/\_mapped//g'`; ~/bin/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta --seqType fq --left ${file}/${id}_mapped_JGIDraftgenome.1 --right ${file}/${id}_mapped_JGIDraftgenome.2 --est_method RSEM --aln_method bowtie --trinity_mode --output_dir ${file}/isoH_cdhit10/; done
```

#### copy all RSEM.isoforms.results and RSEM.genes.results files for each sample to single directory
    mkdir ~/Documents/Data/Assembly/Matrices/isoforms
    mkdir ~/Documents/Data/Assembly/Matrices/genes

#### generate concatenated matrix files for "isoforms" level
```
cd ~/Documents/Data/Assembly/Matrices/isoforms
~/bin/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM 8_6.isoforms.results 8_7.isoforms.results 9_1_1.isoforms.results 9_6.isoforms.results MdEm1.isoforms.results MdEm2.isoforms.results MdEm3.isoforms.results Em1glen.isoforms.results Em2glen.isoforms.results Em3glen.isoforms.results sp_5.isoforms.results sp_6_1.isoforms.results sp_6_2.isoforms.results
```

#### generate concatenated matrix files for "genes" level
```
cd ../genes
~/bin/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM 8_6.genes.results 8_7.genes.results 9_1_1.genes.results 9_6.genes.results MdEm1.genes.results MdEm2.genes.results MdEm3.genes.results Em1glen.genes.results Em2glen.genes.results Em3glen.genes.results sp_5.genes.results sp_6_1.genes.results sp_6_2.genes.results
```

#### generate gene_trans_map for Trinity_isoH_cdhit10.fasta (for later analyses)
    cd ../../
    ~/bin/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity_isoH_cdhit10.fasta > Trinity_isoH_cdhit10.fasta.gene_trans_map

#### generate seq lengths file (for later analyses)
    ~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/fasta_seq_length.pl Trinity_isoH_cdhit10.fasta > Trinity_isoH_cdhit10.fasta.seq_lens

#### generate gene lengths file (for later analyses)
    ~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/TPM_weighted_gene_length.py --gene_trans_map Trinity_isoH_cdhit10.fasta.gene_trans_map --trans_lengths Trinity_isoH_cdhit10.fasta.seq_lens --TPM_matrix Matrices/genes/matrix.TMM.EXPR.matrix > Trinity_isoH_cdhit10.gene_lengths.txt


------------------------------------------------------------
### FUNCTIONAL ANNOTATIONS FOR PREFILTERED TRANSCRIPTOME ###

#### Starting directory
    ~/Documents/Data/Annotations/>

#### blast annotations - custom database
#### blastx (from geneious server)
    nohup blastx -query Trinity_isoH_chdit10.fasta -db combine.uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > Trinity_isoH_cdhit10.blastx.outfmt6
#### blastp (from local uniprot_sprot directory. See: PREPARING LOCAL DATABASES)
    blastp -query ~/Documents/ThesisData/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep -db combine.uniprot_sprot.pep -num_threads 8 -max_target_seqs 1 -outfmt 6 > ~/Documents/ThesisData/Annotations/Trinity_isoH_cdhit10.blastp.outfmt6

#### blast annotations - MEROPS database
#### blastx (from geneious server)
    nohup blastx -query Trinity_isoH_chdit10.fasta -db pepunit.lib -num_threads 8 -max_target_seqs 1 -outfmt 6 > MEROPS.blastx.outfmt6
#### blastp (from local MEROPS directory. See: PREPARING LOCAL DATABASES)
    blastp -query ~/Documents/ThesisData/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep -db pepunit.lib -num_threads 8 -max_target_seqs 1 -outfmt 6 > ~/Documents/ThesisData/Annotations/MEROPS.blastp.outfmt6

#### Pfam annotations
    hmmscan --cpu 8 --domtblout TrinotatePFAM.out ~/Documents/Databases/pfam/Pfam-A.hmm ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep > pfam.log

#### signalP annotations - transdecoder IDs are too long for signalP so custom scripts are used to convert them
```
cd SignalP
~/scripts/sigP_mend_Transpep_for_sigP.py ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep signalp.transdecoder.pep
signalp -f short -n signalp.temp signalp.transdecoder.pep
~/scripts/sigP_mend_temp.py signalp.temp ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep ../signalp.out
```

#### tmhmm annotations
    tmhmm --short < ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.transdecoder.pep > tmhmm.out

#### rnammer annotations
    ~/bin/Trinotate-3.0.2/util/rnammer_support/RnammerTranscriptome.pl --transcriptome ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta --path_to_rnammer ~/bin/rnammer-1.2/rnammer

#### compile database
```
~/bin/Trinotate-3.0.2/util/support_scripts/get_Trinity_gene_to_trans_map.pl ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta > ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.gene_trans_map
Trinotate Emuscae.sqlite init --gene_trans_map ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.gene_trans_map --transcript_fasta ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta --transdecoder_pep ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep
Trinotate Emuscae.sqlite LOAD_swissprot_blastx Trinity_isoH_cdhit10.blastx.outfmt6
Trinotate Emuscae.sqlite LOAD_swissprot_blastp Trinity_isoH_cdhit10.blastp.outfmt6
Trinotate Emuscae.sqlite LOAD_custom_blast --outfmt6 MEROPS.blastx.outfmt6 --prog blastx --dbtype MEROPS
Trinotate Emuscae.sqlite LOAD_custom_blast --outfmt6 MEROPS.blastp.outfmt6 --prog blastp --dbtype MEROPS
Trinotate Emuscae.sqlite LOAD_pfam TrinotatePFAM.out
Trinotate Emuscae.sqlite LOAD_signalp signalp.out
Trinotate Emuscae.sqlite LOAD_tmhmm tmhmm.out
Trinotate Emuscae.sqlite LOAD_rnammer Trinity_isoH_cdhit10.fasta.rnammer.gff
Trinotate Emuscae.sqlite report > Emuscae_annotation_report.xls
```

#### collating metadata from Emuscae_annotation_report.xls
##### number of putative genes
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$5!="." {print $1}' | sort | uniq | wc -l
##### putative genes with blastx annotations - swissprot
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$3!="." {print $1}' | sort | uniq | wc -l
##### putative genes with blastp annotations - swissprot
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$7!="." {print $1}' | sort | uniq | wc -l
##### putative genes with blastx annotations - MEROPS
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$8!="." {print $1}' | sort | uniq | wc -l
##### putative genes with blastp annotations - MEROPS
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$9!="." {print $1}' | sort | uniq | wc -l
##### putative genes with pfam domains
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$10!="." {print $1}' | sort | uniq | wc -l
##### unique pfam domains
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$10!="." {print $10}' | tr "\`" "\n" | cut -d \^ -f1 | sort | uniq | wc -l
##### putative genes with GO annotations
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '{if ($15!="." || $16!=".") {print $1}}' | sort | uniq | wc -l
##### putative genes with KEGG annotations
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$14!="." {print $1}' | sort | uniq | wc -l
##### putative genes with TMHMM annotations
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$12!="." {print $1}' | sort | uniq | wc -l
##### putative genes with SignalP annotations
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$11!="." {print $1}' | sort | uniq | wc -l
##### number of SSPs: 50<aa<300, >3% cysteine, no TMHMM signal (unfiltered for blast and pfam annotations)
```
mkdir Secretome
cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '{if($11!="." && $12==".") {print $5}}' | while read line; do cat ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep | grep -A1 $line; done | cut -d " " -f1 > Secretome/unfiltered.fasta
~/scripts/secretome_parseSSP.py Secretome/unfiltered.fasta 0.03 300 Secretome/unfiltered
cat Secretome/unfiltered.0.03.300aa.fasta | grep -c ">"
```
##### number of SSPs: 50<aa<300, >3% cysteine, no TMHMM signal (filtered for blast and pfam annotations)
```
cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '{if($3=="." && $7=="." && $10=="." && $11!="." && $12==".") {print $5}}' | while read line; do cat ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep | grep -A1 $line; done | cut -d " " -f1 > Secretome/filtered.fasta
~/scripts/secretome_parseSSP.py Secretome/filtered.fasta 0.03 300 Secretome/filtered
cat Secretome/filtered.0.03.300aa.fasta | grep -c ">"
```
##### putative genes with rnammer annotation
    cat Emuscae_annotation_report.xls | tail -n +2 | awk -F "\t" '$4!="." {print $1}' | sort | uniq | wc -l

#### build GO annotation file
    ~/bin/Trinotate-3.0.2/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Emuscae_annotation_report.xls -G --include_ancestral_terms > go_annotations.txt


---------------------------
### ASSEMBLY STATISTICS ###

#### Starting directory
    ~/Documents/RawData/Reads>

#### RNA-seq read representation (geneious server)
```
mkdir stats
bowtie2-build Trinity_full.fasta Trinity_full.fasta
ls | grep -Ev "Archive|isoforms|genes|bowtie2|stats" | while read file; do id=`echo $file | sed 's/\_mapped//g'`; bowtie2 -p 10 -q --no-unal -x bowtie2/Trinity_full.fasta -1 ${file}/${id}_mapped_JGIDraftgenome.1 -2 ${file}/${id}_mapped_JGIDraftgenome.2 | samtools view -@10 -Sb -o bowtie2.bam 2>$1 | tee stats/${id}_stats.txt ; done
```

#### assembly general statistics
```
cd ~/Documents/RawData
mkdir stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl Trinity_full.fasta > stats/Trinity_full.stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl Reads/genes/matrix.TPM.not_cross_norm | tee stats/matrix.TPM.not_cross_norm.counts_by_min_TPM
cd stats
Rscript ~/scripts/stats_counts_by_min_TPM.R 50000 matrix.TPM.not_cross_norm.counts_by_min_TPM

cd ~/Documents/Data/Assembly
mkdir stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl Trinity_isoH_cdhit10.fasta > stats/Trinity_isoH_cdhit10.stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl Matrices/genes/matrix.TPM.not_cross_norm | tee stats/matrix.TPM.not_cross_norm.counts_by_min_TPM
cd stats
Rscript ~/scripts/stats_counts_by_min_TPM.R 25000 matrix.TPM.not_cross_norm.counts_by_min_TPM
```

#### assembly Ex90N50 statistics
```
cd ~/Documents/RawData/stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/contig_ExN50_statistic.pl ../Reads/isoforms/matrix.TMM.EXPR.matrix ../Trinity_full.fasta | tee ExN50.stats
Rscript ~/scripts/stats_ExN50.R ExN50.stats

cd ~/Documents/Data/Assembly/stats
~/bin/trinityrnaseq-Trinity-v2.4.0/util/misc/contig_ExN50_statistic.pl ../Matrices/isoforms/matrix.TMM.EXPR.matrix ../Trinity_isoH_cdhit10.fasta | tee ExN50.stats
Rscript ~/scripts/stats_ExN50.R ExN50.stats
```

#### BUSCO completeness

!!! REMEMBER BUSCO !!!


--------------------------------
### DOWNSTREAM QC STATISTICS ###

#### Starting directory
    ~/Documents/Data/DiffExpr/stats>

#### compare replicates within test conditions (prior to DE analysis)
    ~/bin/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix ~/Documents/Data/Assembly/Matrices/genes/matrix.counts.matrix --samples ../samples.txt --CPM --log2 --min_rowSums 10 --compare_replicates

#### compare replicates across test conditions (prior to DE analysis) - Pearson correlation matrix
    ~/bin/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix ~/Documents/Data/Assembly/Matrices/genes/matrix.counts.matrix --samples ../samples.txt --CPM --log2 --min_rowSums 10 --sample_cor_matrix

#### compare replicates across test conditions (prior to DE analysis) - Principal Component Analysis
    ~/bin/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix ~/Documents/Data/Assembly/Matrices/genes/matrix.counts.matrix --samples ../samples.txt --CPM --log2 --min_rowSums 10 --center_rows --prin_comp 3


----------------------------------------
### DIFFERENTIAL EXPRESSION ANALYSES ###

#### Starting directory
    ~/Documents/Data/DiffExpr

#### run differential expression
```
~/bin/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix ~/Documents/Data/Assembly/Matrices/genes/matrix.counts.matrix --method DESeq2 --samples_file samples.txt --output P0.001
cd P0.001
```

#### analysis of differential expression results & run GOseq - pvalue 0.001 logFC 2-fold
    ~/bin/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ~/Documents/Data/Assembly/Matrices/genes/matrix.TMM.EXPR.matrix --samples ../samples.txt --examine_GO_enrichment --GO_annots ~/Documents/Data/Annotations/go_annotations.txt --gene_lengths ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.gene_lengths.txt

#### sort into directories for each pairwise comparison (coni_vs_glen etc.)
```
ls | cut -d "." -f4 | grep _vs_ | sort | uniq | while read file; do mkdir ${file}; mv *.${file}.* ${file}; done
mkdir unique vs_all
```

#### sort each eg. coni_vs_glen directory into eg. coni-UP and glen-UP
    ls | grep _vs_ | grep -v pdf | while read file; do ls ${file} | grep UP | cut -d "." -f9 | uniq | while read UPfile; do mkdir ${file}/${UPfile}; mv ${file}/*.${UPfile}.* ${file}/${UPfile}; done; done

#### grab all unique transcript IDs upregulated in each condition vs all other conditions
```
cd unique
for i in "coni" "glen" "symp" "para" "spor"; do ls ../ | grep $i | while read file; do cat ../${file}/${i}-UP/matrix.counts.matrix.${file}.DESeq2.DE_results.P0.001_C2.${i}-UP.subset | awk -F "\t" '$1!="sampleA" {print $1}'; done | sort | uniq -c | awk '$1==4 {print $2}' > ${i}-UP.id.txt; done
```

#### get annotations for all previously identified transcripts
    for i in "coni" "glen" "symp" "para" "spor"; do cat ${i}-UP.id.txt | while read line; do echo $line; cat ~/Documents/Data/Annotations/Emuscae_annotation_report.xls | grep -w $line | awk -F "\t" '{if($7!=".") {split($7,blastp,"^"); print blastp[6]} else {print $7}}; {if($10!=".") {split($10,pfam,"^"); print pfam[3]} else {print $10}}; {print $11}; {if($12!=".") {split($12,tmhmm,"^"); print tmhmm[1]" "tmhmm[2]} else {print $12}}' | sed 's/RecName\: Full\=//g'; ls ../ | grep $i | while read file; do cat ../${file}/${i}-UP/matrix.counts.matrix.${file}.DESeq2.DE_results.P0.001_C2.${i}-UP.subset | grep -w $line | cut -f7,11; done; done | paste - - - - - - - - - > ${i}-UP.xls ; done

#### divide DE_results files to give only those upregulated in each individual pairwise condition
```
cd ../
ls | grep _vs_ | grep -v pdf | while read file; do first=`echo $file | cut -d "_" -f1`; second=`echo $file | cut -d "_" -f3`; awk -F "\t" '$7>=0 {print $0}' ${file}/matrix.counts.matrix.${file}.DESeq2.DE_results > ${file}/${first}-UP/${first}-UP.DE_results ; awk -F "\t" '$7<0 {print $0}' ${file}/matrix.counts.matrix.${file}.DESeq2.DE_results > ${file}/${second}-UP/${second}-UP.DE_results ; done
```

#### create annotations for each pairwise comparison against a given condition (ID, blastp, pfam, sigp, tmhmm, logFC, FDR, logFC, FDR, logFC, FDR, logFC, FDR)
```
~/scripts/DE_concatenate_pairs.py ~/Documents/Data/Annotations/Emuscae_annotation_report.xls coni_vs_glen/coni-UP/coni-UP.DE_results coni_vs_symp/coni-UP/coni-UP.DE_results coni_vs_para/coni-UP/coni-UP.DE_results coni_vs_spor/coni-UP/coni-UP.DE_results vs_all/coni.temp.xls

~/scripts/DE_concatenate_pairs.py ~/Documents/Data/Annotations/Emuscae_annotation_report.xls coni_vs_glen/glen-UP/glen-UP.DE_results glen_vs_symp/glen-UP/glen-UP.DE_results glen_vs_para/glen-UP/glen-UP.DE_results glen_vs_spor/glen-UP/glen-UP.DE_results vs_all/glen.temp.xls

~/scripts/DE_concatenate_pairs.py ~/Documents/Data/Annotations/Emuscae_annotation_report.xls coni_vs_symp/symp-UP/symp-UP.DE_results glen_vs_symp/symp-UP/symp-UP.DE_results para_vs_symp/symp-UP/symp-UP.DE_results spor_vs_symp/symp-UP/symp-UP.DE_results vs_all/symp.temp.xls

~/scripts/DE_concatenate_pairs.py ~/Documents/Data/Annotations/Emuscae_annotation_report.xls coni_vs_para/para-UP/para-UP.DE_results glen_vs_para/para-UP/para-UP.DE_results para_vs_symp/para-UP/para-UP.DE_results para_vs_spor/para-UP/para-UP.DE_results vs_all/para.temp.xls

~/scripts/DE_concatenate_pairs.py ~/Documents/Data/Annotations/Emuscae_annotation_report.xls coni_vs_spor/spor-UP/spor-UP.DE_results glen_vs_spor/spor-UP/spor-UP.DE_results spor_vs_symp/spor-UP/spor-UP.DE_results para_vs_spor/spor-UP/spor-UP.DE_results vs_all/spor.temp.xls
```

#### make all DE annotation files more readable
```
cd vs_all
ls | grep temp | while read file; do newfile=`echo $file | sed 's/\.temp//g'`; awk -F "\t" '{split($2,a,"^"); split($3,b,"^"); split($5,c,"^"); {print $1"\t"a[1]" "a[6]"\t"b[1]" "b[3]"\t"$4"\t"c[1]" "c[2]"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}}' ${file} | sed 's/RecName\: Full\=//g' > ${newfile} ; done
```


------------------------------------
### Gene Set Enrichment Analysis ### 

#### Starting directory
    ~/Documents/Data/DiffExpr/GSEA

#### create rank files from DESeq2 results (sign of logFC * inverse of FDR)
```
mkdir Rank_Files
cd ../P0.001
ls | grep _vs_ | grep -v pdf | while read file; do id=`echo $file | sed 's/\_vs\_/\_versus\_/g'`; cat ${file}/matrix.counts.matrix.${file}.DESeq2.DE_results | tail -n +2 | awk -F "\t" '{if($7>=0) {printf "%s\t%e\n", $1,(1/$11)} else {printf "%s\t%e\n", $1,(-1/$11)}}' | sort -grk2 > ../GSEA/Rank_Files/${id}.rnk; done
```

#### create Rank File opposites
```
cd ../GSEA/Rank_Files
ls | while read file; do first=`echo $file | cut -d "_" -f1`; second=`echo $file | cut -d "." -f1 | cut -d "_" -f3`; cat $file | sort -gk2 > ${second}_versus_${first}; done
cd ..
ls | grep -v rnk | while read file; do ~/scripts/GSEA_flip_rank_files.py ${file}; rm ${file}; done
```

#### create GMT file for GO terms prior to GSEA
```
cd ..
cat ~/Documents/Data/Annotations/go_annotations.txt | cut -f2 | tr "," "\n" | sort | uniq | while read line; do echo -ne "${line}\t"; cat ~/Documents/Data/Annotations/go_annotations.txt | grep $line | cut -f1 | tr "\n" "\t"; echo; done > GO.temp.gmt
~/scripts/GSEA_parse_GO_terms.py go.obo GO.temp.gmt GO.term.gmt
```

#### create GMT file for KEGG ORTHOLOGY (KO) terms prior to GSEA (query.ko downloaded from KAAS output)
```
wget "http://www.genome.jp/kegg-bin/download_htext?htext=ko00000.keg&format=htext&filedir=" -O ko00000.keg
cat ko00000.keg | grep -Eo "K[0-9]{5}" | sort | uniq | while read line; do echo -ne "${line}\t"; cat ko00000.keg | grep -m1 $line | cut -d ";" -f2; done > ko00000.txt
cat query.ko | grep -E "K[0-9]{5}" | cut -f2 | sort | uniq | while read line; do echo -ne "${line}\t"; cat query.ko | grep $line | cut -d ":" -f1 | tr "\n" "\t"; echo; done > KO.temp.gmt
cat KO.temp.gmt | cut -f1 | while read line; do cat ko00000.txt | grep $line | tr "\n" "\t"; cat KO.temp.gmt | grep $line | cut -f2- ; done | awk -F "\t" '{printf $1} {for(i=2;i<NF;i++) printf "\t"$i} {printf "\n"}' > KO.term.gmt
```

#### create GMT file for KEGG PATHWAYS prior to GSEA (using previously downloaded ko00000.keg file)
```
~/scripts/GSEA_parse_KEGG_terms.py ko00000.keg KO.term.gmt KEGG.temp.gmt
awk -F "\t" '{printf $1} {for(i=2;i<NF;i++) printf "\t"$i} {printf "\n"}' KEGG.temp.gmt | awk -F "\t" '{if(NF>2) {print $0}}' > KEGG.term.gmt
```

#### create results folders for each pairwise comparison
```
mkdir Results
ls Rank_Files/ | while read file; do id=`echo $file | sed 's/\.rnk//g'`; mkdir Results/${id}; done
cd Results
```

#### run pre-ranked GSEA for all pairwise comparisons between GO terms (classic, 1000 permutations, exclude groups >100 & <8)
    ls | grep "_versus_" | while read file; do java -cp ~/bin/gsea-3.0_beta_3.jar xtools.gsea.GseaPreranked -gmx ~/Documents/Data/DiffExpr/GSEA/GO.term.gmt -norm meandiv -nperm 1000 -rnk ~/Documents/Data/DiffExpr/GSEA/Rank_Files/${file}.rnk -scoring_scheme classic -rpt_label GO.${file} -create_svgs false -make_sets true -plot_top_x 200 -rnd_seed timestamp -set_max 100 -set_min 8 -zip_report false -out ~/Documents/Data/DiffExpr/GSEA/Results/${file} -gui false; done

#### run pre-ranked GSEA for all pairwise comparisons between KO terms (classic, 1000 permutations, exclude groups >1000 & <2)
    ls | grep "_versus_" | while read file; do java -cp ~/bin/gsea-3.0_beta_3.jar xtools.gsea.GseaPreranked -gmx ~/Documents/Data/DiffExpr/GSEA/KO.term.gmt -norm meandiv -nperm 1000 -rnk ~/Documents/Data/DiffExpr/GSEA/Rank_Files/${file}.rnk -scoring_scheme classic -rpt_label KO.${file} -create_svgs false -make_sets true -plot_top_x 200 -rnd_seed timestamp -set_max 1000 -set_min 2 -zip_report false -out ~/Documents/Data/DiffExpr/GSEA/Results/${file} -gui false; done

#### run pre-ranked GSEA for all pairwise comparisons between KEGG pathways (classic, 1000 permutations, exclude groups >1000 & <2)
    ls | grep "_versus_" | while read file; do java -cp ~/bin/gsea-3.0_beta_3.jar xtools.gsea.GseaPreranked -gmx ~/Documents/Data/DiffExpr/GSEA/KEGG.term.gmt -norm meandiv -nperm 1000 -rnk ~/Documents/Data/DiffExpr/GSEA/Rank_Files/${file}.rnk -scoring_scheme classic -rpt_label KEGG.${file} -create_svgs false -make_sets true -plot_top_x 200 -rnd_seed timestamp -set_max 1000 -set_min 2 -zip_report false -out ~/Documents/Data/DiffExpr/GSEA/Results/${file} -gui false; done

#### create annotations for each pairwise comparison against a given condition (repeat for GO, KEGG, KO files)
```
mkdir vs_all vs_all/GO vs_all/KEGG vs_all/KO

ls vs_all | while read file; do ~/scripts/DE_concatenate_pairs.py ../${file}.term.gmt coni_versus_glen/${file}.coni_versus_glen.GseaPreranked.*/gsea_report_for_na_pos_*.xls coni_versus_symp/${file}.coni_versus_symp.GseaPreranked.*/gsea_report_for_na_pos_*.xls coni_versus_para/${file}.coni_versus_para.GseaPreranked.*/gsea_report_for_na_pos_*.xls coni_versus_spor/${file}.coni_versus_spor.GseaPreranked.*/gsea_report_for_na_pos_*.xls vs_all/${file}/coni.xls ; done

ls vs_all | while read file; do ~/scripts/DE_concatenate_pairs.py ../${file}.term.gmt glen_versus_coni/${file}.glen_versus_coni.GseaPreranked.*/gsea_report_for_na_pos_*.xls glen_versus_symp/${file}.glen_versus_symp.GseaPreranked.*/gsea_report_for_na_pos_*.xls glen_versus_para/${file}.glen_versus_para.GseaPreranked.*/gsea_report_for_na_pos_*.xls glen_versus_spor/${file}.glen_versus_spor.GseaPreranked.*/gsea_report_for_na_pos_*.xls vs_all/${file}/glen.xls ; done

ls vs_all | while read file; do ~/scripts/DE_concatenate_pairs.py ../${file}.term.gmt symp_versus_coni/${file}.symp_versus_coni.GseaPreranked.*/gsea_report_for_na_pos_*.xls symp_versus_glen/${file}.symp_versus_glen.GseaPreranked.*/gsea_report_for_na_pos_*.xls symp_versus_para/${file}.symp_versus_para.GseaPreranked.*/gsea_report_for_na_pos_*.xls symp_versus_spor/${file}.symp_versus_spor.GseaPreranked.*/gsea_report_for_na_pos_*.xls vs_all/${file}/symp.xls ; done

ls vs_all | while read file; do ~/scripts/DE_concatenate_pairs.py ../${file}.term.gmt para_versus_coni/${file}.para_versus_coni.GseaPreranked.*/gsea_report_for_na_pos_*.xls para_versus_glen/${file}.para_versus_glen.GseaPreranked.*/gsea_report_for_na_pos_*.xls para_versus_symp/${file}.para_versus_symp.GseaPreranked.*/gsea_report_for_na_pos_*.xls para_versus_spor/${file}.para_versus_spor.GseaPreranked.*/gsea_report_for_na_pos_*.xls vs_all/${file}/para.xls ; done

ls vs_all | while read file; do ~/scripts/DE_concatenate_pairs.py ../${file}.term.gmt spor_versus_coni/${file}.spor_versus_coni.GseaPreranked.*/gsea_report_for_na_pos_*.xls spor_versus_glen/${file}.spor_versus_glen.GseaPreranked.*/gsea_report_for_na_pos_*.xls spor_versus_symp/${file}.spor_versus_symp.GseaPreranked.*/gsea_report_for_na_pos_*.xls spor_versus_para/${file}.spor_versus_para.GseaPreranked.*/gsea_report_for_na_pos_*.xls vs_all/${file}/spor.xls ; done
```

-----------------------------------------------------
### Protein Alignments ### (example with PF00089) ###

#### Starting directory
    ~/Documents/Data/Alignments

#### obtain start and end positions for pfam domains (PF00089) in relevant sequences
```
mkdir PF00089
cd PF00089
cat ~/Documents/Data/Annotations/Emuscae_annotation_report.xls | grep "PF00089" | cut -f1,10 | awk -F "\t" '{split($2,pfam,"\`"); for(i=1; i<=length(pfam); i++) {print $1"\t"pfam[i]}}' | grep "PF00089" | awk -F "\t" '{split($2,pfam,"^"); {print $1"\t"pfam[4]}}' | tr "-" "\t" > PF00089.lengths.txt
```

#### use custom script to slice pfam domain from selected sequences
    ~/scripts/alignments_filter_domains_from_lengths.py PF00089.lengths.txt ~/Documents/Data/Assembly/Trinity_isoH_cdhit10.fasta.transdecoder.pep PF00089.fasta

#### run cd-hit to remove duplicate sequences
    cd-hit-est -i PF00089.fasta -o PF00089_cdhit10.fasta -c 1.0 -n 10

#### carry out protein alignment using hmmalign and downloaded hmm profile
    hmmalign PF00089.hmm PF00089_cdhit10.fasta > PF00089.hmmalign.stk

