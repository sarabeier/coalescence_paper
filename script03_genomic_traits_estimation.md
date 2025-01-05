# Command lines and script for the computation of genomic traits

### by Sara Beier & Angel Rain

The current file contains shell command lines to extraplotate genomic traits from metabarcoding 16s rRNA reads via the PICRUSt2 software (Douglas et al. 2020). Trait values for the PICRUSt2 default reference strains are available from an earlier study (Table S1, Beier et al. 2022)

## Tranforming trait values
The PICRUSt2 software accepts only trait values between 0-1000 with maximal two digits for trait prediction and raw trait values were ranged to meet this criteria in R:

```R
library(tidyverse)

data.range<-function(n)
{
  # transform to range from 0-1000, round by two digits
  round((n-min(n,na.rm=T))/(max(n,na.rm=T)-min(n,na.rm=T))*1000,2)
}

gtraits <- tibble(read.table("Data/TableS1.tsv", header = T, sep = "\t", fill = T)) %>%
  select(PICRUSt.ID,Genome.size,X.TF, Generation.time..gRodon.)

ptraits <- data.frame(assembly=gtraits$PICRUSt.ID,
                      genome.size=data.range(gtraits$Genome.size),
                      TF_perc=data.range(gtraits$X.TF),
                      d.gRodon=data.range(log(gtraits$Generation.time..gRodon.)))
                      

var <-c('genomesize','TF_perc','cub.dRg')
var.min <- c(min(gtraits[2],na.rm=T), min(gtraits[3],na.rm=T),min(gtraits[,4],na.rm=T))
var.max <- c(max(gtraits[2],na.rm=T),max(gtraits[3],na.rm=T), max(gtraits[,4],na.rm=T))
var.trans <- c('NA','NA','log(x)')

pic.traits <- data.frame(var, var.min, var.max,var.trans)
pic.traits
# Export transformation factors
write.table (pic.traits, 'Data/PICRUSt_predictions/trait_info_coal.tsv', sep='\t', row.names=FALSE)

# Export only traits associated to resilience and resistance traits as single files for picrust
## Resistance-related traits
write.table(na.omit(ptraits[c(1,3)]), 'data/picrust2.output/trait.data/p.TF_perc.txt', row.names=FALSE, sep ='\t',quote=FALSE)
write.table(na.omit(ptraits[c(1,2)]), 'data/picrust2.output/trait.data/p.genome.size.txt', row.names=FALSE, sep ='\t',quote=FALSE)
## Resilience-related traits
write.table(na.omit(ptraits[c(1,4)]), 'data/picrust2.output/trait.data/p.d.gRodon.txt', row.names=FALSE, sep ='\t',quote=FALSE)
```

Citation:  
* Beier,S., Werner,J., Bouvier,T., Mouquet,N., Violle, C. (2022). Trait-trait relationships and tradeoffs vary with genome size in prokaryotes. doi: https://doi.org/10.3389/fmicb.2022.985216
* Douglas, G. M., Maffei, V. J., Zaneveld, J. R., Yurgel, S. N., Brown, J. R., Taylor, C. M., et al. (2020). PICRUSt2 for prediction of metagenome functions. Nature Biotechnology 38, 685–688. https://doi.org/10.1038/s41587-020-0548-6.


## Trait-based estimations of genome size, %TF, and generation time using PICRUSt2
Trait values for generation time, genome size and %TF were predicted using the the default reference database from picrust using the ranged data created above. The fasta file with 16s rRNA gene sequence data (dada-chem.seqs.nochim.fasta) as well as the corresponding count table (dada-chem.rcounts.ASV.2.tab) were created as detailed in script02_community.dada.R.


```bash
# Edit fasta.files sequences (Convert to multiline fastafiles)
sed 's/ /&\n/g' dada-cryo.seqs.nochim.fasta |sed 's/ //g'| sed 's/(Marine_group_B)//g'|sed 's/SV/>SV/g' >picr_dada-cryo.seqs.nochim.fasta #remove spaces after seqIDs, add '>' before sedID, insert line brack after seqIDs and remove brackets in seqIDs (produce error!)

# PICRUST2 full pipeline
picrust2_pipeline.py -s picr_dada-cryo.seqs.nochim.fasta -i dada-cryo.rcounts.ASV.2.tab --no_pathways --min_samples 0 -o cryo.picrust_out_full -p 25

# Hidden state prediction (hsp) using custom trait tables
# Genome size
hsp.py --observed_trait_table trait.data/p.genome.size.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.genomesize_predicted -n -p 10
#%TF
hsp.py --observed_trait_table trait.data/p.TF_perc.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.TF_perc_predicted -n -p 10
# generation time (gRodon)
hsp.py --observed_trait_table trait.data/p.d.gRodon.txt -t picrust2_out_full/out.tre -o picrust.traits/pic.dgR_predicted -n -p 10 
```

# 16S copy number (RRN) prediction from custom tree using PICRUSt2
Analyses on the default PICRUSt2 RRN prediction indicated that the trait values assigned to the default reference database might be biased (Figure S4, S5 in Beier et al. 2022). We therefore predicted RRN from the manually curated rrnDB database (Stoddard et al. 2015) using rrnDB trait values (https://rrndb.umms.med.umich.edu/static/download/, rrnDB-5.7.tsv). A custom phylogenetic tree was created based on the rrnDB phylogeny as detailed in Beier et al. 2022 (https://github.com/sarabeier/genomic.traits, scripts_phylogenetic.signal.md). For RRN prediction from a custom reference tree PICRUST2 v2.4.2 was used since reference tree application was implemented after version 2.4.0.

Format rrnDB database in R
```R
## RRN trait table from from rrnDB: rrnDB-5.7.tsv
rrnDBtraits <- read.csv ("rrnDB-5.7.tsv", sep='\t', header=T, fill=T)[-c(1:215),c(1,12)] #rrnDB trait table, exclude first rows, no sequence data available
names(rrnDBtraits)=c("assembly","16S_rRNA_Count")
write.table(rrnDBtraits,file="picrust2.output/trait.data/16S_rrnDB.txt",row.names=FALSE, sep ='\t',quote=FALSE)
```
Run prediction in the bash shell
```bash
#Setting up the custom reference tree
cd picrust2.output/rrnDB_customREF

#Copy rrnDB alignment and tree into rrnDB_customREF directory
cp /rrnDB/muscle2.rrnDB.fasta ../picrust2.output/rrnDB_customREF
cp /rrnDB/rrndb_fasttree.format.tree ../picrust2.output/rrnDB_customREF

# Check aligment sequences
conda activate raxmlng-1.0.3
raxml-ng --check --msa muscle2.rrnDB.fasta --model GTR+G --prefix T1

#Evaluating the Model Parameters (https://github.com/Pbdas/epa-ng#setting-the-model-parameters)
nohup raxml-ng --evaluate --msa muscle2.rrnDB.fasta --tree rrndb_fasttree.format.tree --prefix info --model GTR+G+F --threads 10 &

# Build a hmm profile for the aligned seq
conda activate qiime2-2021.2
nohup hmmbuild muscle2.rrnDB.hmm muscle2.rrnDB.fasta &

#Re-name reference files for prediction (needs to be re-named to rrnDB_customREF*)
mv muscle2.rrnDB.fasta rrnDB_customREF.fasta
mv muscle2.rrnDB.hmm rrnDB_customREF.hmm
mv rrndb_fasttree.format.tree rrnDB_customREF.tre
mv info.raxml.bestModel rrnDB_customREF.model

#Delete uneccesary intermediate files and gzip rrnDB_customREF.fasta
rm info.raxml.*
rm T1.raxml.log
rm nohup.out
gzip rrnDB_customREF.fasta


# Run picrust pipeline with custom tree
## 16S copy number trait from Beier et al., 2022 (16S_rrnDB.txt file, see above)

picrust2_pipeline.py -s picr_dada-cryo.seqs.nochim.fasta -i ../dada.out.coal24/dada-cryo.rcounts.ASV.2.tab --min_samples 0 -o cryo.picrust_out_custom -p 25 --ref_dir /bio/Databases/GenomicTraitsPrediction_PICRUST2/rrnDB_customREF/ --marker_gene_table /bio/Databases/GenomicTraitsPrediction_PICRUST2/rrnDB_customREF/16S_rrnDB.txt --no_pathways --skip_minpath --edge_exponent 0

# Rename and unzip output file
mv cryo.picrust_out_custom/marker_predicted_and_nsti.tsv.gz ./pic.16S_predicted_rrnDB.gz
gunzip pic.16S_predicted_rrnDB.gz
```
Citation:  
* Beier,S., Werner,J., Bouvier,T., Mouquet,N., Violle, C. (2022). Trait-trait relationships and tradeoffs vary with genome size in prokaryotes. Front. Microbiol. doi: https://doi.org/10.3389/fmicb.2022.985216
* Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. K., and Schmidt, T. M. (2015). rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucleic Acids Res. 43, D593–D598. doi:10.1093/nar/gku1201.
