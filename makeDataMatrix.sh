##Date:18/08/2024
##Author: Sarah Grant
##Description: Take the data from the STAR alignment and genecount mode and builds gene x sample matrix
### The "ReadsPerGene.out.tab" files contain four columns the col1:geneID, col2:unstranded counts, col3:1st read strand counts, col4:2nd read strand counts
### for the dUTP method used for the experiment the 4th column is of interest.

## set path & make output directory
cd /HANGER/sarahGrant/RNA-seq/align
mkdir /HANGER/sarahGrant/RNA-seq/deseq2

# retrieve the 4th column of each "ReadsPerGene.out.tab" file + the first column that contains the gene IDs
paste  V350134869*ReadsPerGene.out.tab | grep -v "_" | awk '{printf "%s\t", $1}{for (i=4;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp

# add header: "gene_name" + the name of each of the counts file
sed -e "1igene_name\t$(ls V350134869*ReadsPerGene.out.tab | tr '\n' '\t' | sed 's/.ReadsPerGene.out.tab//g')" tmp | cut -f1-40 > ../deseq2/V350134869_raw_count_matrix.txt

# remove temporary file
rm tmp
