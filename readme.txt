test2.R is the code for identifing lncRNA target genes.

all.lncrna.txt and allgene.txt are the lncRNAs and genes' names.

gene-net1.txt is the gene-gene interaction data which is downloaded from humannetv2.0

The features of lncRNAs and gene have already been extracted as lncrna.feature.txt and gene.feature.txt.

lncrna-gene.net.txt is the potential correlation between lncRNAs and genes. This is obtained based on their interactions with miRNAs.

lncrna_target.txt is the high-throughput known data from LncRNA2Target v2.0

test2.R includes four parts.

1.normalize features of gene and lncRNAs.
2.Encode gene network by GCN
3. generate negative sample set.
4. Classify by CNN

   
