# 4F
4F is a simple R script to perform sliding window-based analyses on four-fold degenerate codons of a mitochondrial genome.
Basically, two statistics are computed: nucleotide composition and A-T content.
The nucleotide composition at fourfold degenerate sites is related with single-strand state duration during mitochondrial DNA
(mtDNA) replication. The more the heavy (H) strand remains unpaired, the more the spontaneous hydrolytic deamination of C to
U and A to hX (hypoxanthine) takes place. Such an increase of T and hX in the H strand leads to a corresponding increase in
the percentages of A and C in the complementary lagging (L) strand where the H strand remains for longer time in the
single-stranded condition, i.e. near to the Origin of Replication (OR). Moreover, single-stranded-guanine may spontaneously
oxidize to 8-hydroxyguanine, which basepairs with adenine: thus, in this case, G decreases and T increases on the H strand.
In a nutshell, T will only tend to accumulate near to the OR of the H strand, while the opposite is true for A and C; finally,
G may behave in either way. This asymmetrical composition can leave a neutral signature in
fourfold degenerate sites, being them under no or weak selection.
The A-T skew at four-fold degenerate sites is known to be correlated with the position of the ORs as well: extreme (i.e.,
closer to +1 and -1) values are associated with Protein Coding Genes (PCGs) located near to the OR of the H strand, while
balanced (i.e., closer to $0$) values are associated with PCGs located near to the OR of the L strand.
