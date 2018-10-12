# Structator: fast index-based search for RNA sequence-structure patterns

Fernando Meyer, Stefan Kurtz, Rolf Backofen, Sebastian Will and Michael Beckstette. BMC Bioinformatics. 2011;12:214. doi: 10.1186/1471-2105-12-214

## Application examples

* Searching a subset of Rfam release 10.0 for RNA family CTV_rep_sig (Rfam Acc. RF00193) by building high-scoring global chains of matches

Files required\
Data: RFAM10_8MB.fa\
Secondary structure descriptor (SSD): RF00193.pat\
Alphabet: rna.alphab\
Watson-Crick and wobble complementarity rules: dna_rna.comp

Sample command for the index construction\
./afconstruct RFAM10_8MB.fa -alph rna.alphab -a -s indexname1

Sample command for the search\
./afsearch indexname1 -pat RF00193.pat -comp dna_rna.comp -a -global -minlen 5


* Searching human chromosome 20 for RNA gene HAR1F (Rfam Acc. RF00635) by building high-scoring local chains of matches

Files required\
Data: ftp://ftp.ensembl.org/pub/release-60/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.60.dna.chromosome.20.fa.gz\
The file can be unpacked with gunzip (http://www.gzip.org/)\
SSD: RF00635.pat\
Alphabet: dna.alphab\
Watson-Crick and wobble complementarity rules: dna_rna.comp

Sample command for the index construction\
./afconstruct Homo_sapiens.GRCh37.60.dna.chromosome.20.fa -alph rna.alphab -a -s indexname2

Sample command for the search\
./afsearch indexname2 -pat RF00635.pat -comp dna_rna.comp -a -local -wf 10 -show