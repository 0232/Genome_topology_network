# Genome_topology_network
1. The draft genomes now can be analyzed.
Initiate the GTN with "gtn_CompleteOnly.pl" if your data consist of complete genomes only, or utilize "gtn_WithDraft.pl". The GTN firstly locates the synteny block in genomes if draft data exist.

2. The Markov Cluster Algorithm is introduced to assign the gene families.
The Markov Cluster (MCL) Algorithm is introduced in this version. This algorithm is also used in many clustering tool, such as orthoMCL,

3. Unique gene connection information will be yielded as a result file.
GTN is an approach which studies closely-related bacterial genomes by analyzing the adjacent gene family pairs in genomes. The evolution distance calculation is essentially based on the different relationship of gene connection. In other words, the relationship of gene connection influences the separation of genomes in phylogenetic tree. The unique gene connection reflects the evolution events, such as gene duplication, gene insertion.

For more detail, please refer to 'GTN_document.pdf' in this website.
