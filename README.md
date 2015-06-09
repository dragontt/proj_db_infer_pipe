# proj_db_infer_pipe

This project is an integration of regulatory network mapping algorithms and motif inference algorithms. The networks are inferred from expression data (microarray / RNA-seq) using NetProphet and BART. The DNA binding motifs are inferred from motif database using integrated information of DBD similarity and network correlation, or are inferred using de novo motif inference algorithm (FIRE). 


### Network Inference Methods

* NetProphet 

Computes LASSO regression of gene co-expression with global shrinkage, and combine the network with gene differnetial expression (DE) network. Check NetProphet package README for details. The DE network is built using *.pipe: LIMMA for microarray data, TopHat and CuffLink for RNA-seq data. Check *.pipe (yeast, crypto, or fly) for details.

* BART

Comptues ensemble of Bayesian regression trees, which are weak learners contrained by prior. Check BART cran.r-project manuscript for details.


### Resources Preparation

* Gene Expression Data

- DMel RNA-seq data generated from gene perturbation (Baranski lab).
- DMel microarray data from 9 external resources, one of which has gene perturbation. 

* DNA Binding Motif Database

- CIS-BP: A database that collects directly evidenced motifs (in forms of PFM) and their associated DBD sequence. Check Weirauch MT paper for details.

* Fly Sequence Data

- FlyBase: A comprehensive fly database. Useful files include:
	- Conversion table between gene name and protein name: FBgn <=> FBtr <=> FBpp IDs.
	- DMel protein fasta sequence. Use bedtools getfasta to identify DBD sequence.
	- Gene annotation and ontology files.
- RSAT Metazoa: DMel promoter regions with specified range. 

* Fly Physical Network (for Network Evaluation)

- FlyNet: Binary ChIP and PWM binding networks.
- ChIP network curated from multiple sources (majorly from modENCODE).
   Motif network built by aligning and scoring the known DMel motifs on target promoters.


### Useful Bioinformatics Tools

* MEME Suite

- FIMO: aligns and scores motif PFM on target promoters.
- TOMTOM: compares motifs.
- MEME motif database

* FIRE (de novo Motif Inference Tool)

Infers motifs by computing mutual information of motif profiles and gene expression profiles. The expression profiles are quantized mapped network edge scores. Check FIRE paper from Tavazoie lab for details.

* Clustal Omega

Multi sequnece alginment tool used to align DBD sequences and compute percent identity scores.
