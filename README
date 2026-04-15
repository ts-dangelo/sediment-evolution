This code documents the analysis in the publication. They are named in order of operation.

00-nf-mag-process-reads.sh

	Quality control raw reads downloaded from the Marine Metagenome DB downloaded from ERR/SRA
	using the SOP of nf-core/mag
	
01-make-bowtie2-bams.nf

	Map the QC reads to our highest quality SAG in the Aerophobota AE-B3B sp005223085: AM-917-D13
	
	Data used to make Figure 1B and assess which samples to analyze downstream
	
02-make_bowtie2-95_allSAGs_for_spades.nf

	In order to capture as many reads corresponding to this species as possible, we remapped
	the same reads to a concatenated fasta of all SAGs from the species cluster
	
03-spades.nf

	Creating SPades assemblies from the mapped reads

04-run_prodigal.sh

	Create translated genome files from the MAGs produced above and the SAGs
	
05-create_mmseqs_clusters.txt

	Creating gene clusters from the prodigal files from the SAGs and MAGs using mmseqs

06-make_hyphy_files.nf

	Takes the gene cluster amino acid and nucelotide fasta files, creates codon aware nucleotide alignments,
	gene trees with iqtree and runs evolutionary tests using several tools in the HyPhy package


MMDB_samples.txt: Marine Metagenome Database metadatafile describing the samples downloaded (displayed in Figure 1#)

bt2_95id_10x-AM-917-D13-mean_coverage_samples.txt: 
	Samples from MMDB_samplex that had greater than 10X coverage to SAG AM-917-D13. These samples were focused on for read-mappiung and binning analysese.
	calculations perfomed by samtools_depth_calc.sh
