BConsScore <- function(seq, eval="1e-4", remote=TRUE, db="refseq_protein"){
	if(remote){
		blast.xml <- paste(system(paste("blastp -remote -db ", db, " -outfmt 5 -evalue ", eval, " -query ", seq), intern=TRUE), collapse="")
	}
	tmp <- tempfile()
	Blast_fasta_Name <- paste(tmp, "blastFASTA", sep = ".")
	ClustalOutName <- paste(tmp, "clustal", sep = ".")
	
	#Writes the information of interest from the blast-query to a fasta-format file
	blast.tree <- xmlInternalTreeParse(blast.xml, asText=TRUE)
	hit_IDs <- xpathSApply(blast.tree, "//Hit/Hit_id", xmlValue)
	hit_seqs <- xpathSApply(blast.tree, "//Hit/Hit_hsps/Hsp/Hsp_hseq", xmlValue)
	write(paste(">", hit_IDs, "\n", hit_seqs, "\n", sep=""), file=Blast_fasta_Name)
	
	#Makes a multiple alignment from the fasta file using clustalOmega
	system(paste("clustalo -i ", Blast_fasta_Name," -o ", ClustalOutName))
	
	#Computes the conservation scores using mstatx
	system(paste("/Users/simon/workspace/SLU/MutImpact_project/MstatX/mstatx -m ", ClustalOutName))
	Scores <- read.delim(paste(tmp, "stat", sep="."), header=FALSE)
}