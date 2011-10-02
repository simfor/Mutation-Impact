BlastSearch <- function(seq, eval="1e-2", remote=TRUE, db="refseq_protein", multal=FALSE){
	#paste(system(paste("blastp -remote -db cdd -outfmt 5 -evalue 1e-2 -query ", seq), intern=TRUE), collapse="")
	if(remote){
		blast <- paste(system(paste("blastp -remote -db ", db, " -outfmt 5 -evalue ", eval, " -query ", seq), intern=TRUE), collapse="")
	}
	
	#Converts the blast output from XML to fasta format
	if(multal){
		tmp <- tempfile()
		blast_XML_Name <- paste(tmp, "blastXML", sep = ".")
		perlName <- paste(tmp, "pl", sep = ".")
		Blast_fasta_Name <- paste(tmp, "blastFASTA", sep = ".")
		ClustalOutName <- paste(tmp, "clustal", sep = ".")
		
		#Creates temp-file with the XML output from blastp
		write(blast, file = blast_XML_Name)
		
		#Creates temp-file with the perl-script for converting XML => fasta
		write(paste("my $blast='", blast_XML_Name, "';", sep = ""), file = perlName, append = TRUE)
		write(paste("open(OUT,'>", Blast_fasta_Name, "');", sep = ""), file = perlName, append = TRUE)
		#write(readLines(file.path(.path.package("MutationImpact"), "scripts", "proteus2")), file = perlName, append = TRUE)
		write(readLines(file.path("/Users/simon/workspace/SLU/MutImpact_project/Mutation-Impact/blastp_XML_to_fasta.pl")), file = perlName, append = TRUE)
		
		#Calls the perl-script which converts the blast XML ouput to fasta format
		system(paste("perl ", perlName))
		#Calls clustal which makes a multiple alignment from the fasta file
		system(paste("clustalo -i ", Blast_fasta_Name," -o ", ClustalOutName))
	}
	#list()
	readLines(ClustalOutName)
}



ConsScore <- function(seq, eval="1e-4", remote=TRUE, db="refseq_protein"){
	if(remote){
		blast <- paste(system(paste("blastp -remote -db ", db, " -outfmt 5 -evalue ", eval, " -query ", seq), intern=TRUE), collapse="")
	}
	
	tmp <- tempfile()
	blast_XML_Name <- paste(tmp, "blastXML", sep = ".")
	perlName <- paste(tmp, "pl", sep = ".")
	Blast_fasta_Name <- paste(tmp, "blastFASTA", sep = ".")
	ClustalOutName <- paste(tmp, "clustal", sep = ".")
	MstatxOutName <- paste(tmp, "MstatX", sep = ".")
	
	#Creates temp-file with the XML output from blastp
	write(blast, file = blast_XML_Name)
	
	#Creates temp-file with the perl-script for converting XML => fasta
	write(paste("my $blast='", blast_XML_Name, "';", sep = ""), file = perlName, append = TRUE)
	write(paste("open(OUT,'>", Blast_fasta_Name, "');", sep = ""), file = perlName, append = TRUE)
	#write(readLines(file.path(.path.package("MutationImpact"), "scripts", "proteus2")), file = perlName, append = TRUE)
	write(readLines(file.path("/Users/simon/workspace/SLU/MutImpact_project/Mutation-Impact/blastp_XML_to_fasta.pl")), file = perlName, append = TRUE)
	
	#Calls the perl-script which converts the blast XML ouput to fasta format
	system(paste("perl ", perlName))
	#Calls clustal which makes a multiple alignment from the fasta file
	system(paste("clustalo -i ", Blast_fasta_Name," -o ", ClustalOutName))
	
	system(paste("/Users/simon/workspace/SLU/MutImpact_project/MstatX/mstatx -m ", ClustalOutName))
	Scores <- read.delim(paste(tmp, "stat", sep="."))
}

