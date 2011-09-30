library(BioSeqClass)
library(Biostrings)
library(tkrplot)
library(XML)

MutationImpact <- function(ff1, ff2, secstr=TRUE, dom=TRUE){
	protA <- readFASTA(ff1, strip.descs=TRUE)
	protB <- readFASTA(ff2, strip.descs=TRUE)
	
	#Aligns the two sequences
	print("Aligning sequences", quote=FALSE)
	aligned <- pairwiseAlignment(pattern = protA[[1]]$seq, subject = protB[[1]]$seq, type="global", substitutionMatrix = "BLOSUM50", gapOpening = -5, gapExtension = -1)
	
	#Retrieves the secondary structure for protA
	if(secstr){
		print("Retrieving secondary structure", quote=FALSE)
		protA_SecondaryStructure <- SecStructure(protA[[1]]$seq, strsplit(toString(pattern(aligned)), "")[[1]] )
	}
	else
		protA_SecondaryStructure <- list("", PSIPRED=list(ConfidenceScore=0))
	
	#Calculates the distance between the two sequences
	D <- Distance(toString(pattern(aligned)), toString(subject(aligned)))
	
	#Calculates the weight matrix and returns the different distance scores
	dist <- Weight(protA_SecondaryStructure[[1]], protA_SecondaryStructure[[2]]$PSIPRED$ConfidenceScore, D)
	
	#Looks for conserved domains in protA
	if(dom){
		print("Looking for conserved domains", quote=FALSE)
		domains <- Conserved_domains(ff1, strsplit(toString(pattern(aligned)), "")[[1]])
	}
	else
		domains <- list(domain_IDs="", domain_description="", domain_pos=rep(0, length(strsplit(toString(pattern(aligned)), "")[[1]])), domain_from="", domain_to="", e_value="", conflict="", all_domain_pos=0)
	
	#Result
	Visualize(strsplit(toString(pattern(aligned)), "")[[1]], strsplit(toString(subject(aligned)), "")[[1]], dist, protA_SecondaryStructure[[1]], domains)
	dist	
}

Distance <- function(seqA, seqB, W){
	#Collects physicochemical descriptors for protA from aaIndex database
	seqA_RADA880102 <- featureAAindex(seqA, aaindex.name="RADA880102")
	seqA_FAUJ880103 <- featureAAindex(seqA, aaindex.name="FAUJ880103")
	seqA_ZIMJ680104 <- featureAAindex(seqA, aaindex.name="ZIMJ680104")
	seqA_GRAR740102 <- featureAAindex(seqA, aaindex.name="GRAR740102")
	seqA_CRAJ730103 <- featureAAindex(seqA, aaindex.name="CRAJ730103")
	seqA_BURA740101 <- featureAAindex(seqA, aaindex.name="BURA740101")
	seqA_CHAM820102 <- featureAAindex(seqA, aaindex.name="CHAM820102")

	#Collects physicochemical descriptors for protB from aaIndex database
	seqB_RADA880102 <- featureAAindex(seqB, aaindex.name="RADA880102") 
	seqB_FAUJ880103 <- featureAAindex(seqB, aaindex.name="FAUJ880103")
	seqB_ZIMJ680104 <- featureAAindex(seqB, aaindex.name="ZIMJ680104")
	seqB_GRAR740102 <- featureAAindex(seqB, aaindex.name="GRAR740102")
	seqB_CRAJ730103 <- featureAAindex(seqB, aaindex.name="CRAJ730103")
	seqB_BURA740101 <- featureAAindex(seqB, aaindex.name="BURA740101")
	seqB_CHAM820102 <- featureAAindex(seqB, aaindex.name="CHAM820102")

	#Concatenates the vectors with physicochemical descriptors into two matrices corresponding to seqA and seqB
	As <- rbind(seqA_RADA880102, seqA_FAUJ880103, seqA_ZIMJ680104, seqA_GRAR740102, seqA_CRAJ730103, seqA_BURA740101, seqA_CHAM820102)
	Bs <- rbind(seqB_RADA880102, seqB_FAUJ880103, seqB_ZIMJ680104, seqB_GRAR740102, seqB_CRAJ730103, seqB_BURA740101, seqB_CHAM820102)

	#Calculates the distance matrix
	D <- As-Bs
	
	#Sets all NA elements in D to a default distance. These values corresponds to gaps in the aligned sequences
	D[is.na(D)] <- 0
	
	#Normalizes the distances
	D <- scale(D)
	D[is.nan(D)] <- 0

	D
}

SecStructure <- function(seq, seq_align){
	#The function predicts 2D-structure from seq. It returns a 2D-structure vector with inserted gaps matching those in seq_align
	Structure <- predictPROTEUS(seq, proteus2.organism="euk")
	seq_2D <- as.character((1:length(seq_align)))
	confidence <- (1:length(seq_align))
	position_nr <- 1
	gap <- 0
	
	#Insertion of gaps in the 2D-structure matching the gaps in seq_align
	for(position in seq_align){
		if(position != '-'){
			seq_2D[position_nr] <- Structure[[1]]$PSIPRED$SecondaryStructure[position_nr - gap]
		}
		else{
			seq_2D[position_nr] <- "-"
			gap <- gap + 1
		}
		position_nr <- position_nr + 1
	}
	list(seq_2D, Structure[[1]])
}

Weight <- function(Structure, Confidence, D){
	#The function takes a a string vector for the protein secondary structure and a vector with confidence scores as arguments. 
	#It returns a vector with weights for each position in the AA-sequence
	W <- matrix(1,nrow=dim(D)[1],ncol=dim(D)[2])
	position_nr <- 1
	gap <- 0

	for(position in Structure){
		if(position == 'E' && Confidence[position_nr - gap] > 5 ){
			W[,position_nr] <- 1.5
		}
		else if(position == 'H' && Confidence[position_nr - gap] > 5 ){
			W[1:5,position_nr] <- 1.2
			W[7,position_nr] <- 1.2
			
			#The weight for "Normalized frequency of alpha-helix"
			if(D[6,position_nr]>=0)
				W[6,position_nr] <- 1
			else
				W[6,position_nr] <- 1.5
		}
		else if(position == '-'){ #Penalty for insertion in mutated seq
			W[position_nr] <- 1
			gap <- gap + 1
		}
		position_nr <- position_nr + 1
	}
	
	#Calculates the distance scores
	property_distances <- D*W
	property_scores <- rowSums(property_distances)
	merged_prop_distances <- colSums(abs(D)*W)
	merged_score <- sum(merged_prop_distances)
	TotalScore <- sqrt(sum(D^2*W))
	
	list(property_distances=property_distances, property_scores=property_scores, merged_prop_distances=merged_prop_distances, merged_score=merged_score, TotalScore=TotalScore)
}

Conserved_domains <- function(ff, seq_align){
	#Performs a blast search using the program blastcl3
	blast.xml <- paste(system(paste("blastp -remote -db cdd -outfmt 5 -evalue 1e-2 -query ", ff), intern=TRUE), collapse="")
	#Creates an R treestructure from the XML-output given by blastcl3
	blast.tree <- xmlInternalTreeParse(blast.xml, asText=TRUE)
	#Collects the information of interest from the tree
	domain_IDs <- xpathSApply(blast.tree, "//Hit/Hit_id", xmlValue)
	domain_description <- xpathSApply(blast.tree, "//Hit/Hit_def", xmlValue)
	domain_from <- as.numeric(xpathSApply(blast.tree, "//Hit/Hit_hsps/Hsp/Hsp_query-from", xmlValue))
	domain_to <- as.numeric(xpathSApply(blast.tree, "//Hit/Hit_hsps/Hsp/Hsp_query-to", xmlValue))
	e_value <- as.numeric(xpathSApply(blast.tree, "//Hit/Hit_hsps/Hsp/Hsp_evalue", xmlValue))
	
	#Modifies the domain boundaries in accordance with the gaps in seq_align
	position_nr <- 1
	for(position in seq_align){
		if(position == '-'){
			i <- 1
			for(from in domain_from){
				if(from > position_nr){
					domain_from[i] <- domain_from[i] + 1
				}
				i <- i + 1
			}
			i <- 1
			for(to in domain_to){
				if(to > position_nr){
					domain_to[i] <- domain_to[i] + 1
				}
				i <- i + 1
			}
		}
		position_nr <- position_nr + 1
	}
	
	#Creates vectors indicating where the domains are located along the sequence and the domain predictions that overlap. These vectors are used by Visualize() and Domain_visualization() to visualize the domains
	domain_pos <- rep(0, length(seq_align))
	conflict <- rep(0, length(domain_from))
	all_domain_pos <- list()
	i <- 1
	for(from in domain_from){
		#Looks for overlaping domains. An overlap can occur in 4 different ways
		case1 <- domain_from[i]<=domain_from & (domain_to[i]>domain_from & domain_to[i]<domain_to)
		case2 <- (domain_from[i]>domain_from & domain_from[i]<domain_to) & domain_to[i]>=domain_to
		case3 <- (domain_from[i]>domain_from & domain_from[i]<domain_to) & (domain_to[i]>domain_from & domain_to[i]<domain_to)
		case4 <- domain_from[i]<=domain_from & domain_to[i]>=domain_to
		overlap <- case1 | case2 | case3 | case4
		
		#If overlap occurs, the domain with the lowest e-value is choosen for visualization
		if(length(which(overlap))>1){
			conflict[i] <- 1 #Indicates that there are conflicting domain predictions
			if(e_value[i] <= min(e_value[which(overlap)])){
				domain_pos[c(domain_from[i]:domain_to[i])] <- i
			}
		}
		else{
			domain_pos[c(domain_from[i]:domain_to[i])] <- i
		}
		
		#Creates a list containing one vector for each domain. This is used in the domain visualization to show the domains that are discarded in the the primary visualization because they overlap with a domain that has a lower e-value
		all_domain_pos[[i]] <- rep(0, length(seq_align))
		all_domain_pos[[i]][c(domain_from[i]:domain_to[i])] <- 1

		i <- i + 1
	}
	
	#Returns a list with the information of interest
	list(domain_IDs=domain_IDs, domain_description=domain_description, domain_pos=domain_pos, domain_from=domain_from, domain_to=domain_to, e_value=e_value, conflict=conflict, all_domain_pos=all_domain_pos)
}
