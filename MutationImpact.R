library(BioSeqClass)
library(Biostrings)
#library(gplots)
library(tkrplot)

Main <- function(ff1, ff2){
	protA <- readFASTA(ff1, strip.descs=TRUE)
	protB <- readFASTA(ff2, strip.descs=TRUE)
	
	#Aligns the two sequences
	print("Aligning sequences")
	aligned <- Align(protA[[1]]$seq, protB[[1]]$seq)
	
	#Retrieves the secondary structure for protA
	print("Retrieving secondary structure")
	protA_SecondaryStructure <- SecStructure(protA[[1]]$seq, strsplit(toString(pattern(aligned)), "")[[1]] )
	
	#Calculates the weight vector
	W <- Weight(protA_SecondaryStructure[[1]], protA_SecondaryStructure[[2]]$PSIPRED$ConfidenceScore)
	
	#Calculates the distance between the two sequences
	#dist <- Distance(protA[[1]]$seq, protB[[1]]$seq, W)
	dist <- Distance(toString(pattern(aligned)), toString(subject(aligned)), W)
	
	#Result
	Visualize(strsplit(toString(pattern(aligned)), "")[[1]], strsplit(toString(subject(aligned)), "")[[1]], dist)
	dist
	
	#print(aligned)
	#paste("The distance between the sequences: ", dist)
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

	#Calculates the absolute distance matrix
	D <- abs(As-Bs)
	
	#Sets all NA elements in D to a default distance. These values corresponds to gaps in the aligned sequences
	D[is.na(D)] <- 1

	#Calculates the distance scores
	property_distances <- t(t(D)*W)
	property_scores <- rowSums(property_distances) #sqrt(rowSums(t(t(D^2)*W)))
	merged_prop_distances <- colSums(property_distances)
	merged_score <- sum(merged_prop_distances) #sqrt(sum(t(t(merged_prop_distances^2))))
	TotalScore <- sqrt(sum(t(D^2)*W))
	
	list(property_distances=property_distances, property_scores=property_scores, merged_prop_distances=merged_prop_distances, merged_score=merged_score, TotalScore=TotalScore)
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

Weight <- function(Structure, Confidence){
	#The function takes a a string vector for the protein secondary structure and a vector with confidence scores as arguments. 
	#It returns a vector with weights for each position in the AA-sequence
	weight <- (1:length(Structure))
	position_nr <- 1
	gap <- 0

	for(position in Structure){
		if((position == 'H' || position == 'E') && Confidence[position_nr - gap] > 5 ){
			weight[position_nr] <- 1.5
		}
		else if(position == '-'){ #Penalty for insertion in mutated seq
			weight[position_nr] <- 1.3
			gap <- gap + 1
		}
		else{
			weight[position_nr] <- 1
		}
		position_nr <- position_nr + 1
	}
	weight
}

Align <- function(seqA, seqB){
	aligned <- pairwiseAlignment(pattern = seqA, subject = seqB, type="global", substitutionMatrix = "BLOSUM50", gapOpening = -5, gapExtension = -1)
	
#	seqA_align <- strsplit(toString(pattern(aligned)), "")
#	seqB_align <- strsplit(toString(subject(aligned)), "")
	
#	list()
}
