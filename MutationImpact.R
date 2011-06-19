library(BioSeqClass)
library(Biostrings)

Main <- function(ff1, ff2){
	protA <- readFASTA(ff1, strip.descs=TRUE)
	protB <- readFASTA(ff2, strip.descs=TRUE)
	
	#Retrieves the secondary structure for protA
	protA_SecondaryStructure <- SecStructure(protA[[1]]$seq)
	
	#Calculates the weight vector
	W = Weight(protA_SecondaryStructure[[1]]$PSIPRED$SecondaryStructure, protA_SecondaryStructure[[1]]$PSIPRED$ConfidenceScore)
	
	#Calculates the distance between the two sequences
	dist <- distance(protA[[1]]$seq, protB[[1]]$seq, W)
}

distance <- function(seqA, seqB, W){
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
	As <- rbind(protA_RADA880102, protA_FAUJ880103, protA_ZIMJ680104, protA_GRAR740102, protA_CRAJ730103, protA_BURA740101, protA_CHAM820102)
	Bs <- rbind(protB_RADA880102, protB_FAUJ880103, protB_ZIMJ680104, protB_GRAR740102, protB_CRAJ730103, protB_BURA740101, protB_CHAM820102)

	#Calculates the absolute distance matrix
	D <- abs(As-Bs)

	#Calculates the distance score
	dist <- sqrt(sum(W*t(D^2)))
}

SecStructure <- function(seq){
	Structure <- predictPROTEUS(seq, proteus2.organism="euk")
}

Weight <- function(Structure, Confidence){
	#The function takes a a string vector for the protein secondary structure and a vector with confidence scores as arguments. 
	#It returns a vector with weights for each position in the AA-sequence
	weight <- (1:length(Structure))
	position_nr <- 1

	for(position in Structure){
		if((position == 'H' || position == 'E') && Confidence[position_nr] > 5 ){
			weight[position_nr] <- 1.5
		}
		else{
			weight[position_nr] <- 1
		}
		position_nr <- position_nr + 1
	}
	weight
}
