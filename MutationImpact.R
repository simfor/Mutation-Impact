library(BioSeqClass)
library(Biostrings)

Main <- function(){

}

distance <- function(ff1, ff2){
	protA <- readFASTA(ff1, strip.descs=TRUE)
	protB <- readFASTA(ff2, strip.descs=TRUE)

	#Collects physicochemical descriptors for protA from aaIndex database
	protA_RADA880102 <- featureAAindex(protA[[1]]$seq, aaindex.name="RADA880102") 
	protA_FAUJ880103 <- featureAAindex(protA[[1]]$seq, aaindex.name="FAUJ880103")
	protA_ZIMJ680104 <- featureAAindex(protA[[1]]$seq, aaindex.name="ZIMJ680104")
	protA_GRAR740102 <- featureAAindex(protA[[1]]$seq, aaindex.name="GRAR740102")
	protA_CRAJ730103 <- featureAAindex(protA[[1]]$seq, aaindex.name="CRAJ730103")
	protA_BURA740101 <- featureAAindex(protA[[1]]$seq, aaindex.name="BURA740101")
	protA_CHAM820102 <- featureAAindex(protA[[1]]$seq, aaindex.name="CHAM820102")

	#Collects physicochemical descriptors for protB from aaIndex database
	protB_RADA880102 <- featureAAindex(protB[[1]]$seq, aaindex.name="RADA880102") 
	protB_FAUJ880103 <- featureAAindex(protB[[1]]$seq, aaindex.name="FAUJ880103")
	protB_ZIMJ680104 <- featureAAindex(protB[[1]]$seq, aaindex.name="ZIMJ680104")
	protB_GRAR740102 <- featureAAindex(protB[[1]]$seq, aaindex.name="GRAR740102")
	protB_CRAJ730103 <- featureAAindex(protB[[1]]$seq, aaindex.name="CRAJ730103")
	protB_BURA740101 <- featureAAindex(protB[[1]]$seq, aaindex.name="BURA740101")
	protB_CHAM820102 <- featureAAindex(protB[[1]]$seq, aaindex.name="CHAM820102")

	#Concatenates the vectors with physicochemical descriptors into two matrices corresponding to protA and protB
	As <- rbind(protA_RADA880102, protA_FAUJ880103, protA_ZIMJ680104, protA_GRAR740102, protA_CRAJ730103, protA_BURA740101, protA_CHAM820102)
	Bs <- rbind(protB_RADA880102, protB_FAUJ880103, protB_ZIMJ680104, protB_GRAR740102, protB_CRAJ730103, protB_BURA740101, protB_CHAM820102)

	#Calculates the absolute distance matrix
	D <- abs(As-Bs)

	#Calculates the distance score
	dist <- sqrt(sum(D^2))
}

SecStructure <- function(seq){
	Structure <- predictPROTEUS(seq, proteus2.organism="euk")
}

Weight <- function(Structure){
	#The function takes a a string vector for the protein secondary structure as argument and returns a vector with weights for each position in the AA-sequence
	weight <- (1:length(Structure))
	position_nr <- 1

	for(position in Structure){
		if(position == 'C'){
			weight[position_nr] <- 1
		}
		else{
			weight[position_nr] <- 1.5
		}
		position_nr <- position_nr + 1
	}
	weight
}
