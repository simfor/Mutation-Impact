library(tkrplot)
library(BioSeqClass)
library(Biostrings)

domain_visualization <- function(pattern, subject, dist, domains){
	tt <- tktoplevel()
	left <- tclVar(1)
	right <- tclVar(100)

	domain_plot <- function(){
		lleft <- as.numeric(tclvalue(left))
	      	rright <- as.numeric(tclvalue(right))
	      	x <- c(lleft:rright)
	      	#par(mfrow=c(2,1))
	      	layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,2))
	
		#The graphical parameters for the distance plot
		#par(fig=c(0,1,0.7,1), cex=0.7, cex.axis=0.9 )#mar=c(3.0, 3.0, 1.5, 1.5)
		#Plots the sum of distances
		plot(x, dist$merged_prop_distances[x], type="l", axes=FALSE, ann=FALSE, xlim=range(x), ylim=c(0, max(dist$merged_prop_distances)), cex=0.7)
		axis(1, at=min(x):max(x), , cex.axis=0.6)
		#axis(2, at=0:max(dist$merged_prop_distances))

		#The graphical parameters for the mutated sequence
		#par(fig=c(0,1,0.71,0.75), cex=0.6, new=TRUE)
		#Adds the mutated sequence
		axis(1, at=min(x):max(x), lab=subject[x], line=2, lty=0, cex.axis=0.6)
		mtext("Mutated sequence", 1, cex=0.6, col="green", line=2)
	
		#The graphical parameters for the reference sequence
		#par(fig=c(0,1,0.7,0.8), cex=0.6, new=TRUE)
		#Adds the reference sequence
		axis(1, at=min(x):max(x), lab=pattern[x], line=3, lty=0,  cex.axis=0.6) #tck=-0.1
		mtext("Ref. sequence", 1, cex=0.6, line=5, col="green")
		
		#domains_in_view <- domains$domain_pos[x]
		plot(rep(10,length(x)), axes=FALSE, ann=FALSE, col="white", xlim=range(x))
		i <- 1
		for(domain in domains$all_domain_pos){
			if(length(which(domain[x]==1)) > 0){
				rect(x[min(which(domain[x]==1))], 15-i,x[max(which(domain[x]==1))], 16-i, col="green", xpd=NA)
				text(x[median(which(domain[x]==1))], 15.5-i, labels=paste("Domain ID:", domains$domain_IDs[i], "e-value:", domains$e_value[i], sep=" "), cex=0.8, xpd=NA)
			}
			i <- i + 1
		}
	}
	
	scroll <- function(...){
		tkrreplot(img)
	}
	right_button <- function(...){
		tclvalue(left) <- as.character(as.numeric(tclvalue(left)) + 10)
		tclvalue(right) <- as.character(as.numeric(tclvalue(right)) + 10)
		tkrreplot(img)
	}
	
	left_button <- function(...){
		tclvalue(left) <- as.character(as.numeric(tclvalue(left)) - 10)
		tclvalue(right) <- as.character(as.numeric(tclvalue(right)) - 10)
		tkrreplot(img)
	}
	save_button <- function(...){
		fileName<-tclvalue(tkgetSaveFile())
		if (!nchar(fileName))
		    tkmessageBox(message="No file was selected!")
		else
		    tkmessageBox(message=paste("Graph saved as: ",fileName))

		pdf(fileName)
		Visualize_plot()
		dev.off()
	}
	
	img <- tkrplot(tt, domain_plot, vscale=1.05, hscale=2.5) 
	s1 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=left, orient="horiz",label='left')
	s2 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=right, orient="horiz",label='right')
	b1 <- tkbutton(tt, text='->', command=right_button)
	b2 <- tkbutton(tt, text='<-', command=left_button)
	b3 <- tkbutton(tt, text='Save as PDF', command=save_button)
	
	tkgrid(img)
	tkgrid(s1)
	tkgrid(s2)
	tkgrid(b1)
	tkgrid(b2)
	tkgrid(b3)
	tkgrid.configure(s1, sticky="ew")
	tkgrid.configure(s2, sticky="ew")
	tkfocus(img)
}
