library(tkrplot)
library(BioSeqClass)
library(Biostrings)

#pattern <- strsplit(toString(pattern(align_O00198_insert)), "")[[1]]
#subject <- strsplit(toString(subject(align_O00198_insert)), "")[[1]]

Visualize <- function(pattern, subject, dist, sec_structure, domains){

	tt <- tktoplevel()
	tkwm.title(tt,"Mutation impact")
	left <- tclVar(1)
	right <- tclVar(100)

	Visualize_plot <- function(){ 
		#The new left and right margins when scrolling
		lleft <- as.numeric(tclvalue(left))
	      	rright <- as.numeric(tclvalue(right))
		
		x <- c(lleft:rright)
		plot_colors <- c("black","blue","red","forestgreen","yellow","green","magenta","burlywood3")
		plot_names <- c("Sum of all properties","Transfer free energy from octanol to water", "Normalized van der Waals volume", "Isoelectric point", "Polarity", "Normalized frequency of turn", "Normalized frequency of alpha-helix", "Free energy of solution in water")
		
		#The graphical parameters for the distance plot
		par(cex=0.7, cex.axis=0.9, mar=c(12, 4.1, 4.1, 2.1))
		#Plots the sum of distances
		plot(x, dist$merged_prop_distances[x], type="l", col=plot_colors[1], axes=FALSE, ann=FALSE, xlim=range(x), ylim=c(min(dist$property_distances), max(dist$merged_prop_distances))) 
	
		axis(1, at=min(x):max(x))
		axis(2, at=round(min(dist$property_distances),digits=2):round(max(dist$merged_prop_distances), digits=2))
#		title(xlab="Position", col.lab=rgb(0,0.5,0))
		title(ylab="Distance", col.lab=rgb(0,0.5,0))
		
		#Adds the different distances to the plot
		lines(x, dist$property_distances[1,][x], lty=2, col=plot_colors[2])
		lines(x, dist$property_distances[2,][x], lty=2, col=plot_colors[3])
		lines(x, dist$property_distances[3,][x], lty=2, col=plot_colors[4])
		lines(x, dist$property_distances[4,][x], lty=2, col=plot_colors[5])
		lines(x, dist$property_distances[5,][x], lty=2, col=plot_colors[6])
		lines(x, dist$property_distances[6,][x], lty=2, col=plot_colors[7])
		lines(x, dist$property_distances[7,][x], lty=2, col=plot_colors[8])

		#The distance between the two legends
		x_legend2 <- length(x)/3
		
		legend(lleft, max(dist$merged_prop_distances), plot_names, cex=0.8, col=plot_colors, title="Properties", lty=c(1,2,2,2,2,2,2,2), bty="n") #"topleft"
		legend(lleft + x_legend2, max(dist$merged_prop_distances), c(dist$merged_score, dist$property_scores), title="Total", cex=0.8, bty="n") #"topright"

		#The graphical parameters for the mutated sequence
#		par(fig=c(0,1,0.71,0.75), cex=0.6, new=TRUE)
		#Adds the mutated sequence
		axis(1, at=min(x):max(x), lab=subject[x], line=2, lty=0)
		mtext("Mutated sequence", 1, cex=0.6, col="green", line=2)
	
		#The graphical parameters for the reference sequence
#		par(fig=c(0,1,0.7,0.8), cex=0.6, new=TRUE)
		#Adds the reference sequence
		axis(1, at=min(x):max(x), lab=pattern[x], line=3, lty=0) #tck=-0.1
		mtext("Ref. sequence", 1, cex=0.6, line=5, col="green")
		
		#The graphical parameters for the 2D structure
#		par(fig=c(0,1,0,0.1), cex=0.6, new=TRUE)#col.lab="blue"
		#Adds the 2D structure
		axis(1, at=min(x):max(x), lab=sec_structure[x], line=6, col.lab="blue", lty=0) #tck=-0.1
		mtext("2D structure", 1, cex=0.6, line=8, col="red")
		
		#The graphical parameters for the conserved domains
		par(fig=c(0,1,0,0.1), new=TRUE)
		#Adds the conserved domains
		domains_in_view <- domains$domain_pos[x]
		for(curr_domain in unique(domains_in_view)){
			if(curr_domain != 0){
				if(domains$conflict[curr_domain] == 0){
					rect(x[min(which(domains_in_view==curr_domain))], 0,x[max(which(domains_in_view==curr_domain))], min(dist$merged_prop_distances[which(dist$merged_prop_distances!=0)]), col="green")
				}
				else{
					rect(x[min(which(domains_in_view==curr_domain))], 0,x[max(which(domains_in_view==curr_domain))], min(dist$merged_prop_distances[which(dist$merged_prop_distances!=0)]), col="red")
				}
				text(x[median(which(domains_in_view==curr_domain))], min(dist$merged_prop_distances[which(dist$merged_prop_distances!=0)])/2, labels=paste("Domain ID:", domains$domain_IDs[curr_domain], "e-value:", domains$e_value[curr_domain], sep=" "))
			}
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
	domain_button <- function(...){
		domain_visualization(pattern, subject, dist, domains)
	}

	#Creates the widgets
	img <- tkrplot(tt, Visualize_plot, vscale=1.05, hscale=2.5) 
	s1 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=left, orient="horiz",label='left')
	s2 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=right, orient="horiz",label='right')
	b1 <- tkbutton(tt, text='->', command=right_button)
	b2 <- tkbutton(tt, text='<-', command=left_button)
	b3 <- tkbutton(tt, text='Save as PDF', command=save_button)
	b4 <- tkbutton(tt, text='Conserved domains view', command=domain_button)

	#Puts all the widgets in a grid
	tkgrid(img)
	tkgrid(s1)
	tkgrid(s2)
	tkgrid(b1)
	tkgrid(b2)
	tkgrid(b3, b4, column=0)
	tkgrid.configure(b4, sticky="e")
	tkgrid.configure(s1, sticky="ew") 
	tkgrid.configure(s2, sticky="ew")
	tkfocus(img)
}
