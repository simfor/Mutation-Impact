tt <- tktoplevel()
left <- tclVar(1)
oldleft <- tclVar(1)
right <- tclVar(100)

Visualize_test <- function(pattern, subject, dist, left, right){
	lleft <- as.numeric(tclvalue(left))
      	rright <- as.numeric(tclvalue(right))
		
	length <- c(lleft:rright)
	plot_colors <- c("black","blue","red","forestgreen","yellow","green","magenta","burlywood3")
	plot_names <- c("Sum of all properties","Transfer free energy from octanol to water", "Normalized van der Waals volume", "Isoelectric point", "Polarity", "Normalized frequency of turn", "Normalized frequency of alpha-helix", "Free energy of solution in water")
	
#	par(mfrow=c(2,1))
#	layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,3) )
	
#	text(1,1,subject)

	par(fig=c(0,1,0.1,1), cex=0.8, cex.axis=0.9)
	plot(length, dist$merged_prop_distances[lleft:rright], type="l", col=plot_colors[1], axes=FALSE, ann=FALSE)
	
	axis(1, at=1:length(dist$merged_prop_distances))
	axis(2, at=1:max(dist$merged_prop_distances))
	title(xlab="Position", col.lab=rgb(0,0.5,0))
	title(ylab="Distance", col.lab=rgb(0,0.5,0))

	lines(length, dist$property_distances[1,][lleft:rright], lty=2, col=plot_colors[2])
	lines(length, dist$property_distances[2,][lleft:rright], lty=2, col=plot_colors[3])
	lines(length, dist$property_distances[3,][lleft:rright], lty=2, col=plot_colors[4])
	lines(length, dist$property_distances[4,][lleft:rright], lty=2, col=plot_colors[5])
	lines(length, dist$property_distances[5,][lleft:rright], lty=2, col=plot_colors[6])
	lines(length, dist$property_distances[6,][lleft:rright], lty=2, col=plot_colors[7])
	lines(length, dist$property_distances[7,][lleft:rright], lty=2, col=plot_colors[8])

	x_legend2 <- length(length)/1.5

	legend(1, 20, plot_names, cex=0.8, col=plot_colors, title="Properties", lty=c(1,2,2,2,2,2,2,2), bty="n") # "topleft"
	legend(x_legend2, 20, c(dist$merged_score, dist$property_scores), title="Total", cex=0.8, bty="n") # "topright"

	par(fig=c(0,1,0.01,0.05), cex=0.4, new=TRUE)
#	textplot(c(pattern,subject), halign="center")
	axis(1, at=1:length(length), lab=subject[lleft:rright], lty=0)
#	title(sub="Mutated sequence")
	mtext("Mutated sequence", 1, cex=0.6)

	par(fig=c(0,1,0,0.1), cex=0.4, new=TRUE)
	axis(1, at=1:length(length), lab=pattern[lleft:rright], tck=-0.1)
#	title(xlab="Ref. sequence")
	mtext("Ref. sequence", 1, cex=0.6)
}

