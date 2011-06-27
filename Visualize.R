pattern <- strsplit(toString(pattern(align_O00198_insert)), "")[[1]]
subject <- strsplit(toString(subject(align_O00198_insert)), "")[[1]]

tt <- tktoplevel()
left <- tclVar(1)
oldleft <- tclVar(1)
right <- tclVar(100)

Visualize_test <- function(){
	lleft <- as.numeric(tclvalue(left))
      	rright <- as.numeric(tclvalue(right))
		
	x <- seq(lleft,rright,by=1)
	plot_colors <- c("black","blue","red","forestgreen","yellow","green","magenta","burlywood3")
	plot_names <- c("Sum of all properties","Transfer free energy from octanol to water", "Normalized van der Waals volume", "Isoelectric point", "Polarity", "Normalized frequency of turn", "Normalized frequency of alpha-helix", "Free energy of solution in water")
	
#	par(mfrow=c(2,1))
#	layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,3) )
	
#	text(1,1,subject)

	par(fig=c(0,1,0.1,1), cex=0.7, cex.axis=0.9)
	plot(x, dist$merged_prop_distances[x], type="l", col=plot_colors[1], axes=FALSE, ann=FALSE)
	
	axis(1, at=1:length(x))
	axis(2, at=1:max(dist$merged_prop_distances))
	title(xlab="Position", col.lab=rgb(0,0.5,0))
	title(ylab="Distance", col.lab=rgb(0,0.5,0))

	lines(x, dist$property_distances[1,][x], lty=2, col=plot_colors[2])
	lines(x, dist$property_distances[2,][x], lty=2, col=plot_colors[3])
	lines(x, dist$property_distances[3,][x], lty=2, col=plot_colors[4])
	lines(x, dist$property_distances[4,][x], lty=2, col=plot_colors[5])
	lines(x, dist$property_distances[5,][x], lty=2, col=plot_colors[6])
	lines(x, dist$property_distances[6,][x], lty=2, col=plot_colors[7])
	lines(x, dist$property_distances[7,][x], lty=2, col=plot_colors[8])

	x_legend2 <- length(x)/1.5

	legend(1, 20, plot_names, cex=0.8, col=plot_colors, title="Properties", lty=c(1,2,2,2,2,2,2,2), bty="n") # "topleft"
	legend(x_legend2, 20, c(dist$merged_score, dist$property_scores), title="Total", cex=0.8, bty="n") # "topright"

	par(fig=c(0,1,0.01,0.05), cex=0.6, new=TRUE)
#	textplot(c(pattern,subject), halign="center")
	axis(1, at=1:length(x), lab=subject[x], lty=0)
#	title(sub="Mutated sequence")
	mtext("Mutated sequence", 1, cex=0.6)

	par(fig=c(0,1,0,0.1), cex=0.6, new=TRUE)
	axis(1, at=1:length(x), lab=pattern[x], lty=0) #tck=-0.1
#	title(xlab="Ref. sequence")
	mtext("Ref. sequence", 1, cex=0.6)
}

img <- tkrplot(tt, Visualize_test)

f2 <- function(...){
        ol <- as.numeric(tclvalue(oldleft))
        tclvalue(oldleft) <- tclvalue(left)
        r <- as.numeric(tclvalue(right))
        tclvalue(right) <- as.character(r + as.numeric(...) - ol)
        tkrreplot(img)
}

f3 <- function(...){
        tkrreplot(img)
}

f4 <- function(...){
        ol <- as.numeric(tclvalue(oldleft))
        tclvalue(left) <- as.character(ol+10)
        tclvalue(oldleft) <- as.character(ol+10)
        r <- as.numeric(tclvalue(right))
        tclvalue(right) <- as.character(r+10)
        tkrreplot(img)
}

s1 <- tkscale(tt, command=f2, from=1, to=length(dist$merged_prop_distances), variable=left, orient="horiz",label='left')
s2 <- tkscale(tt, command=f3, from=1, to=length(dist$merged_prop_distances), variable=right, orient="horiz",label='right')
b1 <- tkbutton(tt, text='->', command=f4)

tkpack(img,s1,s2,b1)