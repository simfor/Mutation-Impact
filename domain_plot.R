library(tkrplot)
library(BioSeqClass)
library(Biostrings)

Domain_visualization <- function(pattern, subject, dist, domains){
	tt <- tktoplevel()
	tkwm.title(tt,"Domain view")
	left <- tclVar(1)
	right <- tclVar(100)

	domain_plot <- function(){
		lleft <- as.numeric(tclvalue(left))
	      	rright <- as.numeric(tclvalue(right))
	      	x <- c(lleft:rright)
	      	layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,2))
	
		#Plots the sum of distances
		plot(x, dist$merged_prop_distances[x], type="l", axes=FALSE, ann=FALSE, xlim=range(x), ylim=c(0, max(dist$merged_prop_distances)), cex=0.7)
		axis(1, at=min(x):max(x), , cex.axis=0.6)
		
		#Adds the mutated sequence
		axis(1, at=min(x):max(x), lab=subject[x], line=2, lty=0, cex.axis=0.6)
		mtext("Mutated sequence", 1, cex=0.6, col="green", line=2)
	
		#Adds the reference sequence
		axis(1, at=min(x):max(x), lab=pattern[x], line=3, lty=0,  cex.axis=0.6) #tck=-0.1
		mtext("Ref. sequence", 1, cex=0.6, line=5, col="green")
		
		#Adds the conserved domains
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
		domain_plot()
		dev.off()
	}
	domain_info_button <- function(...){
		ReturnVal <- modalDialog("Domain info","Enter a sequence position. Additional information about domains on this position will be listed","")
  		if (ReturnVal=="ID_CANCEL")
    			return()
    		
    		info <- ""
    		i <- 1
    		for(dom_info in domains$all_domain_pos){
    			if(ReturnVal %in% which(dom_info==1)){
    				this_domain <- paste(domains$domain_IDs[[i]], domains$domain_description[[i]], sep=":\n")
    				info <- paste(info, this_domain, sep="\n\n")
    			}
    			i <- i+1
    		}
    		if(info=="")
    			info <- "No domains found on this position"
    		
		tkmessageBox(title="Domain info",message=info)
	}
	
	#Creates the widgets
	img <- tkrplot(tt, domain_plot, vscale=1.05, hscale=2.5) 
	s1 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=left, orient="horiz",label='left')
	s2 <- tkscale(tt, command=scroll, from=1, to=length(dist$merged_prop_distances), variable=right, orient="horiz",label='right')
	b1 <- tkbutton(tt, text='->', command=right_button)
	b2 <- tkbutton(tt, text='<-', command=left_button)
	b3 <- tkbutton(tt, text='Save as PDF', command=save_button)
	launchDlg.button <- tkbutton(tt,text="Domain info",command=domain_info_button)
	
	#Puts all the widgets in a grid
	tkgrid(img)
	tkgrid(s1)
	tkgrid(s2)
	tkgrid(b1)
	tkgrid(b2)
	tkgrid(b3, launchDlg.button, column=0)
	tkgrid.configure(launchDlg.button, sticky="e")
	tkgrid.configure(s1, sticky="ew")
	tkgrid.configure(s2, sticky="ew")
	tkfocus(img)
	
	#The modal dialog box which is activated by the domain info button
	modalDialog <- function(title,question,entryInit,entryWidth=20,returnValOnCancel="ID_CANCEL"){
		dlg <- tktoplevel()
		tkwm.deiconify(dlg)
		tkgrab.set(dlg)
		tkfocus(dlg)
		tkwm.title(dlg,title)
		textEntryVarTcl <- tclVar(paste(entryInit))
		textEntryWidget <- tkentry(dlg,width=paste(entryWidth),textvariable=textEntryVarTcl)
		tkgrid(tklabel(dlg,text="       "))
		tkgrid(tklabel(dlg,text=question),textEntryWidget)
		tkgrid(tklabel(dlg,text="       "))
		ReturnVal <- returnValOnCancel
		onOK <- function(){
			ReturnVal <<- tclvalue(textEntryVarTcl)
			tkgrab.release(dlg)
			tkdestroy(dlg)
			tkfocus(tt)
		}
		onCancel <- function(){
			ReturnVal <<- returnValOnCancel
			tkgrab.release(dlg)
			tkdestroy(dlg)
			tkfocus(tt)
		}
		OK.but     <-tkbutton(dlg,text="   OK   ",command=onOK)
		Cancel.but <-tkbutton(dlg,text=" Cancel ",command=onCancel)
		tkgrid(OK.but,Cancel.but)
		tkgrid(tklabel(dlg,text="    "))

		tkfocus(dlg)
		tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
		tkbind(textEntryWidget, "<Return>", onOK)
		tkwait.window(dlg)

		return(ReturnVal)
	}
}
