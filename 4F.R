##############################################################################################################################
#                                                                                                                            #
# 4F is a simple R script to perform sliding window-based analyses on four-fold degenerate codons of a mitochondrial genome. #
#                                                                                                                            #
# Copyright (C) 2015 Federico Plazzi                                                                                         #
#                                                                                                                            #
# This program is free software: you can redistribute it and/or modify                                                       #
# it under the terms of the GNU General Public License as published by                                                       #
# the Free Software Foundation, either version 3 of the License, or                                                          #
# (at your option) any later version.                                                                                        #
#                                                                                                                            #
# This program is distributed in the hope that it will be useful,                                                            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                              #
# GNU General Public License for more details.                                                                               #
#                                                                                                                            #
# You should have received a copy of the GNU General Public License                                                          #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                                      #
#                                                                                                                            #
##############################################################################################################################

#4F version: 1.0

#Loading the orcutt package, which is needed for its function cochrane.orcutt().
library(orcutt)

#Setting defaults.
four.fold.codon <- c("TC","CT","CC","CG","AC","AG","GT","GC","GG")
rc.four.fold.codon <- c("GA","AG","GG","CG","GT","CT","AC","GC","CC")
infile <- "infile"
outfile <- "outfile"
plotfile <- "plot.pdf"
wsize <- 200
wstep <- 50
gene <- "cox1"

#Reading the input files and looking to all genomes by default.
for (i in 1:length(scan("../4F.conf",what=character(),quiet=TRUE))) {
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###infile###",infile <- scan("../4F.conf",what=character(),quiet=TRUE)[i+1],NA)
	}
genome <- character()
for (i in 1:length(scan(paste("../",infile,sep=""),what=character(),quiet=TRUE,sep=">"))) {
	ifelse(i %% 3 == 2,genome <- c(genome,scan(paste("../",infile,sep=""),what=character(),quiet=TRUE,sep=">")[i]),NA)
	}

#Reading genome list.
for (i in 1:(length(scan("../4F.conf",what=character(),quiet=TRUE))-2)) {
	if (scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###genome###") {
		end.genomes <- length(scan("../4F.conf",what=character(),quiet=TRUE))
		for (j in (i+2):length(scan("../4F.conf",what=character(),quiet=TRUE))) {
			if (strsplit(scan("../4F.conf",what=character(),quiet=TRUE)[j],split="")[[1]][1] == "#") ifelse(j-1 > end.genomes,NA,end.genomes <- j-1)
			}
		genome <- scan("../4F.conf",what=character(),quiet=TRUE)[(i+1):end.genomes]
		}
	}
num.genomes <- length(genome)

#Reading gene list.
for (i in 1:(length(scan("../4F.conf",what=character(),quiet=TRUE))-2)) {
	if (scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###gene###") {
		end.genes <- length(scan("../4F.conf",what=character(),quiet=TRUE))
		for (j in (i+2):length(scan("../4F.conf",what=character(),quiet=TRUE))) {
			if (strsplit(scan("../4F.conf",what=character(),quiet=TRUE)[j],split="")[[1]][1] == "#") ifelse(j-1 > end.genes,NA,end.genes <- j-1)
			}
		gene <- scan("../4F.conf",what=character(),quiet=TRUE)[(i+1):end.genes]
		}
	}
ifelse(length(gene) == 1,gene <- rep(gene,length(genome)),NA)

#Defining final figure shape for more than genome.
ncol <- 1
nrow <- 1
if (num.genomes > 1) {
	ifelse(num.genomes %% 2 == 0,upper.bound <- num.genomes/2,upper.bound <- num.genomes %/% 3)
	ncols <- numeric()
	nrows <- numeric()
	mods <- numeric()
	differences <- numeric()
	for (i in 1:upper.bound) {
		ncols <- c(ncols,i)
		ifelse(num.genomes %% i == 0,nrows <- c(nrows,num.genomes/i),nrows <- c(nrows,num.genomes/i+1))
		mods <- c(mods,num.genomes %% i)
		ifelse(num.genomes %% i == 0,differences <- c(differences,abs(i-num.genomes/i)),differences <- c(differences,abs(i-num.genomes/i-1)))
		}
	table.format <- data.frame(ncols=ncols,nrows=nrows,mods=mods,differences=differences)
	if (length(rownames(table.format[table.format$mods == 0,])) == 0) table.format <- table.format[table.format$mods == 0,]
	table.format <- table.format[table.format$mods == max(table.format$mods),]
	table.format <- table.format[table.format$differences == min(table.format$differences),]
	ncol <- table.format[1,1]
	nrow <- table.format[1,2]
	}

#Reading user-defined parameters.
for (i in 1:length(scan("../4F.conf",what=character(),quiet=TRUE))) {
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###outfile###",outfile <- scan("../4F.conf",what=character(),quiet=TRUE)[i+1],NA)
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###plot###",plotfile <- scan("../4F.conf",what=character(),quiet=TRUE)[i+1],NA)
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###ncol###",ncol <- as.numeric(scan("../4F.conf",what=character(),quiet=TRUE)[i+1]),NA)
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###nrow###",nrow <- as.numeric(scan("../4F.conf",what=character(),quiet=TRUE)[i+1]),NA)
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###wsize###",wsize <- as.numeric(scan("../4F.conf",what=character(),quiet=TRUE)[i+1]),NA)
	ifelse(scan("../4F.conf",what=character(),quiet=TRUE)[i] == "###wstep###",wstep <- as.numeric(scan("../4F.conf",what=character(),quiet=TRUE)[i+1]),NA)
	}

par(mfrow=c(nrow,ncol*2))

for (g in 1:num.genomes) {
	#Starting analyzing genomes: open the sequence file and the annotation table.
	for (i in 1:length(scan(paste("../",infile,sep=""),what=character(),quiet=TRUE,sep=">"))) ifelse(scan(paste("../",infile,sep=""),what=character(),quiet=TRUE,sep=">")[i] == genome[g],mtDNA <- scan(paste("../",infile,sep=""),what=character(),quiet=TRUE,sep=">")[i+1],NA)
	annotation <- read.table(paste("../",genome[g],"_annotation.conf",sep=""),quote="",header=TRUE,row.names=1)
	
	#Building a vector for the current genome.
	nuc <- strsplit(mtDNA,split="")[[1]]
	nuc.chosen <- numeric()
	for (i in 1:length(nuc)) ifelse(i %% wstep == 0,nuc.chosen <- c(nuc.chosen,i),NA)
	
	#Identifying four-fold degenerated sites for the current genome.
	four.fold.sites <- numeric()
	four.fold.nuc <- character()
	for (i in 1:length(row.names(annotation))) {
		if (annotation[i,2] > annotation[i,1]) {
			for (j in 1:(annotation[i,2]-annotation[i,1]+1)) {
				if (j %% 3 == 0 && paste(nuc[j+annotation[i,1]-3],nuc[j+annotation[i,1]-2],sep="") %in% four.fold.codon) {
					four.fold.sites <- c(four.fold.sites,j+annotation[i,1]-1)
					four.fold.nuc <- c(four.fold.nuc,nuc[j+annotation[i,1]-1])
					}
				}
			}
		else {
			for (j in 1:(annotation[i,1]-annotation[i,2]+1)) {
				if (j %% 3 == 0 && paste(nuc[annotation[i,1]-j+2],nuc[annotation[i,1]-j+3],sep="") %in% rc.four.fold.codon) {
					four.fold.sites <- c(four.fold.sites,annotation[i,1]-j+1)
					four.fold.nuc <- c(four.fold.nuc,nuc[annotation[i,1]-j+1])
					}
				}
			}
		}
	
	four.fold.data.frame <- data.frame(sites=four.fold.sites,nuc=four.fold.nuc)
	
	#Computing nucleotide frequencies at four-fold degenerated sites.
	As <- numeric()
	Cs <- numeric()
	Gs <- numeric()
	Ts <- numeric()
	ATskews <- numeric()
	nuc.rescaled <- numeric()
	for (i in 1:length(nuc.chosen)) {
		if (nuc.chosen[i] <= (wsize%/%2)) {
			current.data.frame <- four.fold.data.frame[four.fold.data.frame$sites <= (nuc.chosen[i]+wsize%/%2),]
			current.left.tail.data.frame <- four.fold.data.frame[four.fold.data.frame$sites >= (length(nuc)-wsize%/%2+nuc.chosen[i]+1),]
			As <- c(As,100*(length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "A"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"]))/length(c(as.vector(current.left.tail.data.frame$nuc),as.vector(current.data.frame$nuc)))))
			Cs <- c(Cs,100*(length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "C"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "C"]))/length(c(as.vector(current.left.tail.data.frame$nuc),as.vector(current.data.frame$nuc)))))
			Gs <- c(Gs,100*(length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "G"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "G"]))/length(c(as.vector(current.left.tail.data.frame$nuc),as.vector(current.data.frame$nuc)))))
			Ts <- c(Ts,100*(length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "T"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"]))/length(c(as.vector(current.left.tail.data.frame$nuc),as.vector(current.data.frame$nuc)))))
			ATskews <- c(ATskews,(length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "A"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"]))-length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "T"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"])))/length(c(as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "A"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"],as.vector(current.left.tail.data.frame$nuc)[as.vector(current.left.tail.data.frame$nuc) == "T"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"])))
			}
		else if (nuc.chosen[i] > (wsize%/%2) && nuc.chosen[i] < (length(nuc)-wsize%/%2+1)) {
			current.data.frame <- four.fold.data.frame[four.fold.data.frame$sites >= (nuc.chosen[i]-wsize%/%2),]
			current.data.frame <- current.data.frame[current.data.frame$sites <= (nuc.chosen[i]+wsize%/%2),]
			As <- c(As,100*(length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"])/length(as.vector(current.data.frame$nuc))))
			Cs <- c(Cs,100*(length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "C"])/length(as.vector(current.data.frame$nuc))))
			Gs <- c(Gs,100*(length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "G"])/length(as.vector(current.data.frame$nuc))))
			Ts <- c(Ts,100*(length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"])/length(as.vector(current.data.frame$nuc))))
			ATskews <- c(ATskews,(length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"])-length(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"]))/length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"])))
			}
		else {
			current.data.frame <- four.fold.data.frame[four.fold.data.frame$sites >= (nuc.chosen[i]-wsize%/%2),]
			current.right.tail.data.frame <- four.fold.data.frame[four.fold.data.frame$sites <= (wsize%/%2-length(nuc)+nuc.chosen[i]),]
			As <- c(As,100*(length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "A"]))/length(c(as.vector(current.data.frame$nuc),as.vector(current.right.tail.data.frame$nuc)))))
			Cs <- c(Cs,100*(length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "C"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "C"]))/length(c(as.vector(current.data.frame$nuc),as.vector(current.right.tail.data.frame$nuc)))))
			Gs <- c(Gs,100*(length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "G"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "G"]))/length(c(as.vector(current.data.frame$nuc),as.vector(current.right.tail.data.frame$nuc)))))
			Ts <- c(Ts,100*(length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "T"]))/length(c(as.vector(current.data.frame$nuc),as.vector(current.right.tail.data.frame$nuc)))))
			ATskews <- c(ATskews,(length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "A"]))-length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "T"])))/length(c(as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "A"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "A"],as.vector(current.data.frame$nuc)[as.vector(current.data.frame$nuc) == "T"],as.vector(current.right.tail.data.frame$nuc)[as.vector(current.right.tail.data.frame$nuc) == "T"])))
			}
		ifelse(nuc.chosen[i] >= annotation[gene[g],1],nuc.rescaled <- c(nuc.rescaled,nuc.chosen[i]-annotation[gene[g],1]+1),nuc.rescaled <- c(nuc.rescaled,nuc.chosen[i]+length(nuc)-annotation[gene[g],1]+1))
		}

	#Plotting single-nucleotide percentages.
	plot(nuc.rescaled,As,type="p",main=paste(genome[g]," single nucleotides",sep=""),sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),xlim=c(1,length(nuc)),ylim=c(0,100),xlab="mtDNA",ylab="%",pch=21,col="#00FF00",bg="#00FF00")
	lines(nuc.rescaled,Cs,type="p",pch=21,col="#0000FF",bg="#0000FF")
	lines(nuc.rescaled,Gs,type="p",pch=21,col="#000000",bg="#000000")
	lines(nuc.rescaled,Ts,type="p",pch=21,col="#FF0000",bg="#FF0000")
	
	#Computing linear models.
	intercept <- numeric()
	slope <- numeric()
	r2 <- numeric()
	adj.r2 <- numeric()
	p.value <- numeric()
	auto.p.value <- numeric()
	rho <- numeric()
	
	regA <- lm(As~nuc.rescaled)
	abline(regA,col="#00FF00")
	intercept <- c(intercept,summary(regA)$coefficients[1,1])
	slope <- c(slope,summary(regA)$coefficients[2,1])
	r2 <- c(r2,summary(regA)$r.squared)
	adj.r2 <- c(adj.r2,summary(regA)$adj.r.squared)
	p.value <- c(p.value,summary(regA)$coefficients[2,4])
	auto.regA <- cochrane.orcutt(regA)
	auto.p.value <- c(auto.p.value,auto.regA$Cochrane.Orcutt$coefficients[2,4])
	rho <- c(rho,auto.regA$rho)
	regC <- lm(Cs~nuc.rescaled)
	abline(regC,col="#0000FF")
	intercept <- c(intercept,summary(regC)$coefficients[1,1])
	slope <- c(slope,summary(regC)$coefficients[2,1])
	r2 <- c(r2,summary(regC)$r.squared)
	adj.r2 <- c(adj.r2,summary(regC)$adj.r.squared)
	p.value <- c(p.value,summary(regC)$coefficients[2,4])
	auto.regC <- cochrane.orcutt(regC)
	auto.p.value <- c(auto.p.value,auto.regC$Cochrane.Orcutt$coefficients[2,4])
	rho <- c(rho,auto.regC$rho)
	regG <- lm(Gs~nuc.rescaled)
	abline(regG,col="#000000")
	intercept <- c(intercept,summary(regG)$coefficients[1,1])
	slope <- c(slope,summary(regG)$coefficients[2,1])
	r2 <- c(r2,summary(regG)$r.squared)
	adj.r2 <- c(adj.r2,summary(regG)$adj.r.squared)
	p.value <- c(p.value,summary(regG)$coefficients[2,4])
	auto.regG <- cochrane.orcutt(regG)
	auto.p.value <- c(auto.p.value,auto.regG$Cochrane.Orcutt$coefficients[2,4])
	rho <- c(rho,auto.regG$rho)
	regT <- lm(Ts~nuc.rescaled)
	abline(regT,col="#FF0000")
	intercept <- c(intercept,summary(regT)$coefficients[1,1])
	slope <- c(slope,summary(regT)$coefficients[2,1])
	r2 <- c(r2,summary(regT)$r.squared)
	adj.r2 <- c(adj.r2,summary(regT)$adj.r.squared)
	p.value <- c(p.value,summary(regT)$coefficients[2,4])
	auto.regT <- cochrane.orcutt(regT)
	auto.p.value <- c(auto.p.value,auto.regT$Cochrane.Orcutt$coefficients[2,4])
	rho <- c(rho,auto.regT$rho)
	reg.results <- data.frame(intercept=intercept,slope=slope,R2=r2,aR2=adj.r2,p=p.value,auto.p=auto.p.value,rho=rho,row.names=c("%A","%C","%G","%T"))
	
	#Plotting A-T skews.
	plot(nuc.rescaled,ATskews,type="p",main=paste(genome[g]," A-T skew",sep=""),sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),xlim=c(1,length(nuc)),ylim=c(-1,1),xlab="mtDNA",ylab="A-T skew",pch=21,col="#7F7F00",bg="#7F7F00")	

	#Removing NaN.
	nuc.rescaled <- nuc.rescaled[!is.nan(As)]
	As <- As[!is.nan(As)]
	Cs <- Cs[!is.nan(Cs)]
	Gs <- Gs[!is.nan(Gs)]
	Ts <- Ts[!is.nan(Ts)]

	#Testing autocorrelation using autocorrelograms.
	dev.new()
	par(mfrow=c(2,2))
	A.acf <- acf(As,plot=FALSE)
	C.acf <- acf(Cs,plot=FALSE)
	G.acf <- acf(Gs,plot=FALSE)
	T.acf <- acf(Ts,plot=FALSE)
	plot(A.acf,type="p",ylab="acf",ylim=c(-1,1),main="Adenosine",sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),ci.col="#00FF00",ci.type="ma",pch=21,col="#00FF00",bg="#00FF00")
	plot(C.acf,type="p",ylab="acf",ylim=c(-1,1),main="Cytidine",sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),ci.col="#0000FF",ci.type="ma",pch=21,col="#0000FF",bg="#0000FF")
	plot(G.acf,type="p",ylab="acf",ylim=c(-1,1),main="Guanosine",sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),ci.col="#000000",ci.type="ma",pch=21,col="#000000",bg="#000000")
	plot(T.acf,type="p",ylab="acf",ylim=c(-1,1),main="Thymidine",sub=paste("Starting gene: ",gene[g],"; wsize: ",wsize,"; wstep: ",wstep,".",sep=""),ci.col="#FF0000",ci.type="ma",pch=21,col="#FF0000",bg="#FF0000")
	mtext(paste("Autocorrelograms for ",genome[g],sep=""),line=-1,outer=TRUE)
	dev.copy2pdf(file=paste("./",genome[g],"_acf.pdf",sep=""))
	dev.off()

	#Writing results to output.
	if (g == 1) {
		cat(paste(genome[g],"\n",sep="\n"),file=paste("../",outfile,sep=""),sep="",append=FALSE)
		}
	else {
		cat("\n",file=paste("../",outfile,sep=""),sep="",append=TRUE)
		cat(paste(genome[g],"\n",sep="\n"),file=paste("../",outfile,sep=""),sep="",append=TRUE)
		}
	current.line <- ""
	for (i in 1:length(colnames(reg.results))) current.line <- paste(current.line,colnames(reg.results)[i],sep="\t")
	cat(paste(current.line,"\n",sep=""),file=paste("../",outfile,sep=""),sep="",append=TRUE)
	for (i in 1:length(rownames(reg.results))) {
		current.line <- rownames(reg.results)[i]
		for (j in 1:length(colnames(reg.results))) current.line <- paste(current.line,reg.results[i,j],sep="\t")
		cat(paste(current.line,"\n",sep=""),file=paste("../",outfile,sep=""),sep="",append=TRUE)
		}
	}
dev.copy2pdf(file=paste("../",plotfile,sep=""))
