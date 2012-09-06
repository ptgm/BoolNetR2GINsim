#!/bin/sh
##########################################################################
# Author: Pedro T. Monteiro
# WebPage: http://pedromonteiro.org
#
# This script generates a GINsim file (GINML) containing a random model 
# with <nGenes> genes and a number of <nRegs> regulators per gene.
# The model interactions and parameters are computed using BoolNet R tool.
# The generated model is written in <filename>.ginml 
# This script do NOT consider any of the GINML graphical features.
##########################################################################

if [ $# -ne 3 ]; then
	echo "USAGE: ./randomBoolNet.sh <filename> <nGenes> <nRegs>";
	echo " filename - is the name of the file to be written without extension"
	echo " nGenes   - is the number of nodes in the network"
	echo " nRegs    - is the number of regulators per node"
fi
BASE=$1
NGEN=$2
NREGS=$3
if [ $NGEN -lt 1 ]; then
	echo "ERROR: <nGenes> must be greater than 0";
	exit;
fi
if [ $NREGS -lt 1 ]; then
	echo "ERROR: <nRegs> must be greater than 0";
	exit;
fi
if [ $NGEN -lt $NREGS ]; then
	echo "ERROR: <nGenes> cannot be smaller than <nRegs>";
	exit;
fi

GINML=$BASE.ginml
echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > $GINML
echo "<!DOCTYPE gxl SYSTEM \"http://gin.univ-mrs.fr/GINsim/GINML_2_1.dtd\">" >> $GINML
echo "<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">" >> $GINML
NODES=""
for (( i=1; i<=$NGEN; i++ )); do
	NODES="$NODES G$i";
done
echo "  <graph id=\"default_name\" class=\"regulatory\" nodeorder=\"$NODES\">" >> $GINML

R --no-save $GINML $NGEN $NREGS <<EOF
filename <- commandArgs()[3]
ngenes <- as.numeric(commandArgs()[4])
nregs <- as.numeric(commandArgs()[5])

library(BoolNet)
net <- generateRandomNKNetwork(n=ngenes, k=nregs, noIrrelevantGenes=FALSE)

as.binary <- function(v,maxvalue,base=2) {
	ndigits <- 1 + floor(log(max(maxvalue), base))
	r <- rep(0,ndigits)
	for (i in 1:ndigits) {
		pow <- 2^(ndigits-i)
		r[i] <- (v %/% pow) %% 2
	}
	return(r)
}

# writing nodes
for (tg in c(1:length(net\$interactions))) {
	node <- paste("    <node id=\"G",tg,"\" maxvalue=\"1\">", sep="")
	func <- net\$interactions[tg][[1]]\$func %% 2 # [ TRUE ... FALSE ... FALSE ]
  for (parampos in c(1:length(func))) {
		if (func[parampos]) {
			node <- paste(node,"\n      <parameter", sep="")
			inters <- ""
			dec2bin <- as.binary(parampos-1,length(func)-1)
			regs <- net\$interactions[tg][[1]]\$input
			for (ireg in 1:c(length(regs))) {
				if (dec2bin[ireg] == 1)
					inters <- paste(inters," G",regs[ireg],":G",tg, sep="")
			}
			if (nchar(inters)>0) {
				inters <- paste(" idActiveInteractions=\"",inters,"\"", sep="")
			}
			node <- paste(node,inters," val=\"1\"/> ", sep="")
    }
  }
	node <- paste(node,"\n    </node>", sep="")
	write(node, file=filename, append=TRUE)
}

# writing edges
for (tg in c(1:length(net\$interactions))) {
  for (reg in net\$interactions[tg][[1]]\$input) {
    edge <- paste("    <edge id=\"G",reg,":G",tg,"\" from=\"G",reg,"\" to=\"G",tg,"\" ", sep="")
    edge <- paste(edge, "minvalue=\"1\" sign=\"unknown\">\n    </edge>", sep="")
    write(edge, file=filename, append=TRUE)
  }
}
EOF

echo "  </graph>" >> $GINML
echo "</gxl>" >> $GINML

