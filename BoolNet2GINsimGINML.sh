#!/bin/sh
##########################################################################
# Author: Pedro T. Monteiro
# WebPage: http://pedromonteiro.org
#
# Please contribute for this script at:
# https://github.com/ptgm/BoolNetR2GINsim.git
#
# This script generates a GINsim file (GINML) containing a random model 
# with <nGenes> genes and a number of regulators per gene.
# The model interactions and parameters are computed using BoolNet R tool.
# The generated model is written in <filename>.ginml 
# This script does NOT consider any of the GINML graphical features.
##########################################################################

function throwError {
  echo "ERROR Message: $1"
	echo "USAGE: ./randomBoolNet.sh <filename> <nGenes> <reg-type:same|poisson> <arg> <param-type:random|andor>";
	echo " filename   - is the name of the file to be written without extension"
	echo " nGenes     - is the number of nodes in the network"
	echo " reg-type   - the type of regulator selection:"
	echo "    same    - all genes will have the same # regulators"
	echo "    poisson - the # regulators of each gene will follow a poisson distribution"
	echo " arg        - if (type=same) arg <- #reg else arg <- lambdaParam"
	echo " param-type - the type of parameter selection for each node"
	echo "    random  - all parameters are randomly selected (uniform distribution)"
	echo "    andor   - parameters are selected to perform an AND or a OR on regulators"
	exit;
}
if [ $# -ne 5 ]; then
	throwError "You must pass exactly 5 arguments (passed $#)"
fi

BASENAME=$1
NGEN=$2
REGTYPE=$3
ARG=$4
PARAMTYPE=$5

if [ $NGEN -lt 1 ]; then
	throwError "<nGenes> must be greater than 0";
fi
if [ $REGTYPE != "same" -a $REGTYPE != "poisson" ]; then
	throwError "<reg-type> must be one of the following types: same | poisson";
fi
if [ $REGTYPE == "same" ]; then
	if [ $ARG -lt 1 ]; then
	  throwError "<arg> must be greater than 0";
  fi
	if [ $NGEN -lt $ARG ]; then
		throwError "<nGenes> cannot be smaller than <arg>";
	fi
fi
if [ $REGTYPE == "poisson" ]; then
	if [ $ARG -lt 0 ]; then
		throwError "<arg> must be non-negative for a poisson distribution";
	fi
fi
if [ $PARAMTYPE != "random" -a $PARAMTYPE != "andor" ]; then
	throwError "<param-type> must be one of the following types: random | andor";
fi

GINML=$BASENAME.ginml
echo "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" > $GINML
echo "<!DOCTYPE gxl SYSTEM \"http://gin.univ-mrs.fr/GINsim/GINML_2_1.dtd\">" >> $GINML
echo "<gxl xmlns:xlink=\"http://www.w3.org/1999/xlink\">" >> $GINML
NODES=""
for (( i=1; i<=$NGEN; i++ )); do
	NODES="$NODES G$i";
done
echo "  <graph id=\"default_name\" class=\"regulatory\" nodeorder=\"$NODES\">" >> $GINML

/usr/bin/R --no-save $GINML $NGEN $REGTYPE $ARG $PARAMTYPE <<EOF
filename <- commandArgs()[3]
ngenes <- as.numeric(commandArgs()[4])
regtype <- commandArgs()[5]
arg <- as.numeric(commandArgs()[6])
paramtype <- commandArgs()[7]
if (regtype == "same") {
	nregs <- rep(arg, ngenes)
} else {
	# All genes have # regs = poisson(lambda) + 1
  nregs <- rpois(ngenes, lambda=arg) + 1
	for (i in 1:length(nregs)) {
		if (nregs[i] > ngenes) {
			nregs[i] <- ngenes;
		}
	}
}

vect2pos <- function(vect) {
	sz <- length(vect)
	n <- 1
	for (i in c(1:sz)) {
		if (vect[i]) {
			n <- n + 2^(sz-i);
		}
	}
	return(n)
}

library(BoolNet)
net <- generateRandomNKNetwork(n=ngenes, k=nregs, noIrrelevantGenes=FALSE)
if (paramtype != "random") {
	for (v in c(1:ngenes)) {
		net\$interactions[[v]]\$isAnd <- sample(c(TRUE,FALSE),1);
		for (reg in c(1:nregs[v])) {
			net\$interactions[[v]]\$posSign[reg] <- sample(c(TRUE,FALSE),1);
		}
		if (sample(c(TRUE,FALSE),1)) { # AND
			pos <- vect2pos(net\$interactions[[v]]\$posSign);
			valAll <- 0
			valPos <- 1
		} else {
			pos <- vect2pos(!net\$interactions[[v]]\$posSign);
			valAll <- 1
			valPos <- 0
		}
		for (i in 1:c(length(net\$interactions[[v]]\$func))) {
			net\$interactions[[v]]\$func[i] = valAll;
		}
		net\$interactions[[v]]\$func[pos] = valPos;
	}
}

as.binary <- function(v,maxvalue,base=2) {
  if (maxvalue==0)
		maxvalue <- 1
	ndigits <- 1 + floor(log(max(maxvalue), base))
	r <- rep(0,ndigits)
	for (i in 1:ndigits) {
		pow <- 2^(ndigits-i)
		r[i] <- (v %/% pow) %% 2
	}
	return(r)
}

# writing nodes
for (tg in c(1:ngenes)) {
	node <- paste("    <node id=\"G",tg,"\" maxvalue=\"1\">", sep="")
	func <- net\$interactions[[tg]]\$func %% 2
  for (parampos in c(1:length(func))) {
		if (func[parampos]) {
			node <- paste(node,"\n      <parameter", sep="")
			inters <- ""
			dec2bin <- as.binary(parampos-1,length(func)-1)
			regs <- net\$interactions[[tg]]\$input
			for (ireg in 1:length(regs)) {
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
for (tg in c(1:ngenes)) {
	regs <- net\$interactions[[tg]]\$input
	for (ireg in 1:length(regs)) {
		if (paramtype != "random") {
			if (net\$interactions[[tg]]\$posSign[ireg]) {
				sign <- "positive";
			} else {
				sign <- "negative";
			}
		} else {
			sign <- "unknown";
		}
		edge <- paste("    <edge id=\"G",regs[ireg],":G",tg,"\" from=\"G",regs[ireg],"\" to=\"G",tg,"\" minvalue=\"1\" sign=\"", sign, "\">\n    </edge>", sep="")
		write(edge, file=filename, append=TRUE)
  }
}
EOF

echo "  </graph>" >> $GINML
echo "</gxl>" >> $GINML

