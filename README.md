BoolNetR2GINsim
===============

The aim of this project is to provide a collection of scripts 
to generate random Boolean logical models using BoolNetR.


Requirements?
-------------

You should have R statistical software environment installed 
(http://www.r-project.org) together with the BoolNetR package
(http://cran.r-project.org/web/packages/BoolNet/index.html).


How to use it?
--------------

You can run it using the following command line format:

    sh randomBoolNet.sh <filename> <nGenes> <reg-type:same|poisson> <arg> <param-type:random|and|or>

where the arguments have the following semantics:
  * filename   - is the name of the file to be written without extension
  * nGenes     - is the number of nodes in the network
  * reg-type   - the type of regulator selection:
  *    same    - all genes will have the same # regulators
  *    poisson - the # regulators of each gene will follow a poisson distribution
  * arg        - if (type=same) arg <- #reg else arg <- lambdaParam
  * param-type - the type of parameter selection per node
  *    random  - randomly selected following an uniform distribution
  *    and     - parameters are selected to perform an AND on regulators
  *    or      - parameters are selected to perform an OR on regulators

Authors
-------
Pedro T. Monteiro
