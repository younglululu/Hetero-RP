# Hetero-RP
Heterogeneity Rescaling Pursuit: Enhanced Clustering and Classification in Integrative Genomics
===================

Heterogeneity Rescaling Pursuit (Hetero-RP) is a general preprocessing framework for integrative genomic studies. Hetero-RP is formulated as seeking a feature scaling transformation with the aid of implicitly existing auxiliary knowledge. And Hetero-RP is a a scalable and tuning-free algorithm to weigh important features more highly than less important ones. The effectiveness of Hetero-RP is demonstrated in several applications, including metagenomic contig binning, RBP binding site prediction, and cancer subtype stratification. The applications cover both clustering and classification problems. In each case, a combination of Hetero-RP with state-of-the-art approaches improves the yet satisfactory performance further, showing the wide applicability of the framework.


----------
Description
---------------
This package contains the following files and directories.

	code/buildKnnGraph.m -> To build the intrinsic structure of input data when the attraction-link set is sparse, 
	code/heteroRP.m -> calculate the rescaling, the key algorithm
	code/run_metagenome.m -> a demo on using Hetero-RP in metagenome data
	code/run_rbp.m ->  a demo on using Hetero-RP in rbp data
	code/cvx -> Software for Disciplined Convex Programming [1]
	data/ -> datasets directory
	
Please try to execute 'run_metagenome.m' or 'run_rbp.m' to learn how to use Hetero-RP to enhance the data from multiple sources in either clustering or classification scenarios. And please check 'heteroRP.m' for the detailed usage of the algorithm.
	

----------
Setup cvx
---------------

Before using Hetero-RP, the users are expected to setup the enviroment of cvx.

First of all, we change the current working directory to the code folder and unzip the cvx.
```sh
$ cd code/
$ unzip cvx.zip
```
Next, we open the matlab and change the current working directory to the code folder
```matlab
>> cd cvx
>> cvx_setup
>> cvx_startup
>> cd ..
```
Finally, we can run the demo
```matlab
>> run_metagenome
>> run_rbp
```
	
----------
Contacts and bug reports
------------------------
Please send bug reports, comments, or questions to 

Yang Lu: [ylu465@usc.edu](mailto:ylu465@usc.edu)

Prof. Fengzhu Sun: [fsun@usc.edu](mailto:fsun@usc.edu)


----------
Copyright and License Information
---------------------------------
Copyright (C) 2016 University of Southern California, Yang Lu

Authors: Yang Lu

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.


----------
References
---------------------------------
[1] cvxr.com/cvx/



Last update: 30-Dec-2016