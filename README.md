# Hetero-RP
Heterogeneity Rescaling Pursuit: Enhanced Clustering/Classification in Integrative Genomics
===================

High-throughput technologies have led to large collections of different types of biological data that provide unprecedented opportunities to unravel molecular heterogeneity of biological processes. Nevertheless, how to jointly explore data from multiple sources into a holistic, biologically meaningful interpretation remains challenging. 

Heterogeneity Rescaling Pursuit (Hetero-RP) is a scalable and tuning-free preprocessing framework, which weighs important features more highly than less important ones in accord with implicitly existing auxiliary knowledge. 

Hetero-RP is implemented in python, with dependencies on following packages:

	numpy  
	scipy
	pandas
	scikit-learn
	cvxopt
	diptest

----------
Setup Hetero-RP by Anaconda
---------------
Here we describe using Anaconda to run Hetero-RP on Linux/Unix. Anaconda is a tool to isolate your python installation, which allows you to have multiple parallel installations using different versions of different packages, and gives you a very convenient and fast way to install the most common scientific python packages. Anaconda can be downloaded from [here](https://www.continuum.io/downloads) 

After installing Anaconda, create a new environment that will contain the Hetero-RP installation:

```sh
$ conda create -n heterorp_env python=2.7.6
```

After creating the Anaconda environment, run the following command to activate it:

```sh
$ source activate heterorp_env
```

After that, install the Hetero-RP dependencies into this environment:

```sh
$ conda install -c anaconda numpy scipy pandas scikit-learn cvxopt
```

Finally, install diptest from the local package:

```sh
$ gunzip diptest-master.zip
$ cd diptest-master/
$ python setup.py install
$ cd ..
```

----------
Usage of Hetero-RP
------------------------
Hetero-RP uses several command line options, here is a
complete documentation of these. These can also be viewed by typing ``heteroRP.py
-h`` on the command line:

.. program-output:: (echo 'import conf'; cat ./heteroRP.py; echo 'args=arguments()') | python - --help
   :shell:




**Command:  ** heteroRP.py [-h] [-skip_row_header] [-skip_column_header]
                   [-object_in_row] [-sep SEP]
                   [-data_files DATA_FILES [DATA_FILES ...]]
                   [-signed_graph_file SIGNED_GRAPH_FILE]
                   [-check_feature_validity] [-label_file LABEL_FILE]
                   [-p_value P_VALUE] [-use_knn] [-k K] [-gammaVal GAMMAVAL]
                   [-lambdaVal LAMBDAVAL] [-output_dir OUTPUT_DIR] [-viz]

 - Main arguments:

	-D < dist >: Comma-separated list of distance measurements,  **E.g.** -D D2star,Ma,CVtree. The options include: 
				
	-F < fa_Dir >: Folder containing only fasta files with extension '.fasta', '.fa', and '.fna'.
	
	-I < fa_files >: Comma-separated list of sequence fasta files, e.g. -I speciesA.fa,speciesB.fa,speciesC.fa. Pairwise similarity is calculated based upon the sequences specified with this option.
	
	-K < intK >: Kmer Length.
	
----------
Contacts and bug reports
------------------------
Please send bug reports, comments, or questions to 

Yang Lu: [ylu465@usc.edu](mailto:ylu465@usc.edu)

Prof. Fengzhu Sun: [fsun@usc.edu](mailto:fsun@usc.edu)


----------
Copyright and License Information
---------------------------------
Copyright (C) 2017 University of Southern California, Yang Lu

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



Last update: 24-Jul-2017
