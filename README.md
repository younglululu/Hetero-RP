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


    **Usage**:  heteroRP.py [-h] [-skip_row_header] [-skip_column_header]
                   [-object_in_row] [-sep SEP]
                   [-data_files DATA_FILES [DATA_FILES ...]]
                   [-signed_graph_file SIGNED_GRAPH_FILE]
                   [-check_feature_validity] [-label_file LABEL_FILE]
                   [-p_value P_VALUE] [-use_knn] [-k K] [-gammaVal GAMMAVAL]
                   [-lambdaVal LAMBDAVAL] [-output_dir OUTPUT_DIR] [-viz]


Input arguments:
  
- ``-skip_row_header`` is an optional argument to control whether the first column indicates the header. By default the program assumes there is no header explicitly given in the first column.Explicitly set ``-skip_row_header`` to be able to skip the first column as header for downstream analysis.

- ``-skip_column_header`` is an optional argument to control whether the first row indicates the header. By default the program assumes there is no header explicitly given in the first row. Explicitly set ``-skip_column_header`` to be able to skip the first row as header for downstream analysis.

- ``-object_in_row`` is an optional argument to control whether each row represents an object and each column represents a feature. By default the program assumes each row represents a feature and each column represents an object. Explicitly set ``-object_in_row`` to let each row represent an object and each column represent a feature.

- ``-sep SEP``, by default ``SEP`` is set as ``,``, specifies the delimiter to parse the data files.

- ``-data_files DATA_FILES [DATA_FILES ...]`` specifies the input data files (at least one). Each data file contains a table where each row correspond to a feature, and each column corresponds to an object. The row number can vary but the column number must be consistent.All values are separated by ``SEP``.

- ``-signed_graph_file SIGNED_GRAPH_FILE`` specifies the signed graph encoding either the positive-links and negative-links information. Each row represents one edge in the format: ``object_A SEP object_B SEP weight``. The edge is undirected. The index of objects starts from ``1`` instead of ``0``. The weight is ``1`` for positive-links and ``-1`` for negative links, respectively.
                        

- ``help``, by default ``True``, specifies whether the parser should
  automatically print the help message (supplied as ``doc``) and
  terminate, in case ``-h`` or ``--help`` option is encountered
  (options should exist in usage pattern, more on that below). If you
  want to handle ``-h`` or ``--help`` options manually (as other
  options), set ``help=False``.	


	
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
