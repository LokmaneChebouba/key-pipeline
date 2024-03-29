### Readme to be completed..
# key-pipeline
A computational model to identify key protein-complexes associated to tumor progression
### 1 Introduction
We introduce the ```key-pipeline```, a python package implementing a workflow for identifying key protein-complexes associated to tumor progression. Our software allows researchers to (i) find the upstream/downstream paths starting from a
couple of root nodes in a network using an embeded tool which is [pathrider](https://github.com/arnaudporet/pathrider), a tool developed in our team
to this purpose; (ii) then check the consistency of our data sets and provides explanations for inconsistencies using iggy tool; (iii) validate the predictions made by iggy by computing the number of predictions matching the related experimental fold-change from ICGC data; (iv) design a stability test by comparing prediction on subsets of observations with predictions using all observations; and finally (v) plot both precision scores for each sampling, and the evolution of the prediction compared to the entire set of observations.
### 2 Prerequisites
key-pipeline is a python application that uses many libraries and tools. The easiest way to obtain all depencies packages is using Anaconda. First install either Anaconda or Miniconda, download requirements.txt file and then run:

```
$ conda install −−file requirements.txt
```
The ```requirements.txt``` file contain all depencies required for the successful execution of the program.
###### example of requirements file
```
# This file may be used to create an environment using:
# $ conda create −−name <env> −− file < this file >
# platform: linux−64
plotly=3.10.0=py_0
python=3.7.3= h0371630_0
r=3.5.1=r351_0
```
In this case, the packages plotly, python and r will be installed.
### 3 Usage
Our tool provides a command line interface (CLI), it can be run by entering the arguments file (see below). By default, all the steps of the methods will be run, if you want to run just one or more steps, you can enter the number of steps you want by using ```--steps``` followed by the numbers of desired steps separated by ”,”. The ```--help``` provides the help message describing required inputs and available options. It implements all the steps in the workflow described before. Each step will output one or more files. In general, the output of one step corresponds to the input of another one. This enables a straightforward application of the workflow for users without programming expertise. The typical usage is:
```
$ python pipe.py @arguments.txt
```
Where ```arguments.txt``` contains all the arguments needed for the execution of the pipeline.
####### Example of arguments file
Here is an example of an arguments file:
```
--sif=/home/graph.sif
--dir=up
--icgc=/home/ICGC_data.csv
--b=home/black_listed_genes.txt
--start_sampling=10
--stop_sampling=15
--step_sampling=5
--numbers_run=2
```
For more options you can ask for help:
```
$ python pipe.py --help
```
```
usage: Usage : python pipe.py@arguments_file [--steps]
key-pipeline: require a file preceded by ’@’ and must contain all the required arguments cited below :
optional arguments :
-h , --help                         show this help message and exit
optional arguments :
--dir DIR                           follows the up stream ( ”up” ) or the down stream ( ”down” )
--steps STEPS                       specify the number of the steps to run :
                                    [ 1 ] Threshold extraction
                                    [ 2 ] Pathfinder
                                    [ 3 ] Iggy
                                    [ 4 ] Cross-Validation
                                    [ 5 ] Plots
                                    ( I f many ,separate by ’,’). Default run all steps
--start_sampling START_SAMPLING     The start sampling percentage. Default=10
--stop_sampling STOP_SAMPLING       The stop sampling percentage. Default=100
--step_sampling STEP_SAMPLING       The step of sampling. Default=5
--numbers_run NUMBERS_RUN           The number of runs on each step. Default=100
required arguments:
--sif SIF                           influence graph in SIF format
--icgc ICGC                         ICGC file
--b B                               a list of blacklisted genes (weakly expressed)
```
###### Assumptions 
We assume that the network provided in the arguments file must be in the SIF file format
###### Case studie ...
### 4 Inputs/Outputs
Each step of the tool need some inputs, and produce outputs:
