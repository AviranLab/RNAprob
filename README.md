# RNAprob

A tool for data-directed RNA secondary structure prediction. It uses structure probing data to restrain structure prediction, and outputs MFE structure and suboptimal structures. It is based on [RNAstructure] package.

## Citation

Fei Deng, Mirko Ledda, Sana Vaziri and Sharon Aviran. (2016). Data-directed RNA secondary structure prediction using probabilistic modeling. RNA, 22:1109–1119.



## Installation

First of all, enter the directory of RNAprob and type:
```sh
$ make RNAprob
```
The executable will then be placed in the directory RNAprob/exe. Note that you may need to change compiler setting in compiler.h based on your operating system.
Then, specify the location of thermodynamic parameters and training parameters by adding the following line to .bashrc (create a .bashrc file if it does not exist on your computer):
```sh
export DATAPATH=[directory in which RNAprob resides]/data_tables/
```

## Usage

Four variants of RNAprob are available, the general command is 
```sh
$ RNAprob <seq file> <ct file> -sh <shape file> [options]
```
* \<seq file\> : the name of a sequence file containing input sequence.

* \<ct file\> : the name for output ct file.

* \<shape file\> : the name of an input SHAPE profile.

 Specifically, the commands to run each of the four variants of RNAprob are as follows:
 
* RNAprob-3 : RNAprob  \<seq file\>  \<ct file\>  \-sh  \<shape file\>

* RNAprob-2 : RNAprob  \<seq file\>  \<ct file\>  \-sh  \<shape file\>  \-2s

* RNAprob-3s : RNAprob  \<seq file\>  \<ct file\>  \-sh  \<shape file\>  \-smooth

* RNAprob-2s : RNAprob  \<seq file\>  \<ct file\>  \-sh  \<shape file\>  \-2s  \-smooth

### scorer function
A scorer function that measures prediction accuracy of a predicted structure is also included. This function extends the [scorer] function provided in [RNAstructure] by adding the computation of Matthews Correlation Coefficient (MCC). To compile it, enter the directory of RNAprob and type:
```sh
$ make scorer
```
The executable will then be placed in the directory RNAprob/exe. The general command to run the scorer function is as follows:
```sh
$ scorer <predicted ct> <accepted ct> <output file> [options]
```
For more details, please refer to http://rna.urmc.rochester.edu/Text/scorer.html.


## Licence
RNAprob is free for non-commercial research. For commercial use, please contact the authors.


[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [scorer]: http://rna.urmc.rochester.edu/Text/scorer.html
   [RNAstructure]: http://rna.urmc.rochester.edu/RNAstructure.html
