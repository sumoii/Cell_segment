# Cell_segment
========================

Install Cell_segment
------------------------------------
```
$git clone https://github.com/
```
Dependencies
------------------------------------
1.python3.6+  
2.gcc 4.4+
3.numpy,pandas,anndata,tangram,densecrf

Usage
-----------------------------------
```usage: cell_segment.py [-h] -g IN_GEM -s SC_FILE -b BIN_SIZE -d DENSITY -w1
                       WEIGHT1 -w2 WEIGHT2 -o OUT_FILE -p PREFIX -n T_TIME
                       -f FEATURE

Cell segment (V1.0)

optional arguments:
  -h, --help   show this help message and exit
  -g IN_GEM    the gem of stero-seq
  -s SC_FILE   the h5ad of single cell
  -b BIN_SIZE  the initialize size of a bin (spot number = n*n)
  -d DENSITY   the gene density limit of a bin
  -w1 WEIGHT1  the weight of first Gaussian kernel
  -w2 WEIGHT2  the weight of second Gaussian kernel
  -o OUT_FILE  the initialize output file prefix
  -p PREFIX    the prefix of log file
  -n T_TIME    the times of running tangram
  -f FEATURE   the feature input of tangram
```
