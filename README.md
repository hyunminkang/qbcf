# qbcf

Utilities for quantitative analysis for BCF/VCF files

## Overview

`qbcf` is a collection of C++ tools to facilitate quantitative analysis related to BCF/VCF files.

## Installing qbcf

Before installing `qbcf`, you need to install
[htslib](https://github.com/samtools/htslib) and
[qgenlib](https://github.com/hyunminkang/qgenlib) in the same directory you
want to install `qbcf' (i.e. `qbcf`, `htslib`, and `qgenlib` should be
siblings directories). You also need [cmake](https://cmake.org/) installed in your system.

After installing `htslib` and `qgenlib`, you can clone the current snapshot of this repository to install as well

<pre>
$ mkdir build

$ cd build

$ cmake ..
</pre>

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For qgenlib:
  - $ cmake -DQGEN_INCLUDE_DIRS=/qgenlib_absolute_path/include/  -DQGEN_LIBRARIES=/qgenlib_absolute_path/lib/qgenlib.a ..
  
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/htslib_absolute_path/include/  -DHTS_LIBRARIES=/htslib_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

Finally, to build the binary, run

<pre>
$ make
</pre>

### List of tools contained in `qbcf`

`qbcf` contains many in-house C++ tools that are currently under
the hood development phase. To list the available commands of tools, type:

<pre>
qbcf --help
</pre>

To see the usage of individual commands, type:

<pre>
qbcf [command] --help
</pre>
