# [An R package `PIL`](https://github.com/zhangzhiqiang2020/PIL)

## @ Overview

> The `PIL` is an R package enabling multi-omics integration and network-driven drug target prioritisation for Lymphoma.

## @ Installation

### 1. Install R

Please install R (version 4.2.2 or above); see https://cran.r-project.org

If installed on `Ubuntu` (assuming you have a `ROOT (sudo)` privilege), please do so below

```ruby
sudo su
# here enter your password

wget http://www.stats.bris.ac.uk/R/src/base/R-4/R-4.2.2.tar.gz
tar xvfz R-4.2.2.tar.gz
cd ~/R-4.2.2
./configure
make
make check
make install
R # start R
```

### 2. Install R packages

```ruby
R # start R

# if the package 'BiocManager' not installed, please do so
if(!("BiocManager" %in% rownames(installed.packages()))) install.packages("BiocManager")

# first, install basic packages: remotes, tidyverse
BiocManager::install(c('remotes','tidyverse'), dependencies=T)

# then, install the package 'PIL' (now hosted at github)
BiocManager::install("zhangzhiqiang2020/PIL", dependencies=T, force=T)

# check the package 'PIL' successfully installed
library(help=PIL)
```


## @ Contact

> Please drop [email](zqzhang94@sjtu.edu.cn) for bug reports or other enquiries.
