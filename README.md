# Porechop_ABI

Porechop_abi (*ab initio*) is an extension of Porechop that is able to infer the adapter sequence from the Oxford Nanopore reads. Adapters sequences are discovered directly from the reads using approximate k-mers counting and assembly. Inferred sequences can either be automatically added to the adapter database for the current run (adapters.py file), or just displayed.

Trimming can then occur as usual, using all standard Porechop options.

Porechop_ABI is not designed to infer barcoded sequences adapters, but will report several sequences if a mix of adapters is used.
Demultiplexing should be done using standard Porechop commands or more appropriate tools.  

# Table of Content
* [Requirements](#requirements)
   * [Operating System version](#operating-system-version)
   * [SeqAn 2.4 and zlib 1.2](#seqan-24-and-zlib-12)
   * [Compiler and C++17](compiler-and-c17)
   * [Python version](#python-version)
   * [Networkx](#networkx)
* [Installation](#installation)
* [Quick usage examples](#quick-usage-examples)
   * [Test Porechop_ABI](#test-porechop_abi)
* [Usage](#usage)
   * [Simple runs](#simple-runs)
   * [Useful options](#useful-options)
   * [Advanced settings](#advanced-settings)
   * [Consensus options](#consensus-options)
* [Config File](#config-file)
   * [List of possible parameters](#list-of-possible-parameters)
* [Contributors](#contributors)
* [References](#references)
* [License](#license)


## Requirements 

If you want to install Porechop_ABI, your system must satisfy some requirements.<br>
In the case you chose to install using Conda, you can ignore this section (except the [OS](operating_system_version) requirement) and just use an up-to-date version of conda to install Porechop_ABI.
Check the [installation](INSTALL.md) guide for more details.


### Operating System version
* Mac OS >= 10.12
* Linux, any version with compatible python and compiler requirements should work.


### SeqAn 2.4 and zlib 1.2
Seqan 2.4 and zlib 1.2 are both required library for Porechop_ABI to compile properly.
Both are embeded in this repository, and should be automatically linked during compilation.<br>

* Porechop_ABI updated the already embeded SeqAn library used by Porechop to version 2.4 (previously 2.3). This led to some changes in the requirements.
* To match existing Porechop ability to open compressed files, zlib 1.2 was added to the `include` folder.


### Compiler and C++17
Because of ssome specific features integrated in SeqAn 2.4 and modifications required to keep mac OS support, your compiler must be able to work with C++17.

* If you're using [GCC](https://gcc.gnu.org/), version 8 or later is required (check with `g++ --version`).
* Recent versions of [Clang](http://clang.llvm.org/)(>=5.0.0) and [ICC](https://software.intel.com/en-us/c-compilers) should also work

If for some reason you are interested in a C++14 compatible release, please open an issue here. Such install *should be* possible with additionnal work.


### Python Version
Porechop is written mainly in Python3.
Python version must be 3.6 or above.<br>

With [X] being the version of python you want to install:

On Linux
~~~
sudo apt-get install python3.[X]
~~~

On Mac os (under macports):
~~~
sudo port install python3[X]
~~~

### Networkx

With Pip:
~~~
pip install networkx
~~~
If you have permission errors, you may want to install it for the current user only:
~~~
pip install --user networkx 
~~~

On linux:
~~~
sudo apt-get install python-networkx
~~~
On Mac os (under macports):
~~~
sudo port install py3[X]-networkx 
~~~
With [X] being the version of python you use.


## Installation 

Porechop_ABI can be installed by several means. 
The easiest way to install Porechop_ABI is to use the conda package management system and install it from the bioconda channel in a clean environment.
All details about installation procedures are listed in the [installation](INSTALL.md) guide.<br>


## Quick usage examples
__Basic adapter inference and trimming:__<br>
`porechop_abi -abi -i input_reads.fastq -o output_reads.fastq`

__Only display inferred adapters:__<br>
`porechop_abi -abi --guess_adapter_only -i input_reads.fastq.gz -o output_reads.fastq`
`porechop_abi -abi -go -i input_reads.fastq -v 0 -o output_reads.fastq.gz`


__Building a stronger consensus using more core module runs:__<br>
`porechop_abi -abi -nr 20 -cr 30 -i input_reads.fastq.gz -o output_reads.fastq`

### Test Porechop_ABI
Two additionnal test files are provided to test Porechop_ABI<br>

Faster test: __Simulated data (with less core module runs than usual, discarding default database)__<br>
`porechop_abi -abi -go -dd -nr 5 -cr 15 -i test/test_simulated_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`<br>
Slower test: __Real data (standard parameters, discarding default database)__<br>
`porechop_abi -abi -go -dd -i test/test_realdata_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`<br>


[//]: # (TODO: Add expected results for P_ABI tests.)


## Usage

Porechop_abi offers several new options.

### Simple runs
```bash
--ab_initio / -abi
```

This flag allows to first guess the adapters sequences from the reads, add the sequence to the list of Porechop adapters and then run Porechop as usual. It is compatible with all Porechop options, but behave poorly on on barcoded reads. Adapter sequence inferrence is not activated by default.


```bash
--guess_adapter_only / -go
```

This flag allows to only guess and print the adapter sequence from the reads. It then stops the execution of the program, without 
trimming the reads. 

For all other usages and description of the output files, you can refer to the Porechop [documentation](README_PORECHOP.md). 

### Useful options 
We also added several options to tune both the main porechop implementation and the ab-initio step to your needs.
  
```bash
-abc / --ab_initio_config [path_to_file]
```
Allow you to set a custom config file for the ab_initio phase (default file in Porechop folder)
The config file come with it's own set of parameters check the "Config File" section for more details.

```bash
-tmp / --temp_dir [path_to_folder]
```
Path to a writable temporary directory, used to store count file and temporary fasta files.
The provided directory will be created if it does not exists, if the path is writable.
By default, Porechop_ABI creates a ./tmp/ folder in the working directory (where the command is launched).


```bash
--custom_adapters / -cap [path_to_file]
```     

You may now set custom adapters without editing the adapter.py file.
These adapters are only used at runtime, not added to the adapter.py file.
Custom adapters must be stored in a text file following this format:
```    
    line 1: Adapter name
    line 2: Start adapter sequence
    line 3: End adapter sequence
    --- repeat for each adapter pair---
```

If your adapters does not contains the start or end sequence,
just put an empty line

Each adapters will be named using the adapter name as a prefix, like:

```
"adapter_name_Top"
"adapter_name_Bottom"
```


```bash
--discard_database  / -ddb
```                              
In case you supply your own adapters, it can be usefull (and faster) to ignore the adapters from the Porechop database.
This option was added for this situation, and require either ab-initio (-abi) or a custom adapter (-cap) to be set.(default: False)


### Advanced settings

```bash
--export_graph [output_file]
```    
Let you export the assembly graph with some annotations to the provided path.
Graphs are exported in graphml format, and can be visualised with programs such as Gephi.


*The following options are specific to the drop-cut algorithm used to asjust sequences length after the (greedy) assembly step.
Only experienced users shoud try to adjust these settings.*

```bash
--window_size / -ws  [int]
```   
Change the size of the smoothing window used in the drop cut algorithm. (set to 1 to disable smoothing).


```bash
--no_drop_cut / -ndc
```   
Disable the drop cut step entirely (default: False)


### Consensus options:
Porechop_ABI build consensus sequences from multiple samples to increase confidence in the inferred adapter sequence.
Some parameters can be adjusted from the command line.

```bash
--multi_run / -mr  [int]
```   
Number of putative adapters to reconstruct to ensure the inferred adapters sequences are stable.
Each count file is exported separately, and can be reviewed if needed. (default: 10, set to 1 for single sample mode)

```bash
--consensus_run / -cr [int]
```
If using multi-run option above 1, set the number of additional sample performed if no stable consensus is immediatly found.
All inferred adapter sequences are then clustered by similarity and consensus are generated for each cluster.
(default: 20)
  

```bash
--export_consensus / -ec [output_file]
```
Path to export the intermediate adapters found in consensus mode.


```bash
--all_above_x / -aax [int]
```
Only select consensus sequences if they are build using at least x percent of the total adapters.
Default is 10% (3 sequence over 30). You can set it to 0 to keep all consensus.
  
```bash
--best_of_x / -box [int]
```
Only select the best x consensus sequences from all consensus found. (default: 0, reports all consensus)



## Config file

It is also possible to tune the parameters of the core module used to reconstruct the adapter sequence. Those advanced options are accessible in the config file located at /Porechop_ABI/porechop/ab_initio.config. Note that default values work just fine in practice, and it is most likely that you will not need to edit this file.

If you installed Porechop using the setup.py script, take note that the actual config file used when using porechop from command line will be the one stored in your installation folder (example: /usr/local/lib/python3.__X__/dist-packages/porechop/porechop/ab_initio.config).

In the case you still want to use custom options, it is recommanded to pass a config file as a parameter to Porechop_ABI using -abc / --ab_initio_config instead of editing the default config.


```bash
porechop_abi --ab_initio_config YOUR_CONFIG.config ...
```
### List of possible parameters:
```
  // Sampling options
  sn        : Number of read to sample (default: 40 000 reads)
  sl        : Size of the samples (default : 100 bases)
  limit     : Number of k-mer kept after counting. (default: 500)
  solid_km  : Use solid kmers instead of a hard limit. (specify the minimum count for a k-mer, will override "limit" option).
  
  // K-mer options
  k           : K-mer size (between 2 and 32, default 16).
  lc          : Low Complexity filter value (default 1.0).
  forbid_kmer : Forbidden k-mer list. (list of k-mers, one per line).
  exact_out   : Also export the exact count (mainly for debug / control)
  
  
  // Misc
  v           : Verbosity level (0: minimal output, 1: normal, 2: debug).
  nb_thread   : Number of threads used for the approximate counting.
  skip_end    : Do not process end adapters. (mostly for debug)
  nb_of_runs  : Number of time the sampling and count is performed. (Default: 10)
```

**Notes:**
Verbosity, number of threads, and number of runs are copied from the argmuments supplied to Porechop_ABI.
Changing the value of these parameters in the config file will have no effect.
If you need to adjust one of them, refer to porechop -h


## Contributors

The original Porechop program was provided by Ryan Wick.

The ab initio extension is developed by Quentin Bonenfant, Laurent Noé and Hélène Touzet.


## References
If you want more details about the algorithm, internal working and performances of this tool, a preprint is available here:<br>
https://www.biorxiv.org/content/10.1101/2022.07.07.499093

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
