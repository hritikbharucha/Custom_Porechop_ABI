# Porechop_ABI

Porechop_abi (*ab initio*) is an extension of Porechop that is able to infer the adapter sequence from the Oxford Nanopore reads. Adapters sequences are discovered directly from the reads using approximate k-mers counting and assembly. Inferred sequences can either be automatically added to the adapter database for the current run (adapters.py file), or just displayed.

Trimming can then occur as usual, using all standard Porechop options.

Porechop_ABI is not designed to infer barcoded sequences adapters, but will report several sequences if a mix of adapters is used.
Demultiplexing should be done using standard Porechop commands or more appropriate tools.  

# Table of Content
* [Requirements](#requirements)
   * [Porechop Reqirements](#porechop-requirements)
   * [Python version](#python-version)
   * [Networkx](#networkx)
* [Installation](#installation)
   * [Install from source](#install-from-source)
   * [Build and run without installation](#build-and-run-without-installation)
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

The requirement are almost the same as Porechop (Oct 2018 version). You will need an updated version of Python (>= 3.6) and also need to install the graph library [networkx](https://networkx.github.io/).

### Porechop requirements

See Porechop [documentation](README_PORECHOP.md)

### Python Version
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


Installation is similar to Porechop installation.
First, clone the repository using the recursive option.
If you are familiar with standard Porechop, you can just install as described in the Porechop [documentation](README_PORECHOP.md).

### Install from source
Running the `setup.py` script will compile the C++ components of Porechop_ABI and install a `porechop` executable:

```bash
git clone --recursive https://github.com/bonsai-team/Porechop_ABI.git
cd Porechop_ABI
python3 setup.py install
porechop -h
```

Notes:
* If the install command (3rd line) complains about permissions, you may need to run it with `sudo`.
* In case you get "could not find cpp_function.so" or similar error, you may need to add execution permission to compiled objects.
`sudo chmod +x <your_install folder>/cpp_function.so`
* To install just for your user: `python3 setup.py install --user`
    * If you get a strange "can't combine user with prefix" error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Porechop_ABI`
* If you'd like to specify which compiler to use, set the `CXX` variable: `export CXX=g++-6; python3 setup.py install`
* Porechop includes `ez_setup.py` for users who don't have [setuptools](https://pypi.python.org/pypi/setuptools) installed, though that script is [deprecated](https://github.com/pypa/setuptools/issues/581). So if you run into any installation problems, make sure setuptools is installed on your computer: `pip3 install setuptools
* If you had a previous installation of Porechop, it will try to replace it. Setting a different installation folder while already having Porechop installed may lead to conflict when calling `porechop` from command line.
* For retro-compatibility issues, the name of the program remains `porechop`. Since no changes has been made to the internal functions of Porechop, transition to Porechop_ABI should be transparent.

### Build and run without installation

By simply running `make` in Porechop's directory, you can compile the C++ components but not install an executable. The program can then be executed by directly calling the `porechop-runner.py` script.

```bash
git clone --recursive https://github.com/bonsai-team/Porechop_ABI.git
cd Porechop_ABI
make release
./porechop-runner.py -h
```

## Quick usage examples
__Basic adapter inference and trimming:__<br>
`porechop -abi -i input_reads.fastq -o output_reads.fastq`

__Only display inferred adapters:__<br>
`porechop -abi --guess_adapter_only -i input_reads.fastq.gz -o output_reads.fastq`
`porechop -abi -go -i input_reads.fastq -v 0 -o output_reads.fastq.gz`


__Building a stronger consensus using more core module runs:__<br>
`porechop -abi -nr 20 -cr 30 -i input_reads.fastq.gz -o output_reads.fastq`

### Test Porechop_ABI
Two additionnal test files are provided to test Porechop_ABI<br>
__Simulated data (with less core module runs than usual, discarding default database)__<\br>
`porechop -abi -go -dd -nr 5 -cr 15 -i test/test_simulated_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`
__Real data (standard parameters, discarding default database)__<\br>
`porechop -abi -go -dd -i test/test_realdata_10k_read.fasta -tmp /tmp/pabi_temp -o /dev/null`

## Usage

Porechop_abi offers several new options.

### Simple runs
```bash
--ab_initio / -abi
```

This flag allows to first guess the adapters sequences from the reads, add the sequence to the list of Porechop adapters and then run Porechop as usual. It is compatible with all Porechop options, but behave poorly on on barcoded reads.


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
Each count file is exported separately, and can be reviewed if needed. (default: 10, set to 1 for single run mode)

```bash
--consensus_run / -cr [int]
```
If using multi-run option above 1, set the number of additional run performed if no stable consensus is immediatly found.
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
Only select the best x consensus sequences from all consensus found. (default: 0)



## Config file

It is also possible to tune the parameters of the core module used to reconstruct the adapter sequence. Those advanced options are accessible in the config file located at /Porechop_ABI/porechop/ab_initio.config. Note that default values work just fine in practice, and it is most likely that you will not need to edit this file.

If you installed Porechop using the setup.py script, take note that the actual config file used when using porechop from command line will be the one stored in your installation folder (example: /usr/local/lib/python3.__X__/dist-packages/porechop/porechop/ab_initio.config).

In the case you still want to use custom options, it is recommanded to pass a config file as a parameter to Porechop_ABI using -abc / --ab_initio_config instead of editing the default config.


```bash
Porechop --ab_initio_config YOUR_CONFIG.config ...
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
  nb_of_runs  : Number of run to perform and export at once.
```

**Notes:**
Verbosity, number of threads, and number of runs are copied from the argmuments supplied to Porechop_ABI.
Changing the value of these parameters in the config file will have no effect.
If you need to adjust one of them, refer to porechop -h


## Contributors

The original Porechop program was provided by Ryan Wick.

The ab initio extension is developed by Quentin Bonenfant, Laurent Noé and Hélène Touzet.


## References
If you want more details about the algorithm, internal working and performances of this tool, a preprint is available here:<\br>
https://www.biorxiv.org/content/10.1101/2022.07.07.499093

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
