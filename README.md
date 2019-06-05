# Porechop_abi 

Porechop_abi (*ab initio*) is an extension of Porechop that is able to infer the adapter sequence from the Oxford Nanopore reads. It discovers the adapter sequence from the reads using approximate k-mers and assembly, and add the sequence found to the adapter list (adapters.py file). It then runs Porechop and remove adapters with all usual options : adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads. 

Porechop_abi is not specifically designed for demultiplexing and barcoded sequences. In this case, it has the same behaviour as the original Porechop.

## Requirements 

The requirement are the same as Porechop (Oct 2018 version), except you will also need to install the graph library [networkx] (https://networkx.github.io/).

### Porechop requirements

See Porechop [documentation](README_PORECHOP.md) 

### Networkx

~~~
pip install networkx
~~~

## Installation 

First, clone the repository using the recursive option :

```bash
git clone --recursive https://github.com/qbonenfant/Porechop_ABI.git
```

Then, just install as described in the Porechop [documentation](README_PORECHOP.md).  


## Usage

Porechop_abi offers two new options.

```bash
Porechop –ab_initio
```

This flag allows to first guess the adapter sequence from the reads, add the sequence to the list of Porechop adapters and then run Porechop as usual.  It is compatible with all Porechop options.

```bash
Porechop –guess_adapter_only
```

This flag allows to only guess the adapter sequence from the reads. It then stops  the execution of the program, without trimming the reads. 

It is also possible to tune the parameters used in the algorithm used to reconstruct the adapter sequence. Those advanced options are accessible in the config file XXX. Note that default values work just fine in practice, and it is most likely that you will not need to edit this file. 

For all other usages and description of the output files, you can refer to the Porechop [documentation](README_PORECHOP.md). 


## Contributors

The original Porechop program was provided by Ryan Wick.

The ab initio extension is developed by Quentin Bonenfant, Laurent Noé and Hélène Touzet.

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
