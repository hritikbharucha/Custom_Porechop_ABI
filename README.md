# Porechop_abi 

Porechop_abi (*ab initio*) is an extension of Porechop that is able to infer the adapter sequence from the Oxford Nanopore reads. It discovers the adapter sequence from the reads using approximate k-mers and assembly, and add the sequence found to the adapter list (adapters.py file). It then runs Porechop and remove adapters with all usual options : adapters on the ends of reads are trimmed off, and when a read has an adapter in its middle, it is treated as chimeric and chopped into separate reads. 

Porechop_abi is not specifically designed for demultiplexing and barcoded sequences. In this case, it has the same behaviour as the original Porechop.

## Requirements 

The requirement are the same as Porechop (Oct 2018 version), except you will also need to install the graph library [networkx](https://networkx.github.io/).

### Porechop requirements

See Porechop [documentation](README_PORECHOP.md) 

### Networkx

~~~
pip install networkx
~~~

## Installation 


Installation is similar to Porechop installation.

First, clone the repository using the recursive option.
Then, just install as described in the Porechop [documentation](README_PORECHOP.md).  

Running the `setup.py` script will compile the C++ components of Porechop_ABI and install a `porechop` executable:

```bash
git clone --recursive https://github.com/qbonenfant/Porechop_ABI.git
cd Porechop_ABI
python3 setup.py install
porechop -h
```

Notes:
* If the last command complains about permissions, you may need to run it with `sudo`.
* Install just for your user: `python3 setup.py install --user`
    * If you get a strange "can't combine user with prefix" error, read [this](http://stackoverflow.com/questions/4495120).
* Install to a specific location: `python3 setup.py install --prefix=$HOME/.local`
* Install with pip (local copy): `pip3 install path/to/Porechop_ABI`
* If you'd like to specify which compiler to use, set the `CXX` variable: `export CXX=g++-6; python3 setup.py install`
* Porechop includes `ez_setup.py` for users who don't have [setuptools](https://pypi.python.org/pypi/setuptools) installed, though that script is [deprecated](https://github.com/pypa/setuptools/issues/581). So if you run into any installation problems, make sure setuptools is installed on your computer: `pip3 install setuptools
* If you had a previous installation of Porechop, it will try to replace it. Setting a different installation folder while already having Porechop installed may lead to conflict when calling `porechop` from command line.


### Build and run without installation

By simply running `make` in Porechop's directory, you can compile the C++ components but not install an executable. The program can then be executed by directly calling the `porechop-runner.py` script.

```bash
git clone --recursive https://github.com/qbonenfant/Porechop_ABI.git
cd Porechop_ABI
make
./porechop-runner.py -h
```


## Usage

Porechop_abi offers two new options.

```bash
Porechop --ab_initio
```

This flag allows to first guess the adapter sequence from the reads, add the sequence to the list of Porechop adapters and then run Porechop as usual.  It is compatible with all Porechop options, but current version may behave poorly on on barcoded reads.


```bash
Porechop --guess_adapter_only
```

This flag allows to only guess the adapter sequence from the reads. It then stops  the execution of the program, without 
trimming the reads. 

For all other usages and description of the output files, you can refer to the Porechop [documentation](README_PORECHOP.md). 

## Config file

It is also possible to tune the parameters used in the algorithm used to reconstruct the adapter sequence. Those advanced options are accessible in the config file located in the Porechop_ABI/porechop/ab_initio.config. Note that default values work just fine in practice, and it is most likely that you will not need to edit this file. 

If you installed Porechop using the setup.py script, take note that the actual config file used when using porechop from command line will be the one stored in your installation folder (default: /usr/local/lib/python3.__X__/dist-packages/porechop/porechop/ab_initio.config).


## Contributors

The original Porechop program was provided by Ryan Wick.

The ab initio extension is developed by Quentin Bonenfant, Laurent Noé and Hélène Touzet.

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
