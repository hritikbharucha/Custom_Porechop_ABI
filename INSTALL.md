# Installation

Porechop_ABI can be installed using several methods, that can be separated in two categories:
* Using package managers such as Conda or Pip.
* Installing from sources using scripts or make.

Except for conda installation, all requirements described in [README](README.md) must be satisfied before install.


# Table of Content
* [Install with Conda](#install-with-conda)
* [Fetching source code](#fetching-source-code)
* [Install with PIP](#install-with-pip)
* [Install from sources](#install-from-sources)
* [Build and run without installation](#build-and-run-without-installation)
* [Troubleshooting](#troubleshooting)


## Install with Conda
[Conda](https://docs.conda.io/en/latest/) is a cross platform package / environment management system.
Porechop_ABI was added as a package to the [bioconda](https://bioconda.github.io/) channel to make the installation easier.

With a fresh environment, you may install Porechop_ABI using either the Anaconda GUI, or the following command lines:

```bash
conda create -n my_env
conda activate my_env
conda install porechop_abi -c bioconda
porechop_abi -h
```

See Conda [command list](https://docs.conda.io/projects/conda/en/latest/commands.html) for information about managing conda packages and environments.


## Fetching source code

All the following methods require you to either download and unzip this repository, or simply clone it:

```bash
git clone https://github.com/bonsai-team/Porechop_ABI.git
cd Porechop_ABI
````


## Install with pip
If you are not using Conda but have pip installed, we recommand using pip to install Porechop_ABI, as it is the second easiest way.

In the Porechop_ABI folder, use pip as a python module.
```bash
python3.[X] -m pip install .
porechop_abi -h
```
Be sure the python version match the requirement (>= 3.6 and networkx installed)


## Install from sources
You can also run the `setup.py` script yourself. It will compile the C++ components of Porechop_ABI and install a `porechop_abi` executable on your system.
Depending on your permissions, you can install porechop_abi systemwide (root / sudo most likely  required), or for yourself only.

After cloning the repository:

* System wide install
```bash
python3.[x] setup.py install
porechop_abi -h
```

* Local (user) install
```bash
python3.[x] setup.py install --user
porechop_abi -h
```


## Build and run without installation

By simply running `make` in Porechop_ABI's directory, you can compile the C++ components but not install an executable. The program can then be executed by calling the `porechop-abi-runner.py` script.

```bash
cd Porechop_ABI
make release
./porechop-abi-runner.py -h
```


## Troubleshooting
* When installing from sources, you may in some cases get a "could not find cpp_function.so" or similar error. If that happens, you need to manually add execution permission to the impacted compiled objects.
`chmod +x <your_install folder>/<files that are "not found">`

* For user install (from sources), if you get a strange "can't combine user with prefix" error, read [this](http://stackoverflow.com/questions/4495120).

* You may want / need to install Porechop_ABI to a specific location. You can specify the installation path using --prefix:  `python3 setup.py install --prefix=$HOME/.local`

* If you'd like to specify which compiler to use, set the `CXX` environnemnt variable before installation. Ex: `export CXX=g++-6`

* Porechop includes `ez_setup.py` for users who don't have [setuptools](https://pypi.python.org/pypi/setuptools) installed, though that script is [deprecated](https://github.com/pypa/setuptools/issues/581). So if you run into any installation problems, make sure setuptools is installed on your computer: `python3.[x] -m pip install setuptools`

* If you have already the original Porechop installed, Porechop_ABI *will not* replace it. Despite having the same abilities, Porechop_ABI is a different software. It must be run using the `porechop_abi` command, and not just `porechop`.
