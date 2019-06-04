Porechop *ab initio* version.

This version of porechop contains a wrapper for my [adapter finder](https://github.com/qbonenfant/adaptFinder) tool, allowing adapter inference from the reads.

It automatically integrate the infered adapter to the Porechop adapter list.

## REQUIREMENT
The requirement are the same as porechop, except you'll also need to install networkx. 
~~~
pip install networkx
~~~
Networkx is used for the simple assembly part.

## INSTALLATION of ab initio version

First clone the repository using --recursive option


```bash
git clone --recursive https://github.com/qbonenfant/Porechop.git
```

Then just install as described in the Porechop [documentation](#installation).
You can either install it or just build the c++ sources and then run locally. 

```bash
cd Porechop
./setup.py install
porechop -h
```
If there is any complains about permissions, you can either use sudo -H

```bash
sudo -H ./setup.py install
```
or just build the main files and run it using the porechop_runner.py

```bash
make clean
make
./porechop-runner.py -h
```

If you want to use the new feature, use the --ab_initio flag.

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
