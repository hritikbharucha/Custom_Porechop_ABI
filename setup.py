#!/usr/bin/env python3
"""
Run 'python3 setup.py install' to install Porechop.
"""

# Make sure this is being run with Python 3.4 or later.
import sys
if sys.version_info.major != 3 or sys.version_info.minor < 6:
    print('Error: you must execute setup.py using Python 3.6 or later')
    sys.exit(1)

import os
import shutil
from distutils.command.build import build
from distutils.core import Command
import subprocess
import multiprocessing
import fnmatch
import importlib.util


def printerr(msg):
    print(msg, file=sys.stderr)


def get_choice(msg):
    choice = "_"
    while(choice not in "yYnN"):
        choice = input(msg + " (y / n)\n")
    return(choice in "yY")


def install_networkx():
    # Install networkx with the same python executable for compatibility.
    install_command = [sys.executable, "-m", "pip", "install", "networkx"]
    print("We will try to install networkx using pip:")
    print(" ".join(install_command))
    agreed = get_choice("Do you want to install networkx now?")
    if(agreed):
        try:
            subprocess.check_call(install_command)

        except CalledProcessError as e:
            printerr("Could not install networkx:")
            printerr(e)
            exit(1)
    else:
        print("Networkx will not be installed.")
        print("Try again when networkx is installed for your version of python.")
        print(sys.executable)
        exit(0)


# Install setuptools if not already present.
if not importlib.util.find_spec("setuptools"):
    printerr("/!\\Using EZ_setup (setuptools not found)")
    import ez_setup
    ez_setup.use_setuptools()

# Checking if networkx is installed
if not importlib.util.find_spec("networkx"):
    printerr("")
    printerr("##############################")
    printerr("/!\\networkx was not found/!\\")
    printerr("##############################")
    printerr("")
    install_networkx()


from setuptools import setup
from setuptools.command.install import install

# Get the program version from another file.
exec(open('porechop/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCRIPTION = readme.read().decode()


class PorechopBuild(build):
    """
    The build process runs the Makefile to build the C++ functions into a shared library.
    """

    def run(self):
        build.run(self)  # Run original build code

        clean_cmd = ['make', 'distclean']
        try:
            make_cmd = ['make', '-j', str(min(8, multiprocessing.cpu_count()))]
        except NotImplementedError:
            make_cmd = ['make']

        # building in release mode for installation.
        make_cmd += ["release"]

        def clean_cpp():
            subprocess.call(clean_cmd)

        def compile_cpp():
            subprocess.call(make_cmd)

        self.execute(clean_cpp, [], 'Cleaning previous compilation: ' + ' '.join(clean_cmd))
        self.execute(compile_cpp, [], 'Compiling Porechop: ' + ' '.join(make_cmd))


class PorechopInstall(install):
    """
    The install process copies the C++ shared library to the install location.
    """

    def run(self):
        install.run(self)  # Run original install code
        shutil.copyfile(os.path.join('porechop', 'cpp_functions.so'),
                        os.path.join(self.install_lib, 'porechop', 'cpp_functions.so'))
        # adding approx_counter  file
        shutil.copyfile(os.path.join('porechop', 'approx_counter'),
                        os.path.join(self.install_lib, 'porechop', 'approx_counter'))

        # adding  compatibility shared library
        shutil.copyfile(os.path.join('porechop', 'compatibility.so'),
                        os.path.join(self.install_lib, 'porechop', 'compatibility.so'))

        # adding  msa consensus exec
        shutil.copyfile(os.path.join('porechop', 'msa_consensus'),
                        os.path.join(self.install_lib, 'porechop', 'msa_consensus'))

        # moving config file
        shutil.copyfile(os.path.join('porechop', 'ab_initio.config'),
                        os.path.join(self.install_lib, 'porechop', 'ab_initio.config'))


class PorechopClean(Command):
    """
    Custom clean command that really cleans up everything, except for:
      - the compiled *.so file needed when running the programs
      - setuptools-*.egg file needed when running this script
    """
    user_options = []

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd

        delete_directories = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for dir_name in fnmatch.filter(dir_names, '*.egg-info'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, 'build'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, '__pycache__'):
                delete_directories.append(os.path.join(root, dir_name))
        for delete_directory in delete_directories:
            print('Deleting directory:', delete_directory)
            shutil.rmtree(delete_directory)

        delete_files = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for filename in fnmatch.filter(filenames, 'setuptools*.zip'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.o'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.pyc'):
                delete_files.append(os.path.join(root, filename))
        for delete_file in delete_files:
            print('Deleting file:', delete_file)
            os.remove(delete_file)


setup(name='porechop',
      version=__version__,
      description='Porechop',
      long_description=LONG_DESCRIPTION,
      url='http://github.com/rrwick/porechop',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPL',
      packages=['porechop'],
      entry_points={"console_scripts": ['porechop = porechop.porechop:main']},
      zip_safe=False,
      cmdclass={'build': PorechopBuild,
                'install': PorechopInstall,
                'clean': PorechopClean}
      )
