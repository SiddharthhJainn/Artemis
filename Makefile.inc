# PYTHIA configuration file.
# Generated on Monday 17 January 2022 10:42:53 AM IST with the user supplied options:
# --with-rivet=/home/siddharth/Downloads/rivet/scratch/rivet/local

# Install directory prefixes.
PREFIX_BIN=/home/siddharth/Desktop/PYTHIA/newpythia/pythia8304/bin
PREFIX_INCLUDE=/home/siddharth/Desktop/PYTHIA/newpythia/pythia8304/include
PREFIX_LIB=/home/siddharth/Desktop/PYTHIA/newpythia/pythia8304/lib
PREFIX_SHARE=/home/siddharth/Desktop/PYTHIA/newpythia/pythia8304/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/home/siddharth/Desktop/artemis/fastjet/fastjet-install/bin/
FASTJET3_INCLUDE=-I/home/siddharth/Desktop/artemis/fastjet/fastjet-install/include
FASTJET3_LIB=-L/home/siddharth/Desktop/artemis/fastjet/fastjet-install/lib -Wl,-rpath,/home/siddharth/Desktop/artemis/fastjet/fastjet-install/lib -lfastjet

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=false
LHAPDF6_CONFIG=
LHAPDF6_BIN=
LHAPDF6_INCLUDE=
LHAPDF6_LIB=

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=true
RIVET_CONFIG=rivet-config
RIVET_BIN=/home/siddharth/Downloads/rivet/scratch/rivet/local/bin/
RIVET_INCLUDE=-I/home/siddharth/Downloads/rivet/scratch/rivet/local/include
RIVET_LIB=-L/home/siddharth/Downloads/rivet/scratch/rivet/local/lib -Wl,-rpath,/home/siddharth/Downloads/rivet/scratch/rivet/local/lib -lRivet

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/home/siddharth/Downloads/root-6.24.02-install/bin/
ROOT_INCLUDE=-I/home/siddharth/Downloads/root-6.24.02-install/include
ROOT_LIB=-L/home/siddharth/Downloads/root-6.24.02-install/lib -Wl,-rpath,/home/siddharth/Downloads/root-6.24.02-install/lib -lCore

GZIP_USE=false
GZIP_CONFIG=
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
