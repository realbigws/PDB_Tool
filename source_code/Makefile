

# make the subdir first
CXX = g++
#FLAG = -O3 -static 
FLAG = -O3
SHELL :=/bin/bash

PDBToolDir = ./PDB_Tool_src
PDBToolEXE = /PDB_Tool

all: PDB_Tool

PDB_Tool:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../


# clean all the temp files
clean:
