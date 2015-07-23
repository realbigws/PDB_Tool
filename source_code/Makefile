

# make the subdir first
CXX = g++
FLAG = -O3 -static 
SHELL :=/bin/bash

PDBToolDir = ./PDB_Tool_src
PDBToolEXE = /PDB_Tool

all:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../

PDB_Tool:
	cd $(PDBToolDir)&&make
	mv ${PDBToolDir}$(PDBToolEXE) ../


# clean all the temp files
clean:
