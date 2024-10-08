# Chem279_HW2

In this repo, you will find code that was used to calculate analytical and numerical integration for 
shell overlap integrals. 

The numerical integration used the trapezoidal rule in order to find the overlap and the analytical integration was integration done in x, y, and z dimensions. 

This repo contains the bin/ directory where all executable files are moved to once compiled; the sample_input/ directory contains input files for both the analytical and numerical integrations; the sample_output/ directory contains output files for the analytical and numerical integrations; the src/ directory is where source files for both types of integrations are kept. In the root of the directory contains a makefile where both integrations can be run at the same time. Numerical and analytical can be ran separately with 'make P1' and 'make P2'- individually. 