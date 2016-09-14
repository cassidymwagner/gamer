#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_* End_*                        ../../../src/Init/
cp Flu_*                               ../../../src/Fluid/
cp Par_*                               ../../../src/Particle/
cp plot*.gpt                           ../../../bin/Run/
cp Input__Flag_NParPatch               ../../../bin/Run/
cp Input__TestProb                     ../../../bin/Run/

cp Makefile.ELBDM                      ../../../src/Makefile
cp Input__Parameter.ELBDM              ../../../bin/Run/Input__Parameter