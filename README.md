# sunRunner3D-manuscript

Repository containing all codes, scripts, and 2D boundary condition files needed to reproduce the results of the sunRunner3D manuscript.

# Cloning the Repository

To clone this repository use

```
git clone https://github.com/predsci/suRunner3D-manuscript
``` 
# Installing PLUTO

The sunRunner3D tool is based on the PLUTO code. 

Hence PLUTO must first be downloaded and installed from here:

http://plutocode.ph.unito.it/

If you are not familiar with the PLUTO code, we recommend that you take a moment to read the manual and run the examples. To understand how the boundary conditions (BC) can be defined, please refer to Chapter 5 of the PLUTO manual.


# Installing pysunrunner

We use the [pysunrunner](https://github.com/predsci/pysunrunner) package to visualize the results of sunRunner3D. You will need to clone and install the package using:

```
git clone https://github.com/predsci/pysunrunner
pip install -e pysunrunner
```
We recommend that you setup a conda or miniconda environment for the package.

# Heliospheric relaxation

You will start by relaxing the heliospheric domain using the boundary conditions we provide.  Start by navigating to your local copy of the relaxation directory

```
cd relax_dir
```
The `init.c` and `userdef_output.c` in this directory are the only two PLUTO routines we have modified. Together with the `definitions.h` and `pluto.ini` files they define your run configuration.
Although this directory includes the `Makefile` we use, it is safer to generate a `Makefile` that is specific to your computer architecture. This can be done by invoking the  PLUTO `setup.py` script
```
python $PLUTO_DIR/setup.py
```
as explained in Section 1.3 of the PLUTO manual. The `2D` subdirectory includes the boundary files loaded at each timestep by `init.c`.  Once the Makefile is created PLUTO can be compiled:
```
make
```
We recommend compiling PLUTO with MPI support and running it on multiple threads, e.g.,
```
mpirun -np 32 ./pluto -no-x3par 1>log 2>err &
```

The relaxation wall clock time will depend on your platform specifications and the number of threads you use. All output files will be generated in the `output` subdirectory, which is currently empty.
If you change the name of the `output` directory, be sure to update the output_dir line in the `pluto.ini` file.

To track the progression of the relaxation you can use:

```
tail output/pluto.0.log
```

# CME Run
Once the heliospheric relaxation is complete you can proceed to the CME run step. 
Start by navigating to the `cme_dir` directory
```
cd cme_dir
```
here too you see the two PLUTO files we have modified: `init.c` and `userdef_output.c'. Note that the `init.c` in this directory includes the CME and is not the same as the one in the `relax_dir` 
You will need to create your platform-specific Makefile,  compile PLUTO, and copy the results of the relaxation to the current directory

```
python $PLUTO_DIR/setup.py
make
cp -r ../relax_dir/output .
```
The output directory will include all the log files from the relaxation which are not needed and can be deleted:
```
cd output
/bin/rm -rf pluto*log
cd ..
```
The CME run is a continuation of the relaxation phase. You will specify this in your command line when starting the run:
```
mpirun -np 32 ./pluto -no-x3par -restart 1 1>log 2>err &
```
Here too, all your results will be generated in the `output` directory and you can track the CME run progression using:
```
tail output/pluto.0.log
```

# Generating Figures

The `fig_dir` directory contains all the scripts needed to reproduce Figures 7-12 from the manuscript. Each figure has its dedicated script. You will need to update each script with the correct path to your PLUTO CME results. 







