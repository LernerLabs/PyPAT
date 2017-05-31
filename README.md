# Current Notes

The documentation is not currently being updated. It is not *known* to be meaningfully out of date, but no promises are made at the moment.

## Memory/speed issues with correlated dynamics code

Right now, running the correlated dynamics code on even a relatively small number of files is taking up a ton of memory. Way more than it should. So, maybe time to restructure.

First, a note on the fies. step5_1 only has 35 frames in it, so should be ignored. The others have 500 frames, and are comparable.

Version 1, which looks like

```
trajin ./step5_2.nc
trajin ./step5_3.nc
trajin ./step5_4.nc
trajin ./step5_5.nc

strip @H*
strip :WAT
center :1-160 mass origin
image origin center
rms first *

matrix correl out ptraj_files/structure/one_all_atom_correlmat.dat start 100 stop 300 byatom
matrix correl out ptraj_files/structure/one_byres_correlmat.dat start 100 stop 300 byres
average  ptraj_files/structure/one_average.pdb start 100 stop 300 pdb
run
```

runs in 6.2 seconds and takes under 1 percent of memory.

Version 2, which has the same preface, but the following matrix/run lines

```
matrix correl out ptraj_files/structure/two_all_atom_correlmat.dat start 200 stop 400 byatom
matrix correl out ptraj_files/structure/two_byres_correlmat.dat start 200 stop 400 byres
average  ptraj_files/structure/two_average.pdb start 200 stop 400 pdb
run
```

Does the same.

Version 12, which combines the two

```
trajin ./step5_2.nc
trajin ./step5_3.nc
trajin ./step5_4.nc
trajin ./step5_5.nc

strip @H*
strip :WAT
center :1-160 mass origin
image origin center
rms first *

matrix correl out ptraj_files/structure/12one_all_atom_correlmat.dat start 100 stop 300 byatom
matrix correl out ptraj_files/structure/12one_byres_correlmat.dat start 100 stop 300 byres
average  ptraj_files/structure/12one_average.pdb start 100 stop 300 pdb
run

matrix correl out ptraj_files/structure/12two_all_atom_correlmat.dat start 200 stop 400 byatom
matrix correl out ptraj_files/structure/12two_byres_correlmat.dat start 200 stop 400 byres
average  ptraj_files/structure/12two_average.pdb start 200 stop 400 pdb
run
```

Gets up to 18% of memory (5.7g), and runs for significantly longer.

More troubling, the output files for 12one... are the same, but the
ones for 12two are different! Looking at the output and output files,
it looks like the strip, center, image, rms commands did *not* get
applied to the second set of commands! That's kind of a huge problem!

Indeed, putting the strip commands in fixes everything!



# And now the documentation

Documentation may be found in the Doc/ subdirectory. Please refer to
either PyPAT.pdf or the html documentation for installation and usage
instructions.

Quick guide to installing the scripts:

```bash
pip install -e .
```

Quick guide to rebuilding the documentation:

```bash
make clean
make html
```

(optionally)
```bash
make latex
cd ./build/latex
make all-pdf
```

Quick guide for making correlated dynamics movies similar to the documentation:

(Make sure you have structure.parm7 structure.pdb and step5_1.nc in your directory)

```bash
mkdir ptraj_files
mkdir ptraj_files/images
mkdir ptraj_files/structure
~/coding/PyPAT/drivers/write_ptraj_input_files.py --strip-water --strip-hydros --input-dir=. --mdcrd=step5_1.nc --structure-name=structure --start=500 --stop=5500 --ps=5 --window-size=1000 --window-spacing=500
~/coding/PyPAT/drivers/run_ptraj.py --input-dir=. --prmtop=structure.parm7 --structure=structure
~/coding/PyPAT/drivers/make_correlated_dynamics_plots.py --structure-name=structure --start=500 --stop=5500 --window-size=1000 --window-spacing=500
~/coding/PyPAT/drivers/make_movies.py --structure-name=structure --plot-types=ca,avg,max,min,abs,straight,mainheavy
```

```bash
module load amber/16

mkdir ptraj_files
mkdir ptraj_files/images
mkdir ptraj_files/structure

#export START=500
#export STOP=2500
#export START=75500
#export STOP=90500
#export START=500
#export STOP=15500
#export START=500
#export STOP=100500

export START=500
export STOP=10500
export WSIZE=1000
export WSPACE=500

export TRAJS=''
for i in `seq 1 10`; do export TRAJS=$TRAJS,step5_$i.nc; done
export TRAJS=${TRAJS:1}
echo $TRAJS

~/coding/PyPAT/drivers/write_ptraj_input_files.py --strip-water --strip-hydros --other-ptraj-strips=":CLA,:SOD" --input-dir=. --mdcrd=$TRAJS --structure-name=structure --start=$START --stop=$STOP --ps=5 --window-size=$WSIZE --window-spacing=$WSPACE
~/coding/PyPAT/drivers/run_ptraj.py --input-dir=. --prmtop=structure.parm7 --structure=structure
~/coding/PyPAT/drivers/make_correlated_dynamics_plots.py --structure-name=structure --start=$START --stop=$STOP --window-size=$WSIZE --window-spacing=$WSPACE
~/coding/PyPAT/drivers/make_movies.py --structure-name=structure --plot-types=ca,avg,max,min,abs,straight,mainheavy

mv ptraj_files ptraj_files_$START\_$STOP\_$WSIZE
```
