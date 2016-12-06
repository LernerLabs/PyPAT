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
