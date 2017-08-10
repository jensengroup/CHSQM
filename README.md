# CHSQM
Prediction of labile carbon hydrogens using semiempirical methods

CHSQM consists of two bash scripts, `CHSQM1` and `CHSQM2`. CHSQM1 and CHSQM2 call other python scripts and bash scripts and calls OpenBabel. Some of the python scripts use RDKit. The method makes use of the semiempirical quantum chemistry program MOPAC.

The very first time you use the scripts type `chmod +x CHSQM1 CHSQM2`

Try example:

    ./CHSQM1 example pm3_mop rm_proton

"example" refers to the text file `example.smiles` and "pm3_mop" refer to the text file `pm3_mop.header`. You can also use rm_hydride and rm_hydrogen instead of rm_proton

CHSQM creates a folder called example_rm_proton+pm3_mop with the mopac input files and submits to a slurm queue using the bash script `submit_bactches_mopac`.

Once all calculations finish:

    ./CHSQM2 example_rm_proton+pm3_mop

This creates example_rm_proton+pm3_mop_1.pka and example_rm_proton+pm3_mop_3.pka which contains the atom number(s) of the carbon atoms with the most labile hydrogen and the name of the mopac output file for the corresponding isomer using a 1 and 3 pH unit cutoff.

highlight_atoms.py creates svg files with pictures of the molecules where atom centers listed in the .pka file are highlighted. These two svg files are then merged into a single file using merge_svg.py, which also changes the colors of the highlighted atoms.  Note that highlight_atoms.py reguires the installation of ipython.

The resulting svg file (marked combined) can be opened in a browser or converted to pdf by third party software
