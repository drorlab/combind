# ComBind

ComBind integrates data-driven modeling and physics-based docking for
improved binding pose prediction and binding affinity prediction.

Given the chemical structures of several ligands that can bind
a given target protein, ComBind solves for a set of poses, one per ligand, that
are both highly scored by physics-based docking and display similar interactions
with the target protein. ComBind quantifies this vague notion of "similar" by
considering a diverse training set of protein complexes and computing the
overlap between protein–ligand interactions formed by distinct ligands when
they are in their correct poses, as compared to when they are in randomly
selected poses. To predict binding affinities, poses are predicted for
the known binders using ComBind, and then the candidate molecule is scored
according to the ComBind score w.r.t the selected poses.

## Predicting poses for known binders

First, see instructuctions for software installation at the bottom of this page.

Running ComBind can be broken into several components: data curation,
data preparation (including docking), featurization of docked poses,
and the ComBind scoring itself.

Note that if you already have docked poses for your molecules of interest, you
can proceed to the featurization step. If you are knowledgable about your target
protein, you may well be able to get better docking results by manually
preparing the data than would be obtained using the automated procedure
implemented here.

### Curation of raw data

To produce poses for a particular protein, you'll need to provide a 3D structure
of the target protein and chemical structures of ligands to dock.

These raw inputs need to be properly stored so that the rest of the pipeline
can recognize them.

The structure(s) should be stored in a directory `structures/raw`.
Each structure should be split into two files `NAME_prot.mae` and `NAME_lig.mae`
containing only the protein and only the ligand, respectively.

If you'd prefer to prepare your structures yourself, save your
prepared files to `structures/proteins` and `structures/ligands`. Moreover,
you could even just begin with a Glide docking grid which you prepared yourself
by placing it in `docking/grids`.

Ligands should be specified in a csv file with a header line containing at
least the entries "ID" and "SMILES", specifying the ligand name and the ligand
chemical structure.

### Data preparation and docking

Use the following command, to prepare the structural data using Schrodinger's
prepwizard, align the structures to each other, and produce a docking grid.

```
combind structprep
```

In parallel, you can prepare the ligand data using the following command.
By default, the ligands will be written to seperate files (one ligand per file).
You can specify the `--multiplex` flag to write all of the ligands to the same
file.

```
combind ligprep ligands.csv
```

Once the docking grid and ligand data have been prepared, you can run the
docking. The arguments to the dock command are a list of ligand files to be
docked. By default, the docking grid is the alphabetically first grid present
in `structures/grids`; use the `--grid` option to specify a different grid.

```
combind dock ligands/*/*.maegz
```

### Featurization

Note that this is the 

```
combind featurize features docking/*/*_pv.maegz
```

### Pose prediction with ComBind

```
combind pose-prediction features poses.csv
```

## ComBind virtual screening

To run ComBindVS, first use ComBind to 

## Installation

Start by cloning this git repository (likely into your home directory).

ComBind requires access to Glide along with several other Schrodinger tools
and the Schrodinger Python API.

The Schrodinger suite of tools can be accessed on Sherlock by running
`ml chemistry schrodinger`. This will add many of the Schrodinger tools to
your path and sets the SCHRODINGER environmental variable. (Some tools are
not added to your path and you'll need to write out $SCHRODINGER/tool.)
After running this you should be able to run Glide by typing `glide` in the
command line.

You can only access the Schrodinger Python API using their interpretter.
Creating a virtual environment that makes their interpretter the default
python interpretter is the simplest way to do this. To create the environment
and upgrade the relevant packages run the following:

```
cd
$SCHRODINGER/run schrodinger_virtualenv.py schrodinger.ve
source schrodinger.ve/bin/activate
pip install --upgrade numpy sklearn scipy pandas

cd combind
ln -s  ~/schrodinger.ve/bin/activate schrodinger_activate
```

This last line is just there to provide a standardized way to access the
activation script.

Run `source schrodinger_activate` to activate the environment in
the future, you'll need to do this everytime before running ComBind.
This is included in the setup_sherlock script; you can source the
script by running `source setup_sherlock`.

