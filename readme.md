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
according to the ComBind score w.r.t. the selected poses.

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

Use the following command to prepare the structural data using Schrodinger's
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

```
combind featurize features docking/*/*_pv.maegz
```

### Pose prediction with ComBind

```
combind pose-prediction features poses.csv
```

Optionally, you can extract the poses selected by ComBind to a single file.
The resulting file will contain the protein structure followed by one pose (the
one selected by ComBind) for each ligand.

```
combind extract-top-poses poses.csv docking/*/*_pv.maegz
```

## ComBindVS

To run virtual screening using ComBindVS, you must begin with a structure of the
target protein, a set of helper ligands, and a library of compounds to screen.

The first two steps, which can be done in parallel, are to determine poses for
the helper ligands using ComBind and to produce an initial set of docked poses
for the library to be screened. Then, ComBindVS can be 

### Use ComBind to solve for poses of a set of helper ligands

Use ComBind to predict poses for the known binders and extract the selected
poses to a single file, as described above. In the below, we'll assume that this
file is named `helpers_pv.maegz`

### Dock the library to be screened

The library to be screened can be docked the same way as described above,
but here it is highly recommended that you use the `--multiplex` option during
ligprep (to write all the compounds to one file) and the `--screen` option
during docking, which will limit the number of poses per compound to 30 and
not used enhanced pose sampling.

```
combind ligprep library.csv --multiplex
combind dock ligands/library/library.maegz --screen
```

### ComBindVS

To compute the ComBind scores for each pose, we need to compute the pairwise]
features between each candidate pose to the helper ligand poses.

```
combind featurize --no-mcss --screen --max-poses 100000 features_screen docking/library-to-grid/library-to-grid_pv.maegz helpers_pv.maegz
```

With these features in hand, you can then compute the ComBind scores. The ComBind
scores for each pose will be written to the indicated numpy file (here screen.npy).

```
combind screen screen.npy features_screen
```

It is often convenient to apply the scores to the original poseviewer file and
use existing schrodinger utilities to sort the results.

```
combind apply-scores docking/library-to-grid/library-to-grid_pv.maegz screen.npy combind_scores_added_pv.maegz
$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score -o combind_pv.maegz combind_scores_added_pv.maegz
```

## Installation

Start by cloning this git repository.

ComBind requires access to Glide along with several other Schrodinger tools
and the Schrodinger Python API.

First, make sure that you have a SCHRODINGER environmental variable set
pointing to the root of the schrodinger software installation.

You can only access the Schrodinger Python API using their interpretter.
Creating a virtual environment that makes their interpretter the default
python interpretter is the simplest way to do this. To create the environment
and upgrade the relevant packages run the following:

```
$SCHRODINGER/run schrodinger_virtualenv.py schrodinger.ve
source schrodinger.ve/bin/activate
pip install --upgrade numpy sklearn scipy pandas lmdb
```

To setup the environment before each use, run
`source schrodinger.vs/bin/activate` to activate the environment and then
run `source setup.sh` to set combind specific environmental variables.
