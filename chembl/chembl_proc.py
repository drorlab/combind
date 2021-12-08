from schrodinger.structure import StructureReader, StructureWriter
import sys
import pandas as pd
import os

csv = sys.argv[1]
smi = csv.replace('.csv', '.smi')
nonames = csv.replace('.csv', '_nonames.maegz')
noactivities = csv.replace('.csv', '_noactivities.maegz')
maegz = csv.replace('.csv', '.maegz')

df = pd.read_csv(csv)
df = df.drop_duplicates(subset=['ligand_chembl_id'])

if not os.path.exists(smi):
	df[['canonical_smiles', 'ligand_chembl_id']].to_csv(smi, sep=' ', index=False, header=False)

if not os.path.exists(nonames):
	cmd = 'ligprep -epik -ismi {} -omae {} -WAIT -HOST localhost:8'
	cmd = cmd.format(smi, nonames)
	os.system(cmd)

if not os.path.exists(noactivities):
	cmd = 'python {}/dock/ligprep.py {} {}'
	cmd = cmd.format(os.environ['COMBINDHOME'], nonames, noactivities)
	os.system(cmd)

if not os.path.exists(maegz):
	df = df.set_index('ligand_chembl_id')
	with StructureReader(noactivities) as reader, StructureWriter(maegz) as writer:
		for st in reader:
			st.property['r_chembl_activity'] = df.loc[st.title, 'standard_value']
			writer.append(st)
