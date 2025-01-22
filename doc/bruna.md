Rotifer session 15/01/2025
==========================  

```python
#importing rotifer and rodolfo's function
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import rodolfo as rdar
#getting the protein by its accesion
protein = rdbs.sequence['WP_xxxx']
#running a psiblast search
blast = rdar.psiblast(protein)
#running psiblast against uniref50, to get alphafold hits
af = rdar.psiblast(protein, aln=False,db='/databases/fadb/uniprot/uniref50')
```

```python
#demonstration from 20250122
#loading a neighborhood table
ndf = pd.read_pickle('ncbindff.pickle')
#loading neighborhood cursor
gnc = ncbi.GeneNeighborhoodCursor()
#loading Eduardo's function
from rotifer.devel.alpha import epsoares as rdae
#updating the taxonomy
ndfu = rdae.update_lineage(ndf)
#checking pfam annotation
ndfu.pfam.value_counts()
#running a full annotate followed by a taxonomy update
fa = rdae.update_lineage(rdar.add_arch_to_df(gnc.fetchall(sequence.tolist())))
#getting the distribution of queries through the neighborhoods
ndf.groupby('block_id')['query'].sum().value_counts()
#getting the normalized distribution
ndf.groupby('block_id')['query'].sum().value_counts(normalize=True)
#inspecting a randomly picked neighborhood
ndfu[ndfu.block_id.isin(ndfu.block_id.sample())]
#compacting the neighborhoods to better visualize the data
ndfc = ndfu.compact()
#getting the 30 most recurrent archtectures
ndfu.pfam.value_counts().head(30)
#getting the 30 most recurrent neighbors
ndfu.query('query == 0').pfam.value_counts().head(30)
```
