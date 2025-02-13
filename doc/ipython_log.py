# IPython log file

import rotifer.db.ncbi.entrez
get_ipython().run_line_magic('pinfo2', 'rotifer.db.ncbi.entrez')
import rotifer.db.ncbi.entrez as entrez
help(entrez)
from rotifer.db import ncbi
help(ncbi)
help(entrez)
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', 'projects/')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', 'calycins/')
get_ipython().run_line_magic('cd', 'work/')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('cd', '20241202/')
get_ipython().run_line_magic('ls', '')
get_ipython().run_line_magic('less', 'sequence.fa')
get_ipython().run_line_magic('cd', '..')
get_ipython().run_line_magic('cd', '20240924_osso/')
get_ipython().run_line_magic('ls', '')
from pandas as pd
import pandas as pd
osso = pd.read_csv("osso_jackhmmer.tsv", sep="\t")
osso.head()
help(entrez)
fc = entrez.FastaCursor()
fc["WP_258193852.1"]
type(fc["WP_258193852.1"])
fc["WP_258193852.1"].seq
str(fc["WP_258193852.1"].seq)
fc["WP_258193852.1"].id
len(fc["WP_258193852.1"])
print(fc["WP_258193852.1"].format("fasta"))
print(fc["WP_258193852.1"].format("genbank"))
print(fc["WP_258193852.1"].format("gbk"))
sc = entrez.SequenceCursor()
sc["WP_258193852.1"]
sc.fetchone("WP_258193852.1")
next(sc.fetchone("WP_258193852.1"))
sc.fetchall("WP_258193852.1")
osso.head()
sc["WP_220238834.1"]
sc["WP_128382277.1"]
sc["7SFP_E"]
sc["7SFP_A"]
fc["WP_258193852.1"]
osso.sequence.head(20).tolist()
osso.sequence.drop_duplicates().head(20).tolist()
fc.fetchall(osso.sequence.drop_duplicates().head(20).tolist())
from Bio import SeqIO
import sys
s = fc.fetchall(osso.sequence.drop_duplicates().head(20).tolist())
len(s)
fc.missing
type(s)
type(s[0])
[ x.id for x in s ]
fc.missing
fc.tries
SeqIO.write(s, "fasta", sys.stdin)
for x in s:
    print(SeqIO.write(x, "fasta", sys.stdin))
    
s
for x in s:
    print(SeqIO.write(x, "fasta", sys.stdin))
    
for x in s:
    print(SeqIO.write(x, sys.stdin, "fasta"))
    
help(SeqIO.write)
for x in s:
    print(SeqIO.write(x, sys.stdin, "fasta"))
    
for x in s:
    print(SeqIO.write(x, sys.stdout, "fasta"))
    
for x in s:
    SeqIO.write(x, sys.stdout, "fasta")
    
import pandas as pd
from Bio import SeqIO
from rotifer.db import ncbi
osso = pd.read_csv("osso_jackhmmer.tsv", sep="\t")
fc = entrez.FastaCursor()
for x in fc.fetchone(osso.sequence.drop_duplicates().head(20).tolist()):
    SeqIO.write(x, sys.stdout, "fasta")
    
ic = entrez.IPGCursor()
ic["WP_258193852.1"]
from rotifer.pandas import functions as rpf
rpf.print_everything()
ic["WP_258193852.1"]
ic.fetchall(osso.sequence.drop_duplicates().head(20).tolist())
type(ic.fetchall(osso.sequence.drop_duplicates().head(20).tolist()))
ic.fetchall(osso.sequence.drop_duplicates().head(20).tolist())
get_ipython().run_line_magic('page', '')
gnc = ncbi.GeneNeighborhoodCursor()
type(gnc)
help(gnc)
gnc.cursors
gnc = ncbi.GeneNeighborhoodCursor(mirror="/databases/genomes")
gnc.cursors
type(gnc.cursors)
type(gnc.cursors['ftp'])
gnc.readers
gnc.fetchall(osso.sequence.drop_duplicates().head(20).tolist())
get_ipython().run_line_magic('page', '')
import inspect
help(object.mro)
help(rotifer.db)
help(rotifer.db.core)
get_ipython().run_line_magic('pinfo2', 'rotifer.db')
get_ipython().run_line_magic('pinfo2', 'gnc')
help(gnc)
fc = ncbi.FastaCursor()
help(fc)
help(rotifer.db.core.BaseCursor)
help(gnc)
help(rotifer.db.core.BaseCursor)
help(fc)
get_ipython().run_line_magic('pinfo2', 'rotifer.db.ncbi.GeneNeighborhoodCursor')
help(ncbi)
help(ncbi)
help(ncbi.ftp)
help(ncbi)
help(ncbi)
help(ncbi.ftp)
help(fc)
help(rotifer.db.methods)
help(fc)
help(rotifer.db.delegator)
get_ipython().run_line_magic('pinfo2', 'rotifer.db.delegator')
fc.cursors
get_ipython().run_line_magic('logstart', '')
get_ipython().run_line_magic('logstop', '')
