# Com a sequencia passada para todos vocês, o objetivo inicial é fazer análise de cluster e análise de cluster 

# Registre seus resultados de todas as maneiras que puder:
# 1. source.sh (tosource) e ipython_log.py (%logstart)
# 2. Na sua versao desse arquivo, com comentarios e bem organizado
# 3. Em um documento para apresentacao e consulta, que pode ser um doc (Word, Google Doc, Latex, etc) ou apresentacao (Powerpoint, Google Presentation ou Canvas)
# 4. Gere figuras e anotacoes para representar suas conclusoes

#########################################################################################################################################
# Comandos na shell

# Esta parte precisará ser feita no terminal, visto que tem mais recursos do que se fosse executado com Rotifer no Python
# Aqui é apenas um guia inicial, mas existe a liberdade de variar os parâmetros e as análises de modo geral.
# Caso haja alguma mudança, escrever as mudanças neste pipeline e deixar anotado a razão de porque variou o parâmetro ou a análise

# Após conferir, criar o diretório (ou pasta) que serão feitas as realizações das tarefas de bioinformática.

# No diretório, abra o gedit ou outro programa qualquer, coloque a sequência objeto de estudo que foi passado em formato FASTA e salve como "sequence.fa"
# Aqui, neste procedimento, será feito a busca de sequência usando a sequência da yfdq fornecida pelo Davi. No entanto, nas suas análises, substitua a sequência que foi passada a você e faça o mesmo passo, seguindo o mesmo procedimento

>yfdq
MDKSAIEAIQDSSAAAAKATQEQIPAGITHLVAVPNGVTLQNVEISLAGRTRFRGKLVTSSIPD
FVTYVKNREGGHGFIDTDKLGATVFFNLGDADKPGHADHTARLTLQATPAYLAMLAAGGKAFSQ
RDALDFIEDWSHLIGAHQVSDDGAFSPIPLARAIAAIRKVKIKAESTTEQTEGNFQQSRSVLES
VAASSDEGLPDVLSFRTEPYLGLPERTFHMRVSVSTTGREPSLRFRVIALEQEQDEIAKQFKAL
LLGEVGEAASMTIGTFTP

# Busca de sequencias usando psiblast. Digite no terminal:

psiblast -query sequence.fa -out yfdq_psiblast.out -num_threads 4 -num_iterations 5 -inclusion_ethresh 0.001 -evalue 0.01 -db /databases/fadb/nr/nr50

# Após fazer a busca das sequências no psiblast, é necessário converter em formato de .TSV para fazermos as análises no Python. 
# Para isso, execute no terminal:

blast2table -y yfdq_psiblast.out > yfdq_psiblast.tsv

# O mesmo procedimento será feito no jackhmmer. 
# Busque sequências homólogas usando jackhmmer (número de iterações = 5, e-value <=0.001)

jackhmmer -N 5 --incdomE 0.1 --incE 0.001 --cpu 4 sequence.fa /databases/fadb/nr/nr50 > sequence.jackhmmer.out

# O arquivo de saída do jackhmmer precisa também ser convertido em .TSV. Para isso, execute no terminal: 

hmmer2table -a sequence.jackhmmer.out > yfdq_jackhmmer.tsv

#########################################################################################################################################

# A partir daqui, tudo será feito no Python. 
# Inicialmente, será feito:

# 1) Juntar os dados das buscas no psiblast e jackhmmer
# 2) Fazer análise de agrupamento
# 3) Análise de vizinhança para investigar o contexto gênico da proteína de interesse

# Para trabalharmos no Python, é necessário abrir um terminal e carregar o ambiente

module load rotifer mambaforge tass

conda activate py310

conda activate py310 ; tosource

# Observação: se quiser tornar esses módulos permanentes, basta escrever, logo após carregá-los, no terminal: 

module save

# Para conferir se os módulos estão carregados, digite o comando no terminal: 

module list

# Nesse comando, espera-se encontrar, dentre outros módulos (se estiver carregado), os módulos tass, mambaforge, rotifer, nas quais serão importantes para trabalharmos no Ipython3.

# Digite no terminal:

ipython3

# Isso abrirá o Python, na qual faremos toda a análise. Digite:

import io
import os
import sys
import base64
import pickle
import plotly
import importlib
import numpy as np
from Bio import SeqIO
import networkx as nx
import dash_cytoscape as cyto
from dash import Dash, html
import pandas as pd
import networkx as nx
import urllib.request as urlreq
import matplotlib.pyplot as plt
import dash
import dash_bio as dashbio
from dash import Dash, dcc, html, Input, Output, callback
from dash.dependencies import Output, Input
import rotifer
from rotifer.devel.beta import sequence as rdbs
from rotifer.devel.alpha import rodolfo as rdar
from rotifer.devel.alpha import gian_func as gf
import rotifer.pandas.functions as rpf
from rotifer.pandas import functions as rpf
from rotifer.interval import utils as riu
from rotifer.db import ncbi
from rotifer.db.ncbi import entrez

rpf.print_everything()

# Isso carregará todos os módulos e bibliotecas úteis. Vamos considerar os objetivos descritos no início, passo-a-passo. 

# 1) Juntar os dados das buscas no psiblast e jackhmmer

# Neste caso, será usado o ambiente Python para juntar os dados. Caso queira juntar os dados de outra forma, não tem problema.
# Porém, caso não faça o procedimento, é necessário é necesário descrever todo o procedimento no pipeline e entregar.
# Isso é importante porque o procedimento será refeito para garantir a reprodutibilidade do procedimento

# Digite (o eval adiciona a coluna source para sabermos de onde veio o dado):

psiblast_yfdq = pd.read_csv("yfdq_psiblast.tsv", sep='\t').eval('source = "psiblast"')
jackhmmer_yfdq = pd.read_csv("yfdq_jackhmmer.tsv", sep='\t').eval('source = "jackhmmer"')

# Renomear colunas para que correspondam entre os DataFrames (se necessário)

psiblast_yfdq.rename(columns={
    'hit_name': 'sequence',
    'hit_start': 'start',
    'hit_end': 'end',
    'query_name': 'model',
    'query_coverage': 'cov',
    'query_start': 'qstart',
    'query_end': 'qend',
    'bits': 'score',
    'hit_length': 'talilen'
}, inplace=True)

# Selecionar colunas correspondentes para combinar

cols = ['sequence', 'model', 'start', 'end', 'evalue', 'cov', 'qstart', 'qend', 'iteration', 'score', 'talilen', 'source']

# Manter apenas as colunas que têm correspondência

psiblast_yfdq = psiblast_yfdq[cols]
jackhmmer_yfdq = jackhmmer_yfdq[cols]

# Combinar os DataFrames

merged_yfdq = pd.concat([psiblast_yfdq, jackhmmer_yfdq])

# Neste caso, foram obtidas 6116 proteínas em 27758 linhas da tabela. Observe que houve redundância. É necessário remover a redundância. Antes de remover a redundância, é necessário salvar os dados combinados. Para isso:

merged_yfdq.to_csv("combined_yfdq.tsv", sep='\t', index=False)

# Eliminar as linhas duplicadas baseadas na coluna 'sequence'

merged_yfdq.sort_values(['sequence','evalue'], ascending=True, inplace=True)
merged_dropdupl_yfdq = merged_yfdq.drop_duplicates(subset='sequence', keep='first')

# Esse resultado gerou 6116 proteínas não redundantes

# Com isso, nós finalizamos a primeira etapa dos objetivos descritos. Agora é necessário fazer a análise de agrupamento.

# 2) Fazer análise de agrupamento

# Para fazer a análise de agrupamento, faça os procedimentos abaixo usando o módulo add_cluster:

merged_dropdupl_yfdq_filtered = merged_dropdupl_yfdq.query('evalue <= 0.001').copy()
yfdqaln = rdbs.sequence(merged_dropdupl_yfdq_filtered.sequence.tolist())

# Analise das sobreposicoes (multiplos alinhamentos do mesmo par query <-> hit)

# A analise pode continuar apenas produzindo esse dado mas tenha em mente que o hit
# que estiver em merged_yfdq_multicopy precisaria, em principio, ser fragmentado e ter
# duas linhas no alinhamento da familia, uma para cada fragmento
#
# Se quiser fazer isso, venha discutir o codigo ou procedimento conosco
#
# Sem a fragmentacao das duplicacoes, o programa de align() pode ficar confuso

oldLogLevel = rotifer.logger.getEffectiveLevel()
rotifer.logger.setLevel(rotifer.logging.INFO)
merged_yfdq_multicopy = riu.filter_nonoverlapping_regions(merged_dropdupl_yfdq_filtered, reference=['sequence'], start='start', end='end', criteria={'evalue':True,'region_length':False})
rotifer.logger.setLevel(oldLogLevel)
del(oldLogLevel)
multicopy = merged_yfdq_multicopy.sequence.value_counts().where(lambda x: x > 1).dropna().index.tolist()
if multicopy:
    merged_yfdq_multicopy = merged_yfdq_multicopy[merged_yfdq_multicopy.sequence.isin(multicopy)].copy()
else:
    merged_yfdq_multicopy = pd.DataFrame(columns=cols) # Dataframe vazio se nao tiver nada duplicado

# Antes de ver o alinhamento, fique atento ao numero de sequencias:
# So execute os proximos comandos se tiver menos de 4000 sequencias!!!!!
# Caso contrario, prossiga com o clustering

yfdqaln = yfdqaln.align()
yfdqaln.view()

# Adicionando clusters (grupos por similaridade)

yfdqaln.add_cluster(cascade=[(100,100)] + [ (0,x) for x in range(80,-1,-10) ], inplace=True)

# Salve o objeto modificado em formato pickle

yfdqaln.df.to_pickle("MSA_all-clusters.pkl")

# Aqui, é necessário fazer uma análise de grupos. No procedimento abaixo, vamos considerar apenas o agrupamento baseado no critério c0i80 e também a filtragem de sequências que apresentaram e-value menores que 0.001

######################################################
#                                                    #
#      Analisando o critério de agrupamento c80e3    #
#                                                    #
######################################################

# Antes de proceguirmos, é sempre importante observar os grupos formados usando o critério c100i100 para reduzir a redundância.
 
frequencia_absoluta_c100i100 = yfdqaln.df['c100i100'].value_counts() 

# Neste exemplo, foram obtidos 4848 grupos e mostrados os seus representantes. Neste caso essa quantidade corresponde ao total de proteínas
# Posteriormente, vamos determinar a quantidade de grupos na coluna c0i80

frequencia_absoluta_c0i80 = yfdqaln.df['c0i80'].value_counts() 

# Nesta agrupamento foi obtido de 754 grupos

# Vamos analisar a quantidade de grupos usando o critério c0i80

#tabela_frequencia = pd.DataFrame({'c0i80': frequencia_absoluta_c0i80.index, 'Frequência Absoluta': frequencia_absoluta_c0i80.values}) 
tabela_frequencia = frequencia_absoluta_c0i80.reset_index().rename({'index':'c0i80', 'c0i80':'Frequência Absoluta'}, axis=1)

# Os três primeiros grupos que contém mais membros neste exemplo foram: 

# EJK0928014.1: 3338; MBW6805787.1: 309; WP_198262334.1: 72. 

# Esses três grupos representam um total de 3719 de 4848, ou seja, 77 % dos dados não redundantes. Salve todas as informações (todos os grupos) deste resultado

tabela_frequencia.to_csv("frequencia_c80e3.tsv", sep = "\t")

# Uma vez que temos a informação de cada um dos representantes de cada grupo, podemos analisar separadamente grupo por grupo e ver as características que cada grupo possui e que são exclusivas destes (se for). 

# Inicialmente será feito a análise do cluster representado por EJK0928014.1 contendo 3338 proteínas
# Para fazer a análise de cada grupo, é necessário olhar separadamente cada um. Para isso:

#df = yfdqaln.df

# Filtre as linhas que contêm o grupo c0i80 com o maior numero de sequencias
# Esse grupo pode ser identificado com %page tabela_frequencia ou com

yfdq_cluster_maior = tabela_frequencia.head(1).c0i80[0]

# A linha abaixo funciona mas e importante inspecionar a tabela_frequencia
#filtered_df = df[df['c0i80'] == yfdq_cluster_maior]

# Remover colunas que estão em yfdq_cluster_maior_aln antes do merge
#filtered_df = filtered_df.drop(columns=['sequence', 'length', 'type', 'description'])

# Passo 2: Alinhar as sequências
yfdq_cluster_maior_aln = yfdqaln.filter(f'''c0i80 == "{yfdq_cluster_maior}"''').align()
#yfdq_cluster_maior_aln = rdbs.sequence(filtered_df.id.tolist()).align()

# Observe neste exemplo que houve uma melhora significativa no alinhamento. Existem sequências de tamanhos variados e outras contendo falsos positivos
# Para verificar a distribuição do tamanho das sequências, rode o comando:

yfdqaln.to_hist() 
yfdq_cluster_maior_aln.to_hist() 

# Passo 3: importar colunas do add_Cluster usando merge
df = yfdqaln.df.filter(['id','c100i100'] + [ f'c0i{x}' for x in range(80,-1,-10) ]).drop_duplicates()
yfdq_cluster_maior_aln.df = yfdq_cluster_maior_aln.df.merge(df, on='id', how='left')

# A saída desta análise é:

# Total proteins: 3338
# count sequence size
###############################################################################
██                                                                33  38 - 64  
██████████████████                                               212  64 - 90  
███████████████████████████████████████                          454  90 - 115 
████████████████████████████████                                 370  115 - 141
████████████████████████████████████████████████████████         652  141 - 167
██████████████████████████████████████████████                   541  167 - 193
███████████████████████████████████████████████████████████████  727  193 - 219
██████████████████████████████                                   347  219 - 244
                                                                   1  244 - 270
                                                                   1  270 - 296
###################################################################################

# Podemos inspecionar as sequências que possuem entre 38 e 64 resíduos aqueles que possuem entre 244 e 296 e verificar se são de fato proteínas homólogas à yfdq. Se forem, elas precisam ser mantidas no cluster, caso contrário, precisam ser removidas dos dados.

# Nesta distribuição, é sempre interessante estudar o alinhamento. Se houver grandes quantidades de gaps em determinadas porções da proteína, verifique se são falsos positivos usando a predição da estrutura no AlphaFold ou no EMS (https://esmatlas.com/resources?action=fold) seguida da análise no FoldSeek 

# Dessa forma, deverá ser feito essa inspeção de modo a limpar o resultado dos falsos positivos
# Neste exemplo, a sequência vou citar "EHX98206.1" que está no cluster e que é falso positivo. Outro exemplo "KYS01234.1", na qual o foldseek não identificou nenhum homólogo. Para eliminar essas sequências do dataframe, basta fazer o procedimento:

#yfdq_cluster_maior_aln_rm = result_df[~result_df['id'].isin(['EHX98206.1', 'KYS01234.1'])]
aligned_sequences = yfdq_cluster_maior_aln.filter(remove=['EHX98206.1', 'KYS01234.1']).trim(99.9)

# Visualizar o alinhamento
aligned_sequences.view()

# Salve o alinhamento, para usar como entrada para predição da estrutura 3D no AlphaFold

aligned_sequences.to_file("yfdq_cluster_maior_aln_rm.a3m", output_format='a3m')

# Faça a mesma análise para os outros dois grupos. Salve os resultados em diretórios diferentes, renomeando as pastas: cluster_analysis dentro desta pasta deve ter as pastas com o nome do representante_quantidade_de_membros. Por exemplo:

# mkdir cluster_analysis
# cd cluster_analysis
# mkdir EJK0928014.1_3338
# mkdir MBW6805787.1_309
# mkdir WP_198262334.1_72

# Em cada pasta salve os resultados das análises de grupos. Essa análise e resultado será importante para fazer a análise de vizinhança.

# Se por acaso o comprimento estiver associado com as sequências falso positivas. Por exemplo, vamos supor que as sequências que são correspondentes à DUF2303 aplique:

# Seleciona apenas as linhas do dataframe que contem na coluna length valores maiores ou iguais a 70 e menores que 246
filtered_sequences = aligned_sequences.filter('length >= 70 or length < 246').trim(99.9)

# Inspecione o alinhamento caso faça essa filtragem
filtered_sequences.view()

# Ficara ainda melhor usando apenas os representantes do cluster c0i50 e realinhando:
#filtered_sequences.filter('id == c0i50').align().view()

# Salve o alinhamento e o dataframe
#filtered_sequences.filter('id == c0i50').align().to_csv("dados_cluster_filtrado_por_tamanho.tsv", sep = "\t")
pickle.dump(filtered_sequences, open("dados_cluster_filtrado_por_tamanho.pkl", 'wb'))

# Somente use esse critério, caso observe a necessidade deste tipo de filtragem.
#########################################################################################################################################

# 3) Análise de vizinhança para investigar o contexto gênico da proteína de interesse

# Após fazer todas as inspeções, salve o objeto na variável "aligned_sequences". A partir disso, fazer a busca dos genomas e as anotações.

# Ignore as linhas comentadas e vai para o rdar.full_annotate
#genomas = pd.read_csv("lista_de_genomas.tsv", sep = '\t')
#ic = entrez.IPGCursor()
#gnc = ncbi.GeneNeighborhoodCursor(mirror='/databases/genomes')
#ipgs = ic.fetchall(aligned_sequences.df.id.tolist())
#ipgs['prefered'] = ipgs.assembly.isin(genomas.assembly)
#aligned_sequences.ndf = gnc.fetchall(aligned_sequences.df.id.tolist(), ipgs=ipgs)
#aligned_sequences.ndf = rdar.add_arch_to_df(aligned_sequences.ndf)

rdar.full_annotate(aligned_sequences, mirror='/databases/genomes')

# Para trabalhar com tabelas no Pandas ou outro software (como R), será conveniente trabalhar com formatos TSV. Para salvar os resultados brutos, basta digitar o comando abaixo:

aligned_sequences_ndf = aligned_sequences.ndf

aligned_sequences_ndf.to_csv("gene-neighborhood.tsv", sep ="\t")

# Fazer uma análise geral na coluna compact

aligned_sequence_compact = aligned_sequences_ndf.compact()

# Selecionar as linha que contém "DUF2303->YfbR-like" na coluna "compact"
selected_rows = aligned_sequence_compact['compact'].str.contains(r'DUF2303->YfbR-like')

# Filtrando as linhas que atendem ao critério
filtered_df = aligned_sequence_compact[selected_rows]


################################################################################################################################

# Análise de vizinhos mais frequentes

msa = aligned_sequences

co_ocorrencia = gf.cluster_Co_occurrence(msa.ndf, freq_cutoff=0.1)
co_ocorrencia2 = msa.ndf.compact()
co_ocorrencia2.to_csv("compact_groupby_c80e_c80i70_cluster_Co_occurrence_annotation_pfam.csv", sep='\t')

# Função para criar a rede de grafos direcionada
def to_network(df, target=['pfam'], ftype='CDS', interaction=True, ignore = []):
    if isinstance(target, str):
        target = [target]

    w = df.query(f'type == "{ftype}"')[['block_id', 'strand']].copy()
    w['source'] = df[target[0]]

    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), df[df.type == ftype][col], w['source'])

    w['source'] = w['source'].fillna("?").str.split('+')
    w = w.explode(column='source')
    w['target'] = w['source'].shift(-1)
    
    if ignore:
    	w = w[~(w.source.isin(ignore) | w.target.isin(ignore))].copy()
    	
    w = w[(w.block_id == w.block_id.shift(1)) & (w.strand == w.strand.shift(1))].copy()
    w.loc[w.strand == -1, ['source', 'target']] = w.loc[w.strand == -1, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where((w.index.to_series() == w.index.to_series().shift(1)), 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(weight=('strand', 'count'), blocks=('block_id', 'nunique')).reset_index()
    return w

#filtered_ndf = msa.ndf[~msa.ndf['pfam'].str.contains('SP|TM', na=False) & ~msa.ndf['aravind'].str.contains('SP|TM', na=False)]

# Cria a tabela de frequência contendo a contagem par a par de proteínas vizinhas

results_networks = to_network(aligned_sequences_ndf, target=['pfam', 'aravind', 'c80e3'], ignore=['SP', 'TM']).sort_values(['weight'], ascending=False)
results_networks.to_csv("cluster_Co_occurrence.tsv", sep="\t", index=False)

