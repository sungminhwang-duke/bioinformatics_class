### 2026-06-01 FBA ###

#pip install cobra
#pip install cameo

=================================================================

import numpy as np
import cobra
import os
from os.path import join

from pprint import pprint
from cobra.io import validate_sbml_model
report = validate_sbml_model(filename='iMM904.xml') # (2009) https://pubmed.ncbi.nlm.nih.gov/19321003/
pprint(report)

=================================================================
# load a model

from cobra.io import read_sbml_model
model= read_sbml_model('iMM904.xml')
model

# These reations are not found from Bigg: 
# ALCD2x_copy1, ALCD2x_copy2 ,CERH124_copy1, CERH126_copy1, CERH124_copy2, CERH126_copy2, GCC2cm_copy1, GCC2cm_copy2

=================================================================
# check what reactions/metabolites/genes are present in the model

from cameo import fba
import pandas as pd

# reactions from iMM904
model_rxn = model.reactions
model_rxn_df = pd.DataFrame(model_rxn, columns= ['Reactions'])
#model_rxn_df.to_csv('iMM904_rxn.csv', index = False)

# metabolites from iMM904
model_met = model.metabolites
model_met_df = pd.DataFrame(model_met, columns= ['Metabolites'])
#model_met_df.to_csv('iMM904_met.csv', index = False)

# genes from iMM904
model_gene = model.genes
model_gene_df = pd.DataFrame(model_gene, columns=['Genes'])
#model_gene_df.to_csv('iMM904_genes.csv', index = False)                                                

=================================================================
# check media composition
model.medium

=================================================================
# generate in silico strain
#BY4741: MATa his3Δ1 leu2Δ0 met15Δ0 ura3Δ0

from cobra.io import read_sbml_model
model_iMM904_BY4741 = read_sbml_model('iMM904.xml')
model_iMM904_BY4741.solver

#model_iMM904_BY4741.genes.YOR202W.knock_out() #HIS3: if the gene is deleted, cells cannot grow.
model_iMM904_BY4741.genes.YLR303W.knock_out() #MET17, BY4741
model_iMM904_BY4741.genes.YCL018W.knock_out() #LEU2, BY4741
model_iMM904_BY4741.genes.YEL021W.knock_out() #URA3, BY4741
model_iMM904_BY4741.reactions.XYLR.knock_out() # xylose reductase (YHR104W) pathway was deleted!! h_c + nadph_c + xyl__D_c --> nadp_c + xylt_c

=================================================================
# modify medium for the simulation 

medium_iMM904_BY4741 = model_iMM904_BY4741.medium

medium_iMM904_BY4741["EX_ura_e"]= 1000    ### uracil[e] =>
medium_iMM904_BY4741["EX_ade_e"]= 0   # adenine[e] => 
medium_iMM904_BY4741["EX_cit_e"]= 0   # citrate[e] => 
medium_iMM904_BY4741["EX_gly_e"]= 0   # L-glycine[e] => 
medium_iMM904_BY4741["EX_ala__L_e"]= 0   # L-alanine[e] => 
medium_iMM904_BY4741["EX_arg__L_e"]= 0   # L-arginine[e] => 
medium_iMM904_BY4741["EX_asn__L_e"]= 0   # L-asparagine[e] => 
medium_iMM904_BY4741["EX_asp__L_e"]= 0   # L-aspartate[e] => 
medium_iMM904_BY4741["EX_cys__L_e"]= 0   # L-cysteine[e] => 
medium_iMM904_BY4741["EX_glu__L_e"]= 0   # L-glutamate[e] => 
medium_iMM904_BY4741["EX_gln__L_e"]= 0   # L-glutamine[e] => 
medium_iMM904_BY4741["EX_ile__L_e"]= 0   # L-isoleucine[e] => 
medium_iMM904_BY4741["EX_leu__L_e"]= 1000   ### L-leucine[e] => 
medium_iMM904_BY4741["EX_lys__L_e"]= 0   # L-lysine[e] => 
medium_iMM904_BY4741["EX_met__L_e"]= 0   # L-methionine[e] => 
medium_iMM904_BY4741["EX_phe__L_e"]= 0   # L-phenylalanine[e] => 
medium_iMM904_BY4741["EX_pro__L_e"]= 0   # L-proline[e] => 
medium_iMM904_BY4741["EX_ser__L_e"]= 0   # L-serine[e] => 
medium_iMM904_BY4741["EX_thr__L_e"]= 0   # L-threonine[e] => 
medium_iMM904_BY4741["EX_trp__L_e"]= 0   # L-tryptophan[e] => 
medium_iMM904_BY4741["EX_tyr__L_e"]= 0   # L-tyrosine[e] => 
medium_iMM904_BY4741["EX_val__L_e"]= 0   # L-valine[e] => 
medium_iMM904_BY4741["EX_ribflv_e"]= 0   # riboflavin[e] => 
medium_iMM904_BY4741["EX_pnto__R_e"]= 0.78   # (R)-pantothenate[e] => 
medium_iMM904_BY4741["EX_4abz_e"]= 0.78   # 4-aminobenzoate[e] => 
medium_iMM904_BY4741["EX_btn_e"]= 0.78   # biotin[e] => 
medium_iMM904_BY4741["EX_inost_e"]= 0.78   # myo-inositol[e] => 
medium_iMM904_BY4741["EX_nac_e"]= 0.78   # nicotinate[e] => 
medium_iMM904_BY4741["EX_thmpp_e"]= 0.78   # thiamine[e] => 
medium_iMM904_BY4741["EX_nh4_e"]= 1000   # ammonium[e] <=> 
medium_iMM904_BY4741["EX_fe2_e"]= 1000   # iron(2+)[e] <=> 
medium_iMM904_BY4741["EX_pi_e"]= 1000   # phosphate[e] <=> 
medium_iMM904_BY4741["EX_k_e"]= 1000   # potassium[e] <=> 
medium_iMM904_BY4741["EX_na1_e"]= 1000   # sodium[e] <=> 
medium_iMM904_BY4741["EX_so4_e"]= 1000   # sulphate[e] <=> 
medium_iMM904_BY4741["EX_h2o_e"]= 1000  # H2O
medium_iMM904_BY4741["EX_h_e"]= 1000  #H+
#not present: medium_iMM904_BY4741["r_2028"]= 0.78   # pyridoxine[e] => 
#not present: medium_iMM904_BY4741["r_1792"]= 0   # folate[e] => 
#not present: medium_iMM904_BY4741["r_4593"]= 1000   # chloride[e] <=> 
#not present: medium_iMM904_BY4741["r_4597"]= 1000   # Mg(2+)[e] <=> 
#not present: medium_iMM904_BY4741["r_4600"]= 1000   # Ca(2+)[e] <=> 
medium_iMM904_BY4741["EX_o2_e"]= 1000   # oxygen[e] <=>      # if anaerobic, set 0

#1) Glucose only - 20 g <br>
#2) Acetate only - 20 g <br>

medium_iMM904_BY4741["EX_glc__D_e"]= 0   # D-glucose[e] <=> 
medium_iMM904_BY4741["EX_ac_e"]=     0    # acetate[e] => 

model_iMM904_BY4741.medium = medium_iMM904_BY4741
model_iMM904_BY4741.medium


from cameo import fba
import pandas as pd

fba_result_iMM904_BY4741 = fba(model_iMM904_BY4741)
fba_result_iMM904_BY4741
#fba_result_iMM904_BY4741.data_frame
#fba_result_iMM904_BY4741.data_frame.to_csv('./Downloads/BY4741_1_glu.csv')
#fba_result_iMM904_BY4741.data_frame.to_csv('./Downloads/BY4741_2_acet.csv')
=================================================================
# visaulization
import ipywidgets as widgets
from IPython.display import display
import cobra
import escher

solution = model_iMM904_BY4741.optimize()
escher.Builder('iMM904.Central carbon metabolism', reaction_data = solution.fluxes)

=================================================================
# single gene deletion
import pandas
from time import time
from cobra.io import load_model
from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion)

single_gene_deletion(model_iMM904_BY4741, model_iMM904_BY4741.genes[:20])
=================================================================










### 2026-05-29 ###

# check methods such as kegg_info, ...

from Bio.KEGG import REST
dir(REST) 


# In[ ]:


from Bio.KEGG import REST

kegg = REST.kegg_info('kegg')
print(kegg.read())


# In[ ]:


# check pathway information
# https://rest.kegg.jp/info/pathway
from Bio.KEGG import REST

info = REST.kegg_info('pathway')
print(info.read())


# In[ ]:


# check 'human' pathway information
from Bio.KEGG import REST

hsa_pathways = REST.kegg_list("pathway", "hsa")
print(hsa_pathways.read())


# In[ ]:


# check 'human' pathway information and count how many pathways in human from KEGG DB
from Bio.KEGG import REST

hsa_pathways = REST.kegg_list("pathway", "hsa")
hsa_pathways_text = hsa_pathways.read()  # Call read only ONCE

# Now count
pathway_list = hsa_pathways_text.strip().split("\n")
print(f"Number of pathways: {len(pathway_list)}")


# In[ ]:


# check the specific pathway
from Bio.KEGG import REST

entry = REST.kegg_get('hsa00020')
print(entry.read())


# In[ ]:


# pathway image generation
from Bio.KEGG import REST
from IPython.display import Image

img = REST.kegg_get('hsa00020', 'image')
Image(img.read())


# In[ ]:


# Find query in KEGG

from Bio.KEGG import REST

Pd = REST.kegg_find("pathway", "Parkinson+disease")
print(Pd.read())


# In[ ]:


# Map image

from Bio.KEGG import REST
from IPython.display import Image

refPD = REST.kegg_get('map05012', 'image')
Image(refPD.read())


# In[ ]:


# Map image

from Bio.KEGG import REST
from IPython.display import Image

humanPD = REST.kegg_get('hsa05012', 'image')
Image(humanPD.read())


# In[ ]:


# Find query (drug for Parkinson) in KEGG

from Bio.KEGG import REST

Pd_drug = REST.kegg_find("pathway", "antiparkinson+agents")
print(Pd_drug.read())


# In[ ]:


from Bio.KEGG import REST

Pd_drug_check = REST.kegg_get('map07057')
print(Pd_drug_check.read())


# In[ ]:


from Bio.KEGG import REST

dopa = REST.kegg_get('D00059')
print(dopa.read())


# In[ ]:


# Chemical image

from Bio.KEGG import REST
from IPython.display import Image

dopa_chem = REST.kegg_get('map07057', 'image')
Image(dopa_chem.read())


# In[ ]:


# Chemical image

from Bio.KEGG import REST
from IPython.display import Image

dopa_chem_specific = REST.kegg_get('D00059', 'image')
Image(dopa_chem_specific.read())


# In[ ]:


# 연습문제
# Biopython을 사용하여 질병(disease)인 용혈성 요독 증후군 (hemolytic uremic syndrome)에 대한 다음 정보를 확인하시오!
#1) Entry 번호
#2) Pathway 번호






### 2026-05-22 ###

# 단백질 3차 구조 및 PDB 데이터 처리





### 2026-05-15 ###


# In[1]:


# 1) 질의어를 입력받아 UniProtKB DB에서 검색하는 함수 정의
# https://www.uniprot.org/

import requests
import json

def search_uniprot(query):
    url = "https://rest.uniprot.org/uniprotkb/search" 

    params = {'query': query,          # 질의어
              'format': 'json',        # 데이터 반환형식
              'fields': 'accession',   # 검색 대상 (예, accession number)
              'size': 10}              # 최대 검색 개수
    
    response = requests.get(url, params = params)
    
    if response.status_code == 200:    # https://incodom.kr/Status_code
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        return None


# In[6]:


dic_acc = search_uniprot("P53 human")


# In[ ]:


dic_acc


# In[ ]:


list_acc = dic_acc['results']

accession = []
for acc in list_acc:
    accession.append(acc['primaryAccession'])


# In[ ]:


list_acc
accession


# In[ ]:





# In[ ]:


# 2) PROSITE (functional domain, motif) DB에서 검색
# PS로 시작하는 등록번호
# https://prosite.expasy.org/

import requests
from Bio import SeqIO
from io import StringIO

# (1) UniProt에서 단백질 서열 가져오기 함수
def get_uniprot_sequence(protein_id):
    uniprot_url = f"https://www.uniprot.org/uniprot/{protein_id}.fasta"
    response = requests.get(uniprot_url)
    
    if response.status_code == 200:
        fasta_data = response.text
        seq_record = SeqIO.read(StringIO(fasta_data), "fasta")
        return str(seq_record.seq)
    
    
# (2) PROSITE에서 domain 검색하는 함수
def search_prosite(sequence):
    prosite_scan_url = "https://prosite.expasy.org/cgi-bin/prosite/PSScan.cgi"
    params = {"seq": sequence,
              "output": "json"}
    response = requests.post(prosite_scan_url, data=params)
    
    if response.status_code == 200:
        result = response.json()
        domains = []
        
        if "matchset" in result:
            for match in result["matchset"]:
                domains.append({
                    "Prosite ID": match["signature_ac"],
                    "start": match["start"],
                    "stop": match["stop"]
                })
            return domains


# In[ ]:


protein_id = "Q496J9"
sequence = get_uniprot_sequence(protein_id)
domains = search_prosite(sequence)
domains


# In[ ]:


# PROSITE에서 세부 서열패턴 확인-1
from Bio import ExPASy
from Bio.ExPASy import Prosite

handle = ExPASy.get_prosite_raw('PS00217')
record = Prosite.read(handle)
handle.close()

print(record.pattern)


# In[ ]:


# PROSITE에서 세부 서열패턴 확인-2
in_handle = ExPASy.get_prodoc_entry('PDOC00217')
html = in_handle.read()
in_handle.close()
with open("./Downloads/prodocrecord.html", "w") as out_handle:
    out_handle.write(html)


# In[ ]:





# In[ ]:


# 3) 질의어를 입력받아 STRING DB에서 검색하는 함수
# https://string-db.org/
import requests

def get_string_id(uniprot_id):
    url = "https://string-db.org/api/json/get_string_ids"
    params = {"identifiers": uniprot_id,
             "species": 9606} #NCBI taxonomy ID of human
    response = requests.get(url, params=params)
    
    if response.status_code == 200:
        data = response.json()
        if data:
            string_id = data[0]['stringId']
            return string_id


# In[ ]:


uniprot_id = "P69905"
string_id = get_string_id(uniprot_id)
print(f"STRING ID for {uniprot_id}: {string_id}")


# In[ ]:


# STRING에서 단백질의 네트워크 이미지 확인 함수
import requests

#string_id = "9606.ENSP00000251595"
url = f"https://string-db.org/api/image/network?identifiers={string_id}"
print(f"STRING ID for {uniprot_id}: {string_id}")

response = requests.get(url)
if response.status_code == 200:
    with open("./Downloads/network.png", "wb") as file:
        file.write(response.content)
    print(f"Image was created.")


# In[ ]:



#문제 1: “DNA polymerase E. coli”를 검색하여 accession number 5개 출력해 보세요.

#문제 2: 아래 UniProt ID를 이용하여 STRING ID 및 단백질 네트워크 이미지를 확인해 보세요.
        TP53 (P04637), EGFR (P00533), BRCA1 (P38398)





### 2026-05-08 ###

# In[ ]:
# Read newick format

from Bio import Phylo

tree = Phylo.read("./Downloads/9-sample_tree3.nwk","newick")
print(type(tree))
print(tree)


# In[ ]:


# Draw a tree

from Bio import Phylo

tree = Phylo.read("./Downloads/9-sample_tree1.nwk","newick") # (A, B, C);
tree = Phylo.read("./Downloads/9-sample_tree2.nwk","newick") # (A:0.1, B:0.3, C:0.2);
tree = Phylo.read("./Downloads/9-sample_tree3.nwk","newick") # (A, B, (C, D));

Phylo.draw(tree)


# In[ ]:


# Draw a color tree

from Bio import Phylo

tree = Phylo.read("./Downloads/9-sample_tree3.nwk","newick")

tree.rooted = True
tree.root.color = (128,128,128)
print(tree)
print("tree.clade[0]:", tree.clade[1])
print("tree.clade[1]:", tree.clade[1])
print("tree.clade[2,0]:", tree.clade[2,0])
print("tree.clade[2,1]:", tree.clade[2,1])
tree.clade[1].color = "blue"
tree.clade[2,0].color = "red"
Phylo.draw(tree)


# In[ ]:


# Draw a tree with length

from Bio import Phylo

tree = Phylo.read("./Downloads/9-sample_tree4.nwk","newick")
Phylo.draw(tree)


# In[ ]:


# Draw a tree with length showed

from Bio import Phylo

tree = Phylo.read("./Downloads/9-sample_tree4.nwk","newick")
Phylo.draw(tree, branch_labels = lambda c: c.branch_length)


# In[ ]:







### 2026-04-17 ###

# In[ ]:

# 1) MUSCLE can be operated by biopython! First, download the package: https://www.drive5.com/muscle/

from Bio.Align.Applications import MuscleCommandline 

muscle_exe = "./Downloads/muscle-win64.v5.3"  
cmd_line = MuscleCommandline(muscle_exe, input="7-MSA.fasta", out="7-MSA.aln", clw=" ") 
print(cmd_line) 


# In[ ]:


# 2-1) online MUSCLE: https://www.ebi.ac.uk/jdispatcher/
# 2-2) read by biopython!

from Bio import AlignIO 

alignment = AlignIO.read("./Downloads/7-MSA.aln","clustal") 
print(alignment) 


# In[ ]:


#separate the ID info and seuquence

from Bio import AlignIO 

alignment = AlignIO.read("./Downloads/7-MSA.aln","clustal") 
for record in alignment: 
    print("%s - %s" % (record.seq, record.id))


# In[ ]:


from Bio import AlignIO 

alignment = AlignIO.read("./Downloads/7-MSA.aln","clustal") 
for record in alignment: 
    print("%s - %s" % (record.seq[0:10], record.id))


# In[ ]:


# weblogo - 1) online: https://weblogo.threeplusone.com/


# In[ ]:


# weblogo - 2-1) biopython with generating sequence

from Bio.motifs import Motif 
from Bio import motifs 
from Bio.Seq import Seq 
from IPython.display import Image


seqs = [Seq("TACAA"), 
        Seq("TACGC"), 
        Seq("TACAC"), 
        Seq("TACCC"), 
        Seq("AACCC"), 
        Seq("AATGC"), 
        Seq("AATGC"), 
        ]

m = motifs.create(seqs) 
print(m.counts)
Motif.weblogo(m,'./Downloads/7-weblogo.png') # for saving

Image("./Downloads/7-weblogo.png")


# In[ ]:


# weblogo - 2-2) biopython with the aligned file

from Bio import AlignIO, motifs
from Bio.motifs import Motif
from Bio.Seq import Seq 

alignment = AlignIO.read("./Downloads/7-MSA.aln","clustal") 
instance = [] 
for record in alignment: 
    s = Seq(str(record.seq)) 
    instance.append(s) 
m = motifs.create(instance) 

Motif.weblogo(m,'./Downloads/7-weblogo.png') # for saving
Image("./Downloads/7-weblogo.png")


# In[ ]:


import pandas as pd
import logomaker
import matplotlib.pyplot as plt
from Bio import AlignIO

alignment = AlignIO.read("./Downloads/7-MSA.aln", "clustal")
seqs = [str(record.seq) for record in alignment]

df = pd.DataFrame([list(seq) for seq in seqs])
aa = list("ACDEFGHIKLMNPQRSTVWY")
counts_df = pd.DataFrame(0, index=aa, columns=range(df.shape[1]))

for col in df.columns:
    for aa_letter in df[col]:
        if aa_letter in counts_df.index:
            counts_df.loc[aa_letter, col] += 1

logo_df = counts_df.T
logomaker.Logo(logo_df, shade_below=.5, fade_below=.5)
#plt.title("WebLogo: Protein MSA")
#plt.tight_layout()
#plt.savefig("./Downloads/HBA_weblogo_logomaker.png", dpi=300)
#plt.show()


# In[ ]:


import pandas as pd
import logomaker
from Bio import AlignIO

alignment = AlignIO.read("./Downloads/7-MSA.aln", "clustal")
seqs = [str(record.seq) for record in alignment]

df = pd.DataFrame([list(seq) for seq in seqs])
aa = list("ACDEFGHIKLMNPQRSTVWY")
counts_df = pd.DataFrame(0, index=aa, columns=range(df.shape[1]))

for col in df.columns:
    for aa_letter in df[col]:
        if aa_letter in counts_df.index:
            counts_df.loc[aa_letter, col] += 1

logo_df = counts_df.T


# In[ ]:


counts_df
df
seqs
logo_df







### 2026-04-10 ###

#>buccal_swab.unmapped1
#CTTTTGTTAATCGATGATATACAGTCACTCAGCGGAAAAAAAGTCGCAACTCAGGAAGAATTTTTCAATACCTTTAACGCCCTTCATG


#>buccal_swab.unmapped2
#CCAGCCCCCCAGCCTCCCGATCACGGTTTACTACGCCGTGTTGGAGCGCGCCTGCCGCAGCGTGCTCCTAAACGCACCGTCGGAGGCCCCCCAGATTGTCCGC


# In[ ]:


from Bio import Entrez
Entrez.email = "hwangs@kmou.ac.kr"

handle = Entrez.efetch(db="nucleotide", id="CP046379", rettype="gb", retmode="text")
record = handle.read()


# In[ ]:


print(record[0:1000])


# In[ ]:


### 1) with a fasta file

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW 
from Bio import SeqIO 

record = SeqIO.read("./Downloads/6-buccal_swab_unmapped1.fasta", format="fasta") 
handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta")) # https://biopython.org/docs/dev/api/Bio.Blast.NCBIWWW.html


# In[ ]:


print(handle.read())


# In[ ]:


### 2) with DNA sequence

from Bio.Seq import Seq

fasta_file = Seq("CTTTTGTTAATCGATGATATACAGTCACTCAGCGGAAAAAAAGTCGCAACTCAGGAAGAATTTTTCAATACCTTTAACGCCCTTCATG")


# In[ ]:


print(fasta_file)


# In[ ]:


from Bio.Seq import Seq
from Bio.Blast import NCBIWWW 
from Bio import SeqIO 

fasta_file = Seq("CTTTTGTTAATCGATGATATACAGTCACTCAGCGGAAAAAAAGTCGCAACTCAGGAAGAATTTTTCAATACCTTTAACGCCCTTCATG")

handle1 = NCBIWWW.qblast("blastn", "nt", str(fasta_file), format_type="XML") #type 1: XML
handle2 = NCBIWWW.qblast("blastn", "nt", str(fasta_file), format_type="HTML") #type 2: HTML


# In[ ]:


# type 1: XML

with open("./Downloads/6-blast_results.xml", "w") as output_file:
    output_file.write(handle1.read())
    
handle1.close()


# In[ ]:


from Bio import SearchIO

blast_qresult = SearchIO.read('./Downloads/6-blast_results.xml', 'blast-xml')
print(blast_qresult)

#HSP: high-scoring pair


# In[ ]:


# type 2: HTML

with open("./Downloads/6-blast_results.HTML", "w") as output_file:
    output_file.write(handle2.read())
    
handle2.close()


# In[ ]:


print(blast_qresult[0])

# E-value: the number of expected hits of similar quality (score) that could be found just by chance.
# 값이 작으면 작을수록, 서열 매치가 우연히 발생할 확률이 낮다!

# E = m x n  / 2bit-score         (https://www.metagenomics.wiki/tools/blast/evalue)
# m: query sequence length       (https://bio-kcs.tistory.com/entry/BLAST-BLAST-%EC%95%8C%EA%B3%A0%EB%A6%AC%EC%A6%98%EC%97%90-%EB%8C%80%ED%95%B4-%EC%95%8C%EC%95%84%EB%B3%B4%EC%9E%90)
# n: total database length (sum of all sequences)
# bit-score: a normalized score derived from the raw alignment score (S) using the scoring system
# BLAST의 Raw Score(S)를 정규화한 값. Raw Score(S)는 비교하는 서열의 길이와 치환행렬(BLOSUM)에 따라 달라지는데, 이를 보정하기 위해 Bit Score (S')가 사용됨


# In[ ]:


### 3) with amino acid sequence

from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW
from Bio.Seq import Seq

fasta_file = Seq("MKFIEEIVVDAFLPTFRALLAEDLRDRGFTQSEVAEALGISQSAVSKYAHGEVATNERVATDPRVVDLVSRVGDGLATGDMTPVQALVEAEVLIRQLEEGDLLSDLHEDEMPELASHDGFRSIHDPEGRLRTVEQVRSSVRRGLRMLTNTSGFAGLIPNVGSNLVESLPDADSVDDVAAIPGRIFDVKGQATVPGEPEFGVSGHVAGVLLSARAAGADVNAALNIVYDAGVIEDLEAAGYECIEFDPDAPTDPVRELLTARDLPETFVVYQSGGYGIEPITYILGPDAPAVADVVRVLL")

#handle = NCBIWWW.qblast("blastn", "nt", str(fasta_file), format_type="XML")  # for DNA sequence
handle_p = NCBIWWW.qblast("blastp", "nr", str(fasta_file), format_type="XML") # database= nr (non-redundant), swissprot, refseq_protein



# In[ ]:


with open("./Downloads/6-blastp_results.xml", "w") as output_file:
    output_file.write(handle_p.read())
    
handle_p.close()
print("BLASTp is completed and the result is saved as 'blastp_results.xml'! ")

blastp_qresult = SearchIO.read('./Downloads/6-blastp_results.xml', 'blast-xml')
print(blastp_qresult)
print(blastp_qresult[0])


# In[ ]:




### 2026-04-03 ###


# In[ ]:


# 사람 인유두종바이러스(Human Papillomavirus)의 L1 유전자의 염기서열은 HPV 위험성을 판단하는 주요 기준으로 사용된다.
# 이 유전자에 의해 생성되는 단백질은 HPV 캡시드의 약 80%를 구성하며 바이러스의 숙주 세포 침입을 돕는다.

# Entez의 E-utilities를 사용하여 NCBI protein 데이터베이스에서 HPV L1 단백질에 대한
# 다섯 개의 레코드를 GenBank 포맷과 FASTA 포맷으로 가져오는 Biopython 코드를 작성하시오.

from Bio import Entrez
Entrez.email ="hwangs@kmou.ac.kr"

handle = Entrez.esearch(db="protein", term="human papillomavirus AND L1", RetMax=5)
res = Entrez.read(handle)
handle.close()
print(res)

res_ids = res['IdList']
handle_fas = Entrez.efetch(db='protein', id=res_ids, rettype="fasta", retmode='txt')
L1_fas = handle_fas.read()
print(L1_fas)





# In[ ]:

class 클래스이름:
    def 메서드이름(self):
        명령블록


# In[ ]:


class Student:
    def major(self):
        print("My major is")


# In[ ]:


# example

class Student:
    def __init__(self, major):
        self.major = major
    def say(self):
        print(f'My major is {self.major}') # character string(문자열)


# In[ ]:


Bioinfo = Student('Bioinformatics')
Bioinfo.say()


# In[ ]:


Eng = Student('Engineering')
Eng.say()


# In[ ]:





# In[ ]:


class Student:
    def __init__(self, major, year):
        self.major = major
        self.year = year
    def say(self):
        print(f'My major is {self.major} and {self.year} year student.' )


# In[ ]:


Bioinfo = Student('Bioinformatics', 4)
Bioinfo.say()


# In[ ]:


Eng = Student('Engineering', 2)
Eng.say()


# In[ ]:





# In[ ]:


# Bio.seq exercise (https://biopython.org/docs/1.78/api/Bio.Seq.html)

from Bio.Seq import Seq

tatabox_seq = Seq(" tataaaggcAATATGCAGTAG")


# In[ ]:


print(tatabox_seq)


# In[ ]:


tatabox_seq.count("A")


# In[ ]:


tatabox_seq.count("a")


# In[ ]:


dir(tatabox_seq)


# In[ ]:


tatabox_seq.lower()


# In[ ]:


tatabox_seq.upper()


# In[ ]:


print(tatabox_seq)
print(tatabox_seq.strip())


# In[ ]:


# GC contents
from Bio.Seq import Seq

gc = Seq("ATGCATGCATGC")
g_count = gc.count("G")
c_count = gc.count("C")
gc_content = (g_count+c_count)/len(gc)*100


# In[ ]:


gc_content


# In[ ]:


# transcription & translation
from Bio.Seq import Seq

DNA = Seq("ATGAACTAAGTTTAGAAT")


# In[ ]:


mRNA = DNA.transcribe()


# In[ ]:


mRNA


# In[ ]:


AA = DNA.translate()


# In[ ]:


AA


# In[ ]:


AA = DNA.translate(to_stop=True)


# In[ ]:


AA


# In[ ]:


from Bio.Data import CodonTable

codon_table = CodonTable.unambiguous_dna_by_name["Standard"]


# In[ ]:


print(codon_table)


# In[ ]:


# complement & reverse complement

 # 5'- TATAAAGGCAATATGCAGTAG -3'
 # 3'- ATATTTCCGTTATACGTCATC -5'

from Bio.Seq import Seq

seq = Seq("TATAAAGGCAATATGCAGTAG")


# In[ ]:


comp_seq = seq.complement()


# In[ ]:


comp_seq


# In[ ]:


rev_comp_seq = seq.reverse_complement()


# In[ ]:


rev_comp_seq


# In[ ]:


# Bio.SeqUtils module

# GC contents
#gc = Seq("ATGCATGCATGC")
#g_count = gc.count("G")
#c_count = gc.count("C")
#gc_content = (g_count+c_count)/len(gc)*100
# 50%

from Bio.Seq import Seq
from Bio.SeqUtils import GC

exon_seq = Seq("ATGCATGCATGC")
gc_content = GC(exon_seq)

gc_content


# In[ ]:


# Molecular weight
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

seq_mw = Seq("ATGCAGTAG")
molecular_weight(seq_mw)


# In[ ]:


# Translation combination
from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations

trans_six = Seq("ATGCCTTGAAATGTATAG")
six_frame_translations(trans_six)


# In[ ]:


print(six_frame_translations(trans_six))


# In[ ]:


# Tm
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

tm_seq = Seq("AGTCTGGGACGGCGCGGCAATCGCA")

print(mt.Tm_Wallace(tm_seq))
print(mt.Tm_GC(tm_seq))
print(mt.Tm_NN(tm_seq))


# In[ ]:


# AA abbreviations
from Bio.Seq import Seq
from Bio.SeqUtils import seq1

aa_3 = "LeuLysMetValIleThrTrpPhe"
seq1(aa_3)


# In[ ]:


from Bio.SeqUtils import seq3

aa_1 = "LKMVITWF"
seq3(aa_1)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


### 2026-03-27 ###

 
# In[ ]:


from Bio import Entrez
Entrez.email = "자기 이메일 주소"
# 예를 들어, Entrez.email = "hwangs@kmou.ac.kr"


#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi
#https://www.ncbi.nlm.nih.gov/books/NBK25497/

handle = Entrez.einfo()
record = handle.read()


# In[ ]:


print(record)


# In[ ]:


from Bio import Entrez
Entrez.email = "hwangs@kmou.ac.kr"

handle = Entrez.einfo()
record = Entrez.read(handle)


# In[ ]:


#https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi

print(record)
print(len(record["DbList"]))


# In[ ]:





# In[ ]:


from Bio import Entrez
Entrez.email = "hwangs@kmou.ac.kr"

handle = Entrez.esearch(db="pubmed", term="haloferax AND thiN", RetMax=10)
record = Entrez.read(handle)
print(record)
print(record["Count"]) #https://pubmed.ncbi.nlm.nih.gov/


# In[ ]:


from Bio import Entrez
Entrez.email = "hwangs@kmou.ac.kr"

#https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
handle = Entrez.efetch(db="nucleotide", id="NC_002058.3", rettype="gb", retmode="text")


# In[ ]:


print(handle.read())


# In[ ]:


from Bio import Entrez
from Bio import SeqIO
Entrez.email = "hwangs@kmou.ac.kr"

#https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
handle = Entrez.efetch(db="nucleotide", id="NC_002058.3", rettype="gb", retmode="text")
record = SeqIO.parse(handle, "gb")
seq_lists = list(record)

handle.close()
type(seq_lists)

seq_lists[0:4]


# In[ ]:


seq_rec = seq_lists[0]
fasta_format = ">%s\n%s\n"%(seq_rec.id, seq_rec.seq)

# %s: This is a placeholder for a string (text data).
# \n: This represents a newline character (line feed).


# In[ ]:


seq_rec


# In[ ]:


fasta_format


# In[ ]:


print(fasta_format)


# In[ ]:


from Bio import Entrez
from Bio import SeqIO
Entrez.email = "hwangs@kmou.ac.kr"

handle = Entrez.esearch(db="nucleotide", term="thiN", RetMax=3)
record = Entrez.read(handle)
print(record)
record['IdList']

ids = record['IdList']
handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
records = handle.read()
print(records)


# In[ ]:





# In[ ]:




# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


### 2026-03-20 ###

# Scopus: https://www.scopus.com/
# pybliometrics: Python library to pull, cache and extract data from the Scopus database.
#                https://pybliometrics.readthedocs.io/en/stable/index.html
# pip install pybliometrics

import pybliometrics


# In[ ]:


# Scopus API generation at https://dev.elsevier.com/
# Use my API: 810f8e333bce1dee177bcf93d7c79051

pybliometrics.scopus.init()


# In[ ]:


# Document-specific information

from pybliometrics.scopus import AbstractRetrieval

ab = AbstractRetrieval("10.1128/mbio.00633-22")
ab.title
#ab.authors


# In[ ]:


# Terms for searching / Compare the website version https://www.scopus.com/

from pybliometrics.scopus import ScopusSearch

query = ' TITLE-ABS-KEY ( bioplastic  AND bacteria  OR  archaea  AND technology ) '
s = ScopusSearch(query,
                 download = True, # saving the results
                 verbose = True)  # current process


# In[ ]:


s


# In[ ]:


# How many papers?

s.get_results_size()


# In[ ]:


# Make database from the paper info

import pandas as pd

df_s = pd.DataFrame(s.results)


# In[ ]:


df_s


# In[ ]:


# DB check

df_s.head()


# In[ ]:


df_s.shape


# In[ ]:


# DB details
s0 = df_s.loc[37]
s0


# In[ ]:


# Data formatting




df_pubyear
#df_pubyear.shape


# In[ ]:


df_s["aggregationType"].unique()


# In[ ]:


df_s["year"] = df_s["coverDate"].apply(lambda x: x.split("-")[0])


# In[ ]:


df_s["year"] = df_s["year"].astype(int)


# In[ ]:


df_s["num_pub"] = [1] * df_s.shape[0]


# In[ ]:


df_pubyear = df_s.query("aggregationType == 'Journal'").groupby("year").sum()


# In[ ]:


df_pubyear


# In[ ]:


df_pubyear.to_csv("df_pubyear.csv")


# In[ ]:


#Plot using DB

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

# Load the dataset
df_p1 = pd.read_csv("df_pubyear.csv")

# Extract relevant columns and handle missing values
df_p2 = df_p1[["year", "num_pub"]].dropna()

# Convert to appropriate data types
df_p2["year"] = df_p2["year"].astype(int)
df_p2["num_pub"] = df_p2["num_pub"].astype(int)

# Sort by year
df_p2 = df_p2.sort_values(by="year")

# Plot the data
plt.figure(figsize=(10, 5))
plt.plot(df_p2["year"], df_p2["num_pub"], marker="o", linestyle="-")
plt.xlabel("Year")
plt.ylabel("Number of Publications")
plt.title("Number of Publications per Year")
plt.grid(True)
plt.show()


# In[ ]:




