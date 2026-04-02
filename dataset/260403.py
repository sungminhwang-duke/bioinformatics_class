#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# 사람 인유두종바이러스(Human Papillomavirus)의 L1 유전자의 염기서열은 HPV 위험성을 판단하는 주요 기준으로 사용된다.
# 이 유전자에 의해 생성되는 단백질은 HPV 캡시드의 약 80%를 구성하며 바이러스의 숙주 세포 침입을 돕는다.

# Entez의 E-utilities를 사용하여 NCBI protein 데이터베이스에서 HPV L1 단백질에 대한
# 다섯 개의 레코드를 GenBank 포맷과 FASTA 포맷으로 가져오는 Biopython 코드를 작성하시오.


# In[ ]:


### 2026-04-03 ###

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

# biopython packages

# https://biopython.org/docs/latest/api/Bio.html
    


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


# Assignment

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




