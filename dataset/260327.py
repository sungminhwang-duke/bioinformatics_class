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


# Assignment 과제 제출 : 4/2(목) 16시까지!! hwangs@kmou.ac.kr 보내기

사람 인유두종바이러스(Human Papillomavirus)의 L1 유전자의 염기서열은 HPV 위험성을 판단하는 주요 기준으로 사용된다. 이 유전자에 의해 생성되는 단백질은 HPV 캡시드의 약 80%를 구성하며 바이러스의 숙주 세포 침입을 돕는다. 

Entez의 E-utilities를 사용하여 NCBI protein 데이터베이스에서 HPV L1 단백질에 대한 다섯 개의 레코드를 GenBank 포맷과 FASTA 포맷으로 가져오는 Biopython 코드를 작성하시오.




