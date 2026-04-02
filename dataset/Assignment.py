################################################################################################################
########################################### 김민지 2026-04-02 15:57:59 ###########################################
################################################################################################################
from Bio import Entrez, SeqIO

Entrez.email = "kmjkmj1113@naver.com"

info_handle = Entrez.einfo(db="protein")
info = Entrez.read(info_handle)
info_handle.close()

search_handle = Entrez.esearch(
    db="protein",
    term="Human papillomavirus[Organism] AND L1[Gene]",
    retmax=5
)
search_result = Entrez.read(search_handle)
search_handle.close()

ids = search_result["IdList"]

gb_handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="gb",
    retmode="text"
)
genbank_data = list(SeqIO.parse(gb_handle, "genbank"))
gb_handle.close()

SeqIO.write(genbank_data, "HPV_L1.gb", "genbank")

fasta_handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="fasta",
    retmode="text"
)
fasta_data = list(SeqIO.parse(fasta_handle, "fasta"))
fasta_handle.close()

SeqIO.write(fasta_data, "HPV_L1.fasta", "fasta")

print(ids)
print("저장 완료")


################################################################################################################
########################################### 김민석 2026-04-02 15:03:23 ###########################################
################################################################################################################
# In[ ]:
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "kimmin020514@g.kmou.ac.kr"

search_handle = Entrez.esearch(db="protein", term=search_term, retmax=5)

search_results = Entrez.read(search_handle)
search_handle.close()

id_list = search_results["IdList"]
print(id_list)


# ln[ ]: (GenBank 포맷으로 가져오기)
with Entrez.efetch(db="protein", id=id_list, rettype="gb", retmode="text") as handle:
    gb_records = list(SeqIO.parse(handle, "genbank"))

SeqIO.write(gb_records, "hpv_l1_5records.gb", "genbank")
print("GenBank 파일 저장 완료: hpv_l1_5records.gb")

print("\n[GenBank records]")
for record in gb_records:
    print(record.id, record.description)

# ln[ ]: (FASTA 포맷으로 가져오기)
with Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text") as handle:
    fasta_records = list(SeqIO.parse(handle, "fasta"))

SeqIO.write(fasta_records, "hpv_l1_5records.fasta", "fasta")
print("FASTA 파일 저장 완료: hpv_l1_5records.fasta")

print("\n[FASTA records]")
for record in fasta_records:
    print(">" + record.id)
    print(record.seq)


################################################################################################################
########################################### 황소연 2026-04-02 13:11:07 ###########################################
################################################################################################################
#PDF가 이상한데?


################################################################################################################
########################################### 최진호 2026-04-02 11:50:44 ###########################################
################################################################################################################
#pip install biopython

 

from Bio import Entrez

from Bio import SeqIO

Entrez.email = "chlwlsgh134@gmail.com“

 

handle = Entrez.esearch(db="protein", term="Human Papillomavirus L1", retmax=5)

record = Entrez.read(handle)

ids = record["IdList"]

 

handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")

print(handle.read())

 

handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")

records = SeqIO.parse(handle, "gb")

 

print(">%s\n%s" % (rec.id, rec.seq))


################################################################################################################
########################################### 한여원 2026-04-02 11:43:26 ###########################################
################################################################################################################
from Bio import Entrez
Entrez.email = "yuwon410@gmail.com"


# In[ ]:

handle = Entrez.esearch(db="protein", term="Human Papillomavirus L1", retmax=5)
record = Entrez.read(handle)
handle.close()

print(record)
ids = record["IdList"]


# In[ ]:

handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
genbank_data = handle.read()
handle.close()

print(genbank_data)

with open("HPV_L1.gb", "w") as f:
    f.write(genbank_data)


# In[ ]:

handle = Entrez.efetch(db="protein", id=ids,
rettype="fasta", retmode="text")
fasta_data = handle.read()
handle.close()

print(fasta_data)

with open("HPV_L1.fasta", "w") as f:
    f.write(fasta_data)





################################################################################################################
########################################### 김유진 2026-04-02 10:01:59 ###########################################
################################################################################################################
#!pip install biopython

from Bio import Entrez
Entrez.email = "yu20230867@kmou.ac.kr"
handle = Entrez.einfo()
record = handle.read()

print(record)

from Bio import Entrez
Entrez.email = "yu20230867@kmou.ac.kr"
handle = Entrez.einfo()
record = Entrez.read(handle)

print(record)
print(len(record["DbList"]))

from Bio import Entrez
Entrez.email = "yu20230867@kmou.ac.kr"
handle = Entrez.esearch(db="pubmed", term="Human Papillomavirus L1", RetMax=5)
record = Entrez.read(handle)

print(record)

print(record["Count"])


from Bio import Entrez
Entrez.email = "yu20230867@g.kmou.ac.kr"
handle = Entrez.efetch(db="protein", id="NP_041332.2", rettype="gb", retmode="text")

print(handle.read())

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "yu20230867@g.kmou.ac.kr"
handle = Entrez.efetch(db="protein", id="NP_041332.2", rettype="gb", retmode="text")
record = SeqIO.parse(handle, "gb")
seq_lists = list(record)
handle.close()
type(seq_lists)

seq_lists[0:4]

seq_rec = seq_lists[0]
fasta_format = ">%s\n%s\n"%(seq_rec.id, seq_rec.seq)

seq_rec

fasta_format

print(fasta_format)




################################################################################################################
########################################### 최승빈 2026-04-02 08:01:13 ###########################################
################################################################################################################
from Bio import Entrez, SeqIO

Entrez.email = "binseung01@gmail.com“

search_term = '"Human papillomavirus"[Organism] AND L1[Gene]’

with Entrez.esearch(db="protein", term=search_term, retmax=5) as handle:
search_results = Entrez.read(handle)

id_list = search_results["IdList"]

print("검색된 ID 5개:")
print(id_list)

with Entrez.efetch(db="protein", id=id_list, rettype="gb", retmode="text") as handle:
gb_records = list(SeqIO.parse(handle, "genbank"))

SeqIO.write(gb_records, "hpv_l1_5records.gb", "genbank")
print("GenBank 파일 저장 완료: hpv_l1_5records.gb")

with Entrez.efetch(db="protein", id=id_list, rettype="fasta", retmode="text") as handle:
fasta_records = list(SeqIO.parse(handle, "fasta"))

SeqIO.write(fasta_records, "hpv_l1_5records.fasta", "fasta")
print("FASTA 파일 저장 완료: hpv_l1_5records.fasta")

print("\n[GenBank records]")
for record in gb_records:
print(record.id, record.description)

print("\n[FASTA records]")
for record in fasta_records:
print(record.id, record.description)




################################################################################################################
########################################### 이지효 2026-04-02 07:00:02 ###########################################
################################################################################################################
from Bio import Entrez 
Entrez.email = "zyo271@gmailcom"

handle = Entrez.esearch( 
  db="protein",
  term="Human papillomavirus L1",
  retmax=5
)
record = Entrez.read(handle) 
handle.close()

print(record)

id_list = record["IdList"] 
print(id_list)

import time

ids = ",".join(id_list) 
time.sleep(1)

handle = Entrez.efetch( 
  db="protein",
  id=ids, 
  rettype="gb", 
  retmode="text"
)

gb_data = handle.read()
handle.close()

print(gb_data)

with open("hpv_L1.gb", "w") as f:
  f.write(gb_data)


time.sleep(1)

handle = Entrez.efetch( 
  db="protein",
  id=ids, 
  rettype="fasta", 
  retmode="text"
)

fasta_data = handle.read()
handle.close()
print(fasta_data)
with open("hpv_L1.fasta", "w") as f: 
  f.write(fasta_data)




################################################################################################################
########################################### 유혜원 2026-04-01 22:01:26 ###########################################
################################################################################################################
# In[ ]:
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "yhw3954@gmail.com"


# In[ ]:
handle = Entrez.esearch(
    db="protein",
    term="Human papillomavirus[Organism] AND L1",
    retmax=5
)
record = Entrez.read(handle)
handle.close()

ids = record["IdList"]

print("=== ID List ===")
print(ids)


# In[ ]:
handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="gb",
    retmode="text"
)
gb_data = handle.read()
handle.close()

print("\n=== GenBank Format ===")
print(gb_data)


# In[ ]:
handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="fasta",
    retmode="text"
)
fasta_data = handle.read()
handle.close()

print("\n=== FASTA Format ===")
print(fasta_data)


# In[ ]:
handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="gb",
    retmode="text"
)
records = SeqIO.parse(handle, "gb")

print("\n=== FASTA (Parsed) ===")
for rec in records:
    print(">%s" % rec.id)
    print(rec.seq)




################################################################################################################
########################################### 이도현 2026-04-01 18:20:40 ###########################################
################################################################################################################
# In[ ]:
pip install biopython

# In[ ]:
from Bio import Entrez
Entrez.email="dlehgus04@gmail.com"

handle=Entrez.esearch(db="protein", term="Human Papillomavirus AND L1", RetMax=5)
record=Entrez.read(handle)

print(record)
print(record["Count"])

# In[ ]:
print(record["IdList"])

# In[ ]:
from Bio import Entrez
Entrez.email="dlehgus04@gmail.com"

ids=record["IdList"]

handle=Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")

# In[ ]:
print(handle.read())

# In[ ]:
from Bio import Entrez
Entrez.email="dlehgus04@gmail.com"

g_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")

# In[ ]:
print(g_handle.read())




################################################################################################################
########################################### 윤이정 2026-04-01 14:36:21 ###########################################
################################################################################################################
#pip install biopython (재시작하기)

from Bio import Entrez, SeqIO

Entrez.email = "ylj16287@gmail.com"

handle = Entrez.esearch(
    db="protein",
    term="Human Papillomavirus L1",
    retmax=5
)
record = Entrez.read(handle)
ids = record["IdList"]

print("IDs:", ids)

handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="gb",
    retmode="text"
)
genbank_data = handle.read()
handle.close()

print("\n=== GenBank format ===\n")
print(genbank_data)

handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="fasta",
    retmode="text"
)
fasta_data = handle.read()
handle.close()

print("\n=== FASTA format ===\n")
print(fasta_data)





################################################################################################################
########################################### 손혁진 2026-03-31 16:49:33 ###########################################
################################################################################################################
#pip install biopython

from Bio import Entrez

from Bio import SeqIO

Entrez.email = "bo3647@daum.net"

handle = Entrez.esearch(db="protein", term = "HPV L1", retmax=5)

record = Entrez.read(handle)

ids = record['IdList']

gb_handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")

gb_record = gb_handle.read()

print(gb_record)



fasta_handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")

fasta_record = gb_handle.read()

print(fasta_record)


################################################################################################################
########################################### 임예지 2026-03-29 23:30:23 ###########################################
################################################################################################################
[1] 
#!pip install biopython

[2]
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "dladpwl5223@naver.com"

[3]
handle = Entrez.einfo()
record = Entrez.read(handle)
handle.close()

print(record["DbList"])
print(len(record["DbList"]))

[4]
handle = Entrez.esearch(
    db="protein",
    term="Human papillomavirus[Organism] AND L1[Gene]",
    retmax=5
)

record = Entrez.read(handle)
handle.close()

print(record)
print(record["IdList"])

[5]
ids = record["IdList"]

handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="gb",
    retmode="text"
)

gb_data = handle.read()
handle.close()

print(gb_data)

[6]
with open("HPV_L1.gb", "w") as f:
    f.write(gb_data)

print("GenBank 저장 완료")

[7]
handle = Entrez.efetch(
    db="protein",
    id=ids,
    rettype="fasta",
    retmode="text"
)

records = SeqIO.parse(handle, "fasta")
seq_lists = list(records)

handle.close()

seq_lists

[8]
for seq_rec in seq_lists:
    fasta_format = ">%s\n%s\n" % (seq_rec.id, seq_rec.seq)
    print(fasta_format)

[9]
with open("HPV_L1.fasta", "w") as f:
    for seq_rec in seq_lists:
        fasta_format = ">%s\n%s\n" % (seq_rec.id, seq_rec.seq)
        f.write(fasta_format)

print("FASTA 저장 완료")


################################################################################################################
########################################### 정서현 2026-03-29 22:45:02 ###########################################
################################################################################################################
from Bio import Entrez
Entrez.email = "blacklu8280@gmail.com"

def fetch_hpv_l1_data():
search_term = "Human Papillomavirus L1 protein"

print(f"Searching for: {search_term}...")

handle = Entrez.esearch(db="protein", term=search_term, retmax=5)
record = Entrez.read(handle)
handle.close()

id_list = record["IdList"]
    
if not id_list:
print("검색 결과가 없습니다.")
return

formats = ["gb", "fasta"]  
for fmt in formats:
print(f"\n--- Fetching in {fmt.upper()} format ---")
fetch_handle = Entrez.efetch(
db="protein", 
id=id_list, 
rettype=fmt, 
retmode="text"
        )
data = fetch_handle.read()
fetch_handle.close()

print(data[:500] + "\n... (중략) ...")

with open(f"hpv_l1_records.{fmt}", "w") as f:
f.write(data)

if __name__ == "__main__":
    fetch_hpv_l1_data()




################################################################################################################
########################################### 김아원 2026-03-28 12:28:52 ###########################################
################################################################################################################
from Bio import Entrez
Entrez.email = "aone0420@gmail.com"

handle = Entrez.esearch(db="protein", term="Human Papillomavirus AND L1", RetMax=5)
record = Entrez.read(handle)
print(record)
record['IdList']

ids = record['IdList']
handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
records = handle.read()
print(records)

################################################################################################################
########################################### 주상민 2026-03-27 11:49:07 ###########################################
################################################################################################################
Entrez.email = "wntkdals0808@gmail.com"

handle = Entrez.esearch(db="protein", term="HPV L1", RetMax=5)
record = Entrez.read(handle)
print(record)
record['IdList']

ids = record['IdList']
handle = Entrez.efetch(db="protein", id=ids, rettype="gb", retmode="text")
records = handle.read()
print(records)



Entrez.email = "wntkdals0808@gmail.com"

handle = Entrez.esearch(db="protein", term="HPV L1", RetMax=5)
record = Entrez.read(handle)
print(record)
record['IdList']

ids = record['IdList']
handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
records = handle.read()
print(records)


################################################################################################################
########################################### 김다진 0/30점 ###########################################
################################################################################################################


################################################################################################################
########################################### 전주희 0/30점 ###########################################
################################################################################################################

