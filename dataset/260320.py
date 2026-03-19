# In[ ]:


### 2026-03-20 ###

# Scopus: https://www.scopus.com/
# pybliometrics: Python library to pull, cache and extract data from the Scopus database.
#                https://pybliometrics.readthedocs.io/en/stable/index.html
# pip install pybliometrics

import pybliometrics


# In[ ]:


# Scopus API generation at https://dev.elsevier.com/
# Use my API: ???

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


# How many papers?
s
s.get_results_size()


# In[ ]:


# Make database from the paper info

import pandas as pd
df_s
df_s = pd.DataFrame(s.results)


# In[ ]:


# DB check

df_s.head()


# In[ ]:


# DB details
s0 = df_s.loc[37]
s0


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




