#!/usr/bin/env python
# coding: utf-8

# Import the modules needed

# In[334]:


import pandas as pd
import scipy as sp
from IPython.display import HTML, display
import scipy.stats as stats
import itertools

# set the column display maxium to reduce clutter
pd.options.display.max_rows = 5
pd.options.display.max_columns = 10

# for progress bar
from ipywidgets import IntProgress
from IPython.display import display
import time

# for download
def create_download_link(title = "Download CSV file", filename = "pairwise_genes.csv"):  
    html = '<a download="{filename}" href="{filename}" target="_blank">{title}</a>'
    html = html.format(title=title,filename=filename)
    return HTML(html)


# Open the excel file and import the data into a dataframe. _Note: you can just add more data as long as it's in the same format, and it will rebuild everything._

# In[335]:


df = pd.read_excel("combined-correlation.xlsx")
df


# Like I was talking about, in programming, often the columns don't represent data in the same way we use them in excel. Instead of having a second dimension as a column, we just add all that data as more rows. Pandas even has a built-in function to do this called ``melt``. First, we need to make a column for the patient id.

# In[336]:


df["patient_id"] = [*range(1,125)]
df.set_index("patient_id")


# We can then pivot thedata from column to a per-row property while maintaining the `patient_id` while using `melt`. We will assign that to a different dataframe called mutations dataframe, or `mdf` for short

# In[337]:


mdf = df.melt(id_vars=["patient_id"],var_name="gene",value_name="is_mutated")
(gene,is_mutated) =  (mdf['gene'], mdf['is_mutated'])
mdf


# 

# In[ ]:





# Now, we will make a new dataframe with each gene-gene combination as a row. 
# 
# First, let's get all the names of the genes and put that in a list.

# In[338]:


gene_names = df2.gene.unique()
# this is a list of two list
gene_matrix = [gene_names, gene_names]
# this is all the possible combinations of gene 1 and gene 2
gene_product = list(itertools.product(*gene_matrix))
# make that a dataframe
gene_df = pd.DataFrame(gene_product, columns=["gene1","gene2"])
gene_df


# Now we can add those other counting rows with blank default values

# In[339]:


gene_df["Both mutated"] = ""
gene_df["Gene 1 Mutated"] = ""
gene_df["Gene 2 Mutated"] = ""
gene_df["Niether mutated"] = ""
gene_df["Odds ratio"] = ""
gene_df["P value"] = ""
gene_df


# Now for the hard part. We'll walk through each row of `gene_df` and add up how many instances of each of those examples there are.

# In[340]:


# an example, this shows all the patients who have ASXL1 mutated

asx = mdf[(is_mutated) & (gene == "ASXL1")].patient_id.to_list()
asx


# In[341]:


# This shows all the patients with TET2 mutated
tet = mdf[(is_mutated) & (gene == "TET2")].patient_id.to_list()
tet


# In[342]:


# only one patient, number 41, has both
common = list(set(tet).intersection(asx))
common


# In[344]:


# progress bar
max_count = len(gene_df)
f = IntProgress(min=0, max=max_count) # instantiate the bar
display(f) # display the bar
count = 0

# now, or each row in that dataframe, do the following code, 
# where index is the row number, and row is an array/list 
# with all the value from that row

num_total_pt = len(mdf.patient_id.unique())

for index, row in gene_df.iterrows():
    pt_1_mutated = mdf[(is_mutated) & (gene == row["gene1"])].patient_id.to_list()
    pt_2_mutated = mdf[(is_mutated) & (gene == row["gene2"])].patient_id.to_list()  
    pt_both_mutated = list(set(pt_1_mutated).intersection(pt_2_mutated))
    # now get sizes
    num_pt_both_mutated = len(pt_both_mutated)
    num_pt_only_1_mutated = len(pt_1_mutated) - num_pt_both_mutated
    num_pt_only_2_mutated = len(pt_2_mutated) - num_pt_both_mutated
    num_pt_neither = num_total_pt - num_pt_both_mutated - num_pt_only_1_mutated - num_pt_only_2_mutated
    
    # store these in dataframe  
    gene_df.loc[index, "Both mutated"] = num_pt_both_mutated
    gene_df.loc[index, "Gene 1 Mutated"] = num_pt_only_1_mutated
    gene_df.loc[index, "Gene 2 Mutated"] = num_pt_only_2_mutated
    gene_df.loc[index, "Niether mutated"] = num_pt_neither
    
    # build fisher table
    fisher_table = [[num_pt_both_mutated, num_pt_only_1_mutated],[num_pt_only_2_mutated, num_pt_neither]]
    (odds_radio,pvalue) = stats.fisher_exact(fisher_table)
    gene_df.loc[index, "Odds ratio"] = odds_radio
    gene_df.loc[index, "P value"] = pvalue
    f.value += 1 # signal to increment the progress bar

gene_df


# Can download as CSV below

# In[345]:


gene_df.to_csv("pairwise_genes.csv")
create_download_link(filename="pairwise_genes.csv")

