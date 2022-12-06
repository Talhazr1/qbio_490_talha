###constructings umap plots
import os

print("Current working directory: {0}".format(os.getcwd()))

os.chdir('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data')

print("New working directory: {0}".format(os.getcwd()))

import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns
import cptac as cp
import sklearn.cluster as cluster

import umap
import numpy as np
import pandas as pd
%matplotlib inline 
from sklearn.cluster import KMeans
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

# download the GBM dataset
cp.download(dataset="Gbm")
gbm = cp.Gbm()

# extract the data
clinical_data = gbm.get_clinical()
transcriptomic_data = gbm.get_transcriptomics()
protein_data = gbm.get_proteomics()
protein_data.columns = protein_data.columns.get_level_values(0) 
protein_data = protein_data.dropna(axis=1)

wcss = []
for i in range(1, 11):
    kmeans = KMeans(n_clusters=i, init='k-means++', max_iter=300, n_init=10, random_state=0)
    kmeans.fit(protein_data)
    wcss.append(kmeans.inertia_)
plt.plot(range(1, 11), wcss)
plt.title('Fig.1 Elbow Method')
plt.xlabel('Number of clusters')
plt.ylabel('WCSS')
plt.show()
#plateu begins to occur around 4

gbm_rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)

import sklearn.cluster as cluster
import umap
import random
rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)

from sklearn.preprocessing import OneHotEncoder
onehotencoder = OneHotEncoder(handle_unknown='ignore')
rna_clinical = onehotencoder.fit_transform(rna_clinical).toarray()

mapper = umap.UMAP().fit_transform(rna_clinical)

gene_labels = rna_clinical
plt.scatter(mapper[:, 0], mapper[:, 1], s=10, cmap="rainbow");
plt.title("Figure 2. UMAP of a single gene's expression")
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)
import sklearn.cluster as cluster
import umap
import random
from sklearn.preprocessing import OneHotEncoder
#onehotencoder = OneHotEncoder(handle_unknown='ignore')
#rna_clinical = onehotencoder.fit_transform(rna_clinical).toarray()

mapper = umap.UMAP().fit_transform(RNA_data)

#kmeans_labels = cluster.KMeans(n_clusters=4).fit_predict(rna_clinical)
plt.scatter(mapper[:, 0], mapper[:, 1], c=color_list,s=10);
plt.title("Figure 3. UMAP of RNA Expression by Patient Sex ")
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.show()

gbm_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/hcc_clinical.csv',index_col=0)
rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)

RNA_data = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_counts.csv',index_col=0)
RNA_data = RNA_data.T


RNA_data.index = rna_clinical.patient


RNA_data = RNA_data[~RNA_data.index.duplicated(keep='first')]

patients = np.intersect1d(gbm_clinical.Tumor_Sample_Barcode, RNA_data.index)

colors= []
for i in RNA_data.index:
    colors.append(gbm_clinical.loc[gbm_clinical.loc[:,'Tumor_Sample_Barcode'] == i, 'gender'])

test = []
for i in colors:
    try:
        test.append(i.values[0])
    except:
        print()

test = pd.Series(test)
    
RNA_data = RNA_data.loc[RNA_data.index.isin(patients),:]
color_list = np.where(test=='MALE', 'blue', 'red')
RNA_data = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_counts.csv',index_col=0)
RNA_data = RNA_data.T


RNA_data.index = rna_clinical.patient


RNA_data = RNA_data[~RNA_data.index.duplicated(keep='first')]

patients = np.intersect1d(gbm_clinical.Tumor_Sample_Barcode, RNA_data.index)

colors= []
for i in RNA_data.index:
    colors.append(gbm_clinical.loc[gbm_clinical.loc[:,'Tumor_Sample_Barcode'] == i, 'race_list'])

test = []
for i in colors:
    try:
        test.append(i.values[0])
    except:
        print()

test = pd.Series(test)
    
RNA_data = RNA_data.loc[RNA_data.index.isin(patients),:]
color_list = np.where(test=='ASIAN', 'blue', (np.where(test=="WHITE", 'red','yellow')) )
rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)
import sklearn.cluster as cluster
import umap
import random
from sklearn.preprocessing import OneHotEncoder
#onehotencoder = OneHotEncoder(handle_unknown='ignore')
#rna_clinical = onehotencoder.fit_transform(rna_clinical).toarray()

mapper = umap.UMAP().fit_transform(RNA_data)

#kmeans_labels = cluster.KMeans(n_clusters=4).fit_predict(rna_clinical)
plt.scatter(mapper[:, 0], mapper[:, 1], c=color_list,s=10);
plt.title("Figure 3. UMAP of RNA Expression by Patient Race ")
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.show()
RNA_data = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_counts.csv',index_col=0)
RNA_data = RNA_data.T


RNA_data.index = rna_clinical.patient


RNA_data = RNA_data[~RNA_data.index.duplicated(keep='first')]

patients = np.intersect1d(gbm_clinical.Tumor_Sample_Barcode, RNA_data.index)

colors= []
for i in RNA_data.index:
    colors.append(gbm_clinical.loc[gbm_clinical.loc[:,'Tumor_Sample_Barcode'] == i, 'age_class'])

test = []
for i in colors:
    try:
        test.append(i.values[0])
    except:
        print()

test = pd.Series(test)
    
RNA_data = RNA_data.loc[RNA_data.index.isin(patients),:]
color_list = np.where(test=='Young', 'blue', 'red' )
rna_clinical = pd.read_csv('/Users/talhaaccount/Desktop/Qbio_490_R/qbio_490_talha/analysis_data/gbm_rna_clinical_data.csv',index_col=0)
import sklearn.cluster as cluster
import umap
import random
from sklearn.preprocessing import OneHotEncoder
#onehotencoder = OneHotEncoder(handle_unknown='ignore')
#rna_clinical = onehotencoder.fit_transform(rna_clinical).toarray()

mapper = umap.UMAP().fit_transform(RNA_data)

#kmeans_labels = cluster.KMeans(n_clusters=4).fit_predict(rna_clinical)
plt.scatter(mapper[:, 0], mapper[:, 1], c=color_list,s=10);
plt.title("Figure 3. UMAP of RNA Expression by Patient Age ")
plt.xlabel('UMAP 1')
plt.ylabel('UMAP 2')
plt.show()

###constructing heatmaps
# download the Gbm dataset
cptac.download(dataset="Gbm")
gbm = cptac.Gbm()

# extract the data
clinical_data = gbm.get_clinical()
transcriptomic_data = gbm.get_transcriptomics()
protein_data = gbm.get_proteomics()
protein_data.columns = protein_data.columns.get_level_values(0) 
protein_data.columns

from scipy import stats # we are using the stats package in particular

# 1. Identify the genes (RNA, protein) shared between the two data sets 
shared_rna_prot = np.intersect1d(transcriptomic_data.columns.get_level_values(0), protein_data.columns)

# 2. Create the two data frames
rna_shared = transcriptomic_data[shared_rna_prot]
prot_shared = protein_data[shared_rna_prot]
prot_shared.drop(prot_shared.tail(1).index,inplace=True)

ncomparisons = 10 # define this variable in case we want to change the number of correlations to test
                  # this makes it less likely you'll forget to change a number, e.g. in the data frame shape
gene_names = ["IDH1", "TP53", "ATRX", "PIK3CA", "PIK3C2B", "PIK3R1", "NF1", "EGFR", "PTEN", "MAPRE3"]


# Don't worry about this code
# It's good practice to declare your data frame beforehand (it's much faster than appending to a list)
# We fill everything in with 0 just as a placeholder
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)


# 2. fill in the data frame!
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        # then, use .loc[] to store the correlation in corr_df
        corr, pval = stats.spearmanr(rna_shared[g1], prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1, g2] = corr

# 3. create the heat map
plot = sns.heatmap(
    corr_df,  
    cmap='mako',
)
plot.set_xlabel('Protein', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)

# 4. interpret!
#There is a fairly evenly spread distribution of correlations, ranging from
#strong negative correaltions denoted by black squares to strong positive
#correlations denoted by light green scores (predominant across diagonal 
#since RNA-Protein expression levels are likely to be similar for the same
#gene.

young_clinical_data = clinical_data[clinical_data['age'] <= 50]

patients = np.intersect1d(young_clinical_data.index, transcriptomic_data.index) 
patients = np.intersect1d(patients, protein_data.index)

young_transcriptomic_data = transcriptomic_data.loc[patients, :] #filter transcriptomic data by matching to Patient ID's from clinical 
young_proteomic_data = protein_data.loc[patients, :] #match with patient iDs from clinical
# 1. Identify the genes (RNA, protein) shared between the two data sets 
young_shared_rna_prot = np.intersect1d(young_transcriptomic_data.columns.get_level_values(0), young_proteomic_data.columns)

# 2. Create the two data frames
young_rna_shared = young_transcriptomic_data[young_shared_rna_prot]
young_prot_shared = young_proteomic_data[young_shared_rna_prot]

young_shared_rna_prot

ncomparisons = 10 # define this variable in case we want to change the number of correlations to test
                  # this makes it less likely you'll forget to change a number, e.g. in the data frame shape
gene_names = ["IDH1", "TP53", "ATRX", "PIK3CA", "PIK3C2B", "PIK3R1", "NF1", "EGFR", "PTEN", "MAPRE3"]


# Don't worry about this code
# It's good practice to declare your data frame beforehand (it's much faster than appending to a list)
# We fill everything in with 0 just as a placeholder
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)


# 2. fill in the data frame!
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        # then, use .loc[] to store the correlation in corr_df
        corr, pval = stats.spearmanr(young_rna_shared[g1], young_prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1, g2] = corr

# 3. create the heat map
plot = sns.heatmap(
    corr_df,  
    cmap='mako',
)
plot.set_xlabel('Protein (Young Patient)', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)

# 4. interpret!
#There is a fairly evenly spread distribution of correlations, ranging from
#strong negative correaltions denoted by black squares to strong positive
#correlations denoted by light green scores (predominant across diagonal 
#since RNA-Protein expression levels are likely to be similar for the same
#gene.

old_clinical_data = clinical_data[clinical_data['age'] > 50]

patients = np.intersect1d(old_clinical_data.index, transcriptomic_data.index) 
patients = np.intersect1d(patients, protein_data.index)

old_transcriptomic_data = transcriptomic_data.loc[patients, :] #filter transcriptomic data by matching to Patient ID's from clinical 
old_proteomic_data = protein_data.loc[patients, :] #match with patient iDs from clinical

# 1. Identify the genes (RNA, protein) shared between the two data sets 
old_shared_rna_prot = np.intersect1d(old_transcriptomic_data.columns.get_level_values(0), old_proteomic_data.columns)

# 2. Create the two data frames
old_rna_shared = old_transcriptomic_data[old_shared_rna_prot]
old_prot_shared = old_proteomic_data[old_shared_rna_prot]

old_shared_rna_prot

ncomparisons = 10 # define this variable in case we want to change the number of correlations to test
                  # this makes it less likely you'll forget to change a number, e.g. in the data frame shape
gene_names = ["IDH1", "TP53", "ATRX", "PIK3CA", "PIK3C2B", "PIK3R1", "NF1", "EGFR", "PTEN", "MAPRE3"]


# Don't worry about this code
# It's good practice to declare your data frame beforehand (it's much faster than appending to a list)
# We fill everything in with 0 just as a placeholder
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)


# 2. fill in the data frame!
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        # then, use .loc[] to store the correlation in corr_df
        corr, pval = stats.spearmanr(old_rna_shared[g1], old_prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1, g2] = corr

# 3. create the heat map
plot = sns.heatmap(
    corr_df,  
    cmap='mako',
)
plot.set_xlabel('Protein (old Patient)', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)

# 4. interpret!
#There is a fairly evenly spread distribution of correlations, ranging from
#strong negative correaltions denoted by black squares to strong positive
#correlations denoted by light green scores (predominant across diagonal 
#since RNA-Protein expression levels are likely to be similar for the same
#gene.

White_clinical_data = clinical_data[clinical_data['ethnicity_self_identify'] == 'Caucasian']

patients = np.intersect1d(White_clinical_data.index, transcriptomic_data.index) 
patients = np.intersect1d(patients, protein_data.index)

White_transcriptomic_data = transcriptomic_data.loc[patients, :] #filter transcriptomic data by matching to Patient ID's from clinical 
White_proteomic_data = protein_data.loc[patients, :] #match with patient iDs from clinical

# 1. Identify the genes (RNA, protein) shared between the two data sets 
White_shared_rna_prot = np.intersect1d(White_transcriptomic_data.columns.get_level_values(0), White_proteomic_data.columns)

# 2. Create the two data frames
White_rna_shared = White_transcriptomic_data[White_shared_rna_prot]
White_prot_shared = White_proteomic_data[White_shared_rna_prot]
ncomparisons = 10 # define this variable in case we want to change the number of correlations to test
                  # this makes it less likely you'll forget to change a number, e.g. in the data frame shape
gene_names = ["IDH1", "TP53", "ATRX", "PIK3CA", "PIK3C2B", "PIK3R1", "NF1", "EGFR", "PTEN", "MAPRE3"]


# Don't worry about this code
# It's good practice to declare your data frame beforehand (it's much faster than appending to a list)
# We fill everything in with 0 just as a placeholder
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)


# 2. fill in the data frame!
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        # then, use .loc[] to store the correlation in corr_df
        corr, pval = stats.spearmanr(White_rna_shared[g1], White_prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1, g2] = corr

# 3. create the heat map
plot = sns.heatmap(
    corr_df,  
    cmap='mako',
)
plot.set_xlabel('Protein' '             (Race = White)', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)

# 4. interpret!
#There is a fairly evenly spread distribution of correlations, ranging from
#strong negative correaltions denoted by black squares to strong positive
#correlations denoted by light green scores (predominant across diagonal 
#since RNA-Protein expression levels are likely to be similar for the same
#gene.

asian_clinical_data = clinical_data[clinical_data['ethnicity_self_identify'] == 'Han nationality']

patients = np.intersect1d(asian_clinical_data.index, transcriptomic_data.index) 
patients = np.intersect1d(patients, protein_data.index)

asian_transcriptomic_data = transcriptomic_data.loc[patients, :] #filter transcriptomic data by matching to Patient ID's from clinical 
asian_proteomic_data = protein_data.loc[patients, :] #match with patient iDs from clinical
# 1. Identify the genes (RNA, protein) shared between the two data sets 
asian_shared_rna_prot = np.intersect1d(asian_transcriptomic_data.columns.get_level_values(0), asian_proteomic_data.columns)

# 2. Create the two data frames
asian_rna_shared = asian_transcriptomic_data[asian_shared_rna_prot]
asian_prot_shared = asian_proteomic_data[asian_shared_rna_prot]

asian_shared_rna_prot

ncomparisons = 10 # define this variable in case we want to change the number of correlations to test
                  # this makes it less likely you'll forget to change a number, e.g. in the data frame shape
gene_names = ["IDH1", "TP53", "ATRX", "PIK3CA", "PIK3C2B", "PIK3R1", "NF1", "EGFR", "PTEN", "MAPRE3"]


# Don't worry about this code
# It's good practice to declare your data frame beforehand (it's much faster than appending to a list)
# We fill everything in with 0 just as a placeholder
corr_df = pd.DataFrame(np.ndarray(shape=(ncomparisons, ncomparisons), dtype=np.float16),
                      index = gene_names,
                      columns = gene_names)


# 2. fill in the data frame!
for g1 in gene_names:
    for g2 in gene_names:
        # calculate the correlations between protein and RNA
        # then, use .loc[] to store the correlation in corr_df
        corr, pval = stats.spearmanr(asian_rna_shared[g1], asian_prot_shared[g2], nan_policy="omit")
        corr_df.loc[g1, g2] = corr

# 3. create the heat map
plot = sns.heatmap(
    corr_df,  
    cmap='mako',
)
plot.set_xlabel('Protein' '             (Race = Asian)', fontsize=10)
plot.set_ylabel('RNA', fontsize=10)

# 4. interpret!
#There is a fairly evenly spread distribution of correlations, ranging from
#strong negative correaltions denoted by black squares to strong positive
#correlations denoted by light green scores (predominant across diagonal 
#since RNA-Protein expression levels are likely to be similar for the same
#gene.