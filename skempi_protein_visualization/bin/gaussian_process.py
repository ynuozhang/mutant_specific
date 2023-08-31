import torch
import esm
import pandas as pd
import numpy as np
from fuzzywuzzy import fuzz
from umap import UMAP
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from sklearn.gaussian_process.kernels import ConstantKernel as C
from sklearn.model_selection import train_test_split

def add_scaled_columns(dataframe, col1, col2, boundary=10000.):
    # Extract the two columns
    sub_df = dataframe[[col1, col2]].copy()

    # Create a combined series for calculating scaling factor
    combined = sub_df[col1].append(sub_df[col2])

    # 1. Handle NaNs
    combined_filled = combined.fillna(1e9)

    # 2. Determine Scaling Factor
    # Find the maximum value excluding the NaN placeholders
    max_val = combined_filled[combined_filled != 1e9].max()

    # Compute the scaling factor so the max_val is just below 10000
    scaling_factor = (boundary - 1) / max_val

    # 3. Normalize Data
    # Scale the two columns and add them to the original dataframe
    dataframe[f"scaled_{col1}"] = sub_df[col1].mul(scaling_factor).fillna(boundary)
    dataframe[f"scaled_{col2}"] = sub_df[col2].mul(scaling_factor).fillna(boundary)
    return dataframe

# concatenate wt1+wt2 and wt1+mut
def concatenate_vectors(row):
    return np.concatenate([row['wild_seq_1_embeddings'], row['wild_seq_2_embeddings']], axis=1)

# concatenate wt1+wt2 and wt1+mut
def concatenate_vectors_mut(row):
    return np.concatenate([row['wild_seq_1_embeddings'], row['mutant_seq_embeddings']], axis=1)

df = pd.read_hdf('reduced_proteins_embeddings_meta.hdf', key='df')
df['wt1_wt2_concat'] = df.apply(concatenate_vectors, axis=1)
df['wt1_mut_concat'] = df.apply(concatenate_vectors_mut, axis=1)

# alternative method of normalization
df['Affinity_mut_parsed'] = df['Affinity_mut_parsed'] * (10 ** 9)
df['Affinity_wt_parsed'] = df['Affinity_wt_parsed'] * (10 ** 9)
df = add_scaled_columns(df, 'Affinity_mut_parsed', 'Affinity_wt_parsed')

max_molar = df[['scaled_Affinity_mut_parsed', 'scaled_Affinity_wt_parsed']].max().max()
dt = df.copy()
dt.drop_duplicates(subset=['wild_seq_1', 'wild_seq_2'], inplace=True)
