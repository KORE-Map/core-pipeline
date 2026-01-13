# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 09:07:25 2026

@author: KIOM_User
"""

import os

os.chdir(r"your directory")

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

xlsx = "figure2_table.xlsx"
sheet_names = ["A549", "HepG2", "HT29", "SW1783"]

df_list = []
for s in sheet_names:
    df = pd.read_excel(xlsx, sheet_name=s)
    df["Cell"] = s
    df_list.append(df)


data = pd.concat(df_list, ignore_index=True)
data.columns = ["Sample_ID", "A260_A280", "rRNA_28S_18S", "RIN", "Cell"]

# palette
palette = ["#9A3DFF", "#00C7FF", "#007C7A", "#FFD300"]
plot_params = dict(palette=palette, linewidth=0.8, edgecolor="black", size=6)


# A260/A280
plt.figure(figsize=(4, 4), dpi=600)
sns.boxplot(data=data, x="Cell", y="A260_A280", palette=palette,
            width=0.5, fliersize=0, linewidth=1.2, boxprops=dict(alpha=0.4))
sns.stripplot(data=data, x="Cell", y="A260_A280", jitter=True, **plot_params)
plt.ylim(1.0, 2.02)
plt.yticks(np.arange(1.0, 2.02, 0.25), fontsize=6)

for y in [1.5, 2.0]:
    plt.axhline(y=y, color='gray', linestyle=':', linewidth=1)

for spine in ['top', 'right']:
    plt.gca().spines[spine].set_visible(False)

plt.tick_params(axis='x', which='both', top=False, length=0)
plt.xticks(fontsize=10)
plt.xlabel(None)
plt.ylabel("A260/A280 ratio", fontsize=8, weight='bold')
plt.title("A260/A280 ratio", fontsize=13, pad=10)
plt.tight_layout()


# 28S/18S
plt.figure(figsize=(4, 4), dpi=600)
sns.boxplot(data=data, x="Cell", y="rRNA_28S_18S", palette=palette,
            width=0.5, fliersize=0, linewidth=1.2, boxprops=dict(alpha=0.4))
sns.stripplot(data=data, x="Cell", y="rRNA_28S_18S", jitter=True, **plot_params)

plt.ylim(-0.1, 4.1)
plt.yticks(np.arange(0, 4.1, 1), fontsize=6)

plt.ylabel("rRNA 28S/18S ratio", fontsize=8, weight='bold')
plt.axhline(y=2, color='gray', linestyle=':', linewidth=1)

for spine in ['top', 'right']:
    plt.gca().spines[spine].set_visible(False)

plt.tick_params(axis='x', which='both', top=False, length=0)
plt.xticks(fontsize=10)
plt.xlabel(None)
plt.title("rRNA 28S/18S ratio", fontsize=13, pad=10)
plt.tight_layout()


# RIN
plt.figure(figsize=(4, 4), dpi=600)
sns.boxplot(data=data, x="Cell", y="RIN", palette=palette,
            width=0.5, fliersize=0, linewidth=1.2, boxprops=dict(alpha=0.4))
sns.stripplot(data=data, x="Cell", y="RIN", jitter=True, **plot_params)
plt.ylim(4.9, 10.1)
plt.yticks(np.arange(5, 10.1, 1), fontsize=6)

for spine in ['top', 'right']:
    plt.gca().spines[spine].set_visible(False)
    
plt.ylabel("RNA Integrity Number (RIN)", fontsize=8, weight='bold')
plt.axhline(y=7, color='gray', linestyle=':', linewidth=1)
plt.tick_params(axis='x', which='both', top=False, length=0)
plt.xticks(fontsize=10)
plt.xlabel(None)
plt.title("RNA integrity number", fontsize=13, pad=10)
plt.tight_layout()

