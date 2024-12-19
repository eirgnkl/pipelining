import pandas as pd
from snakemake.utils import Paramspace
from scripts.utils import create_tasks_df
from pprint import pprint
import numpy as np

# Generate tasks DataFrame and load configuration
# configfile: 'config.yaml'
tasks_df = create_tasks_df('config.yaml', save='data/tasks.tsv')
tasks_df = pd.read_csv('data/tasks.tsv', sep='\t')

# Extract unique task details
hashes = tasks_df['hash'].unique()
methods = tasks_df['method'].unique()
tasks = tasks_df['task'].unique()

# rule all:
#     input:
#         expand('data/reports/{task}/{method}/{hash}/accuracy.tsv', 
#                task=tasks_df['task'].unique(), 
#                method=tasks_df['method'].unique(), 
#                hash=tasks_df['hash'].unique()),
#         'data/reports/best_results.tsv'        

rule all:
    input:
        [
            expand(
                'data/reports/{task}/{method}/{hash}/accuracy.tsv',
                task=tasks_df['task'].unique(),
                method=[method],
                hash=tasks_df[tasks_df['method'] == method]['hash'].unique()
            )
            for method in tasks_df['method'].unique()
        ],
        'data/reports/best_results.tsv'
        
#add new rule for dataset preprocessing
#1. add feature selection as a column in the create tasks df
#2. identify unique modes
#3. prepare and store datasets accordingly
#4. add correct input ds on the downstream tasks


rule run_method:
    output:
        tsv='data/reports/{task}/{method}/{hash}/accuracy.tsv'
    params:
        thisparam=lambda wildcards: tasks_df.loc[tasks_df['hash'] == wildcards.hash, :].iloc[0, :].to_dict()
    script:
        'scripts/run_method.py'


rule merge:
    input:
        expand(
            'data/reports/{task}/{method}/{hash}/accuracy.tsv',
            zip,
            task=tasks_df['task'].values,
            method=tasks_df['method'].values,
            hash=tasks_df['hash'].values
        )
    output:
        tsv='data/reports/merged_results.tsv'
    run:
        dfs = [pd.read_csv(file, sep='\t') for file in input if os.path.exists(file)]
        merged_df = pd.concat(dfs)
        merged_df.to_csv(output.tsv, sep='\t', index=False)

rule find_best:
    input:
        tsv='data/reports/merged_results.tsv'
    output:
        tsv='data/reports/best_results.tsv'
    run:
        df = pd.read_csv(input.tsv, sep='\t')
        best_row = df.loc[df['corr'].idxmax()]
        best_row.to_frame().T.to_csv(output.tsv, sep='\t', index=False)

