import time
# start_time = time.time()

import scanpy as sc # type: ignore
import pandas as pd # type: ignore
import ast


# Import all method scripts
from ridge import run_ridge_reg

# Dictionary mapping method names to their functions - connect config->smk->run_methods->function.py 
METHOD_MAP = {
    'ridge': dict(function=run_ridge_reg, mode='paired')
}

# Load parameters from Snakemake
params = snakemake.params.thisparam 
 # Task parameters
method = params['method']
task = params['task']
hash_id = params['hash']
input_rna = params['input_rna']
input_metabolomics = params['input_metabolomics']
method_params = ast.literal_eval(params['params'])

output_file = snakemake.output.tsv

# Load the appropriate method function
method_mode = METHOD_MAP[method]['mode']
method_function = METHOD_MAP[method]['function']

# Load data based on method mode
if method_mode == 'paired':
    adata_rna = sc.read_h5ad(input_rna)
    adata_metabolomics = sc.read_h5ad(input_metabolomics)
    result_df = method_function(
        adata_rna=adata_rna,
        adata_metabolomics=adata_metabolomics,
        params=method_params
        )

# Add metadata to the results and save to output file
result_df['hash'] = hash_id
result_df['task'] = task
result_df['method_params'] = method_params
result_df.to_csv(output_file, sep='\t', index=False)
