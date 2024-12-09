import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.linear_model import Ridge
from scipy.stats import pearsonr


#pass each parameter seperately for the actual function
def run_ridge_reg(
        adata_rna,
        adata_metabolomics, 
        params, 
        **kwargs):
    # Extract data matrices
    
    split = adata_rna.obs['split']
    X_rna = adata_rna.X
    Y_metabolomics = adata_metabolomics.X

    # Train-test split based on 'split' column
    train_idx = np.where(split == 'train')[0]
    test_idx = np.where(split == 'test')[0]

    X_train, X_test = X_rna[train_idx], X_rna[test_idx]
    Y_train, Y_test = Y_metabolomics[train_idx], Y_metabolomics[test_idx]

    # Ridge regression with specified alpha
    alpha = float(params['alpha'] ) # have as function parameter
    ridge = Ridge(alpha=alpha)
    ridge.fit(X_train, Y_train)

    # Predictions and evaluation
    Y_pred = ridge.predict(X_test)
    pearson_corr = pearsonr(Y_pred.flatten(), Y_test.flatten())[0]

    # Save results to a DataFrame
    results = pd.DataFrame({
    'corr': [pearson_corr], 
    'alpha': [alpha]})

    return results