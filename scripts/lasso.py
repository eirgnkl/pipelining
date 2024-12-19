import scanpy as sc
import pandas as pd
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr, pearsonr
import numpy as np

def run_lasso(
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

    # Lasso regression with specified alpha
    alpha = float(params['alpha'] ) # have as function parameter
    lassi = Lasso(alpha=alpha)
    lassi.fit(X_train, Y_train)

    # Predictions and evaluation
    Y_pred = lassi.predict(X_test)
    pearson_corr = pearsonr(Y_pred.flatten(), Y_test.flatten())[0]

    # Save results to a DataFrame
    results = pd.DataFrame({
    'corr': [pearson_corr], 
    'alpha': [alpha]})

    return results