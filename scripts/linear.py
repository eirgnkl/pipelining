import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import spearmanr, pearsonr
import numpy as np
from scipy.sparse import issparse



def convert_to_dense(matrix):
    """Converts a sparse matrix to dense if necessary."""
    if issparse(matrix):
        return matrix.toarray()
    return matrix

def run_linreg(
        adata_rna,
        adata_metabolomics, 
        **kwargs):
    # Extract data matrices
    
    split = adata_rna.obs['split']
    X_rna = adata_rna.X
    Y_metabolomics = adata_metabolomics.X

    #convert if needed
    X_rna = convert_to_dense(X_rna)
    Y_metabolomics = convert_to_dense(Y_metabolomics)

    # Train-test split based on 'split' column
    train_idx = np.where(split == 'train')[0]
    test_idx = np.where(split == 'test')[0]

    X_train, X_test = X_rna[train_idx], X_rna[test_idx]
    Y_train, Y_test = Y_metabolomics[train_idx], Y_metabolomics[test_idx]

    # Linear regression
    lin = LinearRegression()
    lin.fit(X_train, Y_train)

    # Predictions and evaluation
    Y_pred = lin.predict(X_test)
    pearson_corr = pearsonr(Y_pred.flatten(), Y_test.flatten())[0]

    # Save results to a DataFrame
    results = pd.DataFrame({
    'corr': [pearson_corr]})

    return results