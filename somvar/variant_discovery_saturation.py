from cyvcf2 import VCF
import pandas as pd
import numpy as np
from typing import List, Set, Tuple, Optional
from scipy.optimize import curve_fit


def michaelis_menten(x, Vmax, K):
    """
    Michaelis-Menten model function for variant discovery saturation.
    
    Args:
        x: Number of samples.
        Vmax: Maximum number of discoverable variants (asymptote).
        K: Half-saturation constant (samples needed to reach half of Vmax).
        
    Returns:
        Expected cumulative variant count.
    """
    return (Vmax * x) / (K + x)

def neg_exponential(x, a, b):
    """
    Negative exponential model for variant discovery saturation.
    
    Args:
        x: Number of samples.
        a: Asymptotic number of variants.
        b: Discovery rate constant.
        
    Returns:
        Expected cumulative variant count.
    """
    return a * (1 - np.exp(-b * x))

def hybrid_saturation(x, a, b, c):
    """
    Hybrid variant discovery saturation model.

    Models the cumulative number of unique variants discovered as a function of
    the number of individuals sampled. It combines an exponential saturation 
    term (for shared variants) with a linear term (for singletons or rare variants).

    Parameters
    ----------
    x : array-like or float
        Number of individuals sampled (can be scalar or array).
    a : float
        Asymptotic number of shared (non-singleton) variants.
    b : float
        Discovery rate constant for shared variants. Higher values saturate faster.
    c : float
        Average number of singleton (or private) variants discovered per sample.

    Returns
    -------
    array-like or float
        Estimated cumulative number of unique variants discovered.

    Formula
    -------
    V(n) = a * (1 - exp(-b * n)) + c * n
    """
    return a * (1 - np.exp(-b * x)) + c * x
    
def hybrid_power(x, a, b, c, alpha):
    return a * (1 - np.exp(-b * x)) + c * (x ** alpha)

def load_gt_matrix(vcf_path, threads=1, chr='all'):
    """
    Load genotype data from a VCF file and return a binary presence/absence matrix.

    Genotypes are encoded using cyvcf2's `gt_types` with `gts012=True`:
        - 0: Homozygous reference (0/0)
        - 1: Heterozygous (0/1 or 1/0)
        - 2: Homozygous alternate (1/1)
        - 3: Missing (./.)

    The returned matrix is transformed as follows:
        - HET (1) and HOM_ALT (2) → 1 (variant present)
        - HOM_REF (0) and MISSING (3) → 0 (variant absent or missing)

    Args:
        vcf_path (str): Path to the input VCF file.
        threads (int): Number of threads to use when reading the VCF.
        chr (str): Chromosome name to filter on, or 'all' to include all chromosomes.

    Returns:
        np.ndarray: A 2D numpy array of shape (variants, samples) with binary values:
                    1 = variant present (HET or HOM_ALT),
                    0 = absent or missing.
    """
    vcf = VCF(vcf_path, threads=20,gts012=True)
    gt_all = []
    chrname=chr
    if not chrname=='all':
        for v in vcf:
            #print(v.CHROM)
            if v.CHROM==chrname:
                gt_all.append(list(v.gt_types))

    else:
        for v in vcf:
            gt_all.append(list(v.gt_types))
    gtarr = np.array(gt_all)
    gtarr[gtarr == 2] = 1
    gtarr[gtarr == 3] = 0
    return gtarr

def simulate_binary_variant_matrix(
    num_bees: int,
    genome_size: int,
    total_variants: int,
    avg_variants_per_bee: int
) -> np.ndarray:
    """
    Simulate a binary variant presence/absence matrix for a set of individuals.

    Args:
        num_bees: Number of bee samples to simulate.
        genome_size: Size of genome (in base pairs).
        total_variants: Total unique variant positions in the population.
        avg_variants_per_bee: Average number of variants observed per bee.

    Returns:
        A (total_variants x num_bees) binary numpy array:
        - Rows = variants
        - Columns = samples
        - Values = 1 if variant present in individual, 0 otherwise
    """
    # Create the shared variant pool (indexed from 0 to total_variants-1)
    variant_indices = np.arange(total_variants)

    # Initialize the binary matrix
    matrix = np.zeros((total_variants, num_bees), dtype=np.uint8)

    # Fill in the matrix with 1s for sampled variants
    for i in range(num_bees):
        sampled_indices = np.random.choice(variant_indices, avg_variants_per_bee, replace=False)
        matrix[sampled_indices, i] = 1

    return matrix

def compute_discovery_curve_binary_matrix(matrix: np.ndarray, seed=42) -> np.ndarray:
    """
    Compute discovery saturation curve from a binary matrix.
    
    Each row is a variant, each column a sample.
    A value of 1 means the sample carries the variant.
    
    Args:
        matrix (np.ndarray): 2D array of shape (variants, samples), 0/1 values.
        seed (int): Random seed for reproducibility.
        
    Returns:
        np.ndarray: 1D array of cumulative number of variants discovered as more samples are added.
    """
    np.random.seed(seed)
    
    # Random permutation of columns (samples)
    perm = np.random.permutation(matrix.shape[1])
    permuted_matrix = matrix[:, perm]

    # Cumulative max across columns to track if variant has been seen in any sample so far
    cum_discovery = np.maximum.accumulate(permuted_matrix, axis=1)

    # Count how many variants have been discovered (have at least one 1 so far) at each step
    discovered = np.sum(cum_discovery, axis=0)
    return discovered

def run_discovery_iter(matrix: np.ndarray, n_iter: int) -> list():
    """
    Perform repeated variant discovery curve estimation over multiple random orderings.

    This function simulates the accumulation of discovered variants by incrementally
    adding individuals (columns) from a binary variant presence/absence matrix
    in randomized order across multiple iterations.

    Args:
        matrix (np.ndarray): A 2D binary numpy array of shape (variants, samples),
                             where 1 indicates the presence of a variant in a sample,
                             and 0 indicates absence or missing.
        n_iter (int): Number of randomized iterations to perform.

    Returns:
        tuple:
            - meanarr (np.ndarray): A 2D numpy array of shape (samples, 2), where the first
              column is the number of samples added and the second column is the average
              number of discovered variants.
            - disc_arr (np.ndarray): A 2D numpy array of shape (n_iter, samples) containing
              the discovery curve for each iteration.
    """
    disc_iter = []
    for i in range(0,n_iter,1):
        print(f'running iteration{i}...')
        disc = compute_discovery_curve_binary_matrix(matrix, seed=i)
        disc_iter.append(disc)
    disc_arr = np.array(disc_iter)
    meandf = pd.DataFrame(list(range(1,disc_arr.shape[1]+1,1)),np.mean(disc_arr, axis=0)).reset_index()
    meandf = meandf[[0,'index']]
    meanarr = np.array(meandf)
    return meanarr, disc_arr


def fit_saturation_curve(
    data: np.ndarray,
    fit: str,
    bounds: Optional[dict[str, Tuple[Tuple[float, ...], Tuple[float, ...]]]] = None
) -> Tuple[float, ...]:
    """
    Fit a discovery saturation model to cumulative variant data.

    Args:
        data (np.ndarray): 2D array of shape (n_samples, 2), where the first column
            represents the number of samples and the second column the cumulative
            number of unique variants observed.
        fit (str): Type of model to fit:
            - 'mm': Michaelis-Menten
            - 'ne': Negative Exponential
            - 'hs': Hybrid Saturation (exp + linear)
            - 'hp': Hybrid Power (exp + power-law)
        bounds (dict, optional): Dictionary with keys matching the chosen `fit`
            and values as (lower_bounds, upper_bounds).
            If not provided, sensible defaults are used.

    Returns:
        Tuple of fitted parameters.
    """
    x, y = data[:, 0], data[:, 1]

    # ---- default bounds ----
    default_bounds = {
        "mm": ((0, 0), (1e6, 1e4)),
        "ne": ((0, 0), (1e6, 1.0)),
        "hs": ((0, 0, 0), (1e7, 1.0, 1e4)),   # c upper bound increased to 1e4
        "hp": ((0, 0, 0, 0), (1e7, 1.0, 1e4, 1.0)),  # alpha constrained [0,1]
    }

    # pick bounds (user-supplied overrides defaults)
    if bounds and fit in bounds:
        lower, upper = bounds[fit]
    else:
        lower, upper = default_bounds[fit]

    # ---- initial guesses ----
    if fit == 'mm':
        p0 = [max(y), np.median(x)]
        popt, _ = curve_fit(michaelis_menten, x, y, p0=p0, bounds=(lower, upper))
        return tuple(popt)

    elif fit == 'ne':
        p0 = [max(y), 0.05]
        popt, _ = curve_fit(neg_exponential, x, y, p0=p0, bounds=(lower, upper))
        return tuple(popt)

    elif fit == 'hs':
        p0 = [max(y), 0.01, 10.0]
        popt, _ = curve_fit(hybrid_saturation, x, y, p0=p0, bounds=(lower, upper))
        return tuple(popt)

    elif fit == 'hp':
        p0 = [max(y), 0.01, 10.0, 0.8]
        popt, _ = curve_fit(hybrid_power, x, y, p0=p0, bounds=(lower, upper), maxfev=20000)
        return tuple(popt)

    else:
        raise ValueError("fit must be one of: 'mm', 'ne', 'hs', 'hp'")


def evaluate_fit(x: np.ndarray, y: np.ndarray, model_func, params: tuple) -> dict:
    """
    Evaluate a fitted saturation model using AIC and BIC.

    Args:
        x (np.ndarray): Input sample sizes (n individuals).
        y (np.ndarray): Observed cumulative variants.
        model_func (callable): The model function used for fitting.
        params (tuple): Parameters returned by fit_saturation_curve.

    Returns:
        dict: {
            "RSS": Residual sum of squares,
            "AIC": Akaike Information Criterion,
            "BIC": Bayesian Information Criterion,
            "n_params": Number of fitted parameters
        }
    """
    yhat = model_func(x, *params)
    n = len(y)
    k = len(params)
    rss = np.sum((y - yhat) ** 2)

    # log-likelihood under Gaussian errors
    ll = -n/2 * np.log(rss/n) - n/2 * np.log(2*np.pi) - n/2

    aic = 2 * k - 2 * ll
    bic = k * np.log(n) - 2 * ll

    return {"RSS": rss, "AIC": aic, "BIC": bic, "n_params": k}

def compare_models(data: np.ndarray, models: list) -> tuple:
    # model lookup
    model_lookup = {
    "mm": michaelis_menten,
    "ne": neg_exponential,
    "hs": hybrid_saturation,
    "hp": hybrid_power,
}
    fitd = {}
    x, y = data[:,0], data[:,1]
    for model in models:
        fitd.setdefault(model,{'params':fit_saturation_curve(data, fit=model)})
        fitd[model]['eval'] =evaluate_fit(x,y, model_lookup[model], fitd[model]['params'])
    eval_dfs = {}     
    for eval_crit in ['AIC', 'BIC', 'RSS']:
        evaldf = pd.DataFrame(index=models, columns=models)
        for model1 in models:
            for model2 in models:
                evaldf.loc[model1, model2] = fitd[model1]['eval'][eval_crit]-fitd[model2]['eval'][eval_crit]
        eval_dfs[eval_crit] = evaldf    
    return fitd, eval_dfs

def get_predicted_max_var(discarr: np.ndarray, quantile: float, method:Optional[str] = 'hs', bounds: Optional[dict] = None, verbose=True):
    if not method in ['mm', 'ne', 'hs', 'hp']:
        raise IOError('method must be one of "mm", "ne", "hs" or "hp"')
    # model lookup
    model = {
    "mm": michaelis_menten,
    "ne": neg_exponential,
    "hs": hybrid_saturation,
    "hp": hybrid_power,
    }

    fit_val_max = []
    discdf = pd.DataFrame(discarr)
    for i,k in discdf.iterrows():
        data = np.array((np.array(discdf.columns),np.array(k))).T
        params = fit_saturation_curve(data, fit=method)
        if method in ['mm', 'ne']:
            fit_val_max.append(model[method](discarr.shape[1], params[0], params[1]))
        elif method=='hs':
            fit_val_max.append(model[method](discarr.shape[1], params[0], params[1], params[2]))
        else:
            fit_val_max.append(model[method](discarr.shape[1], params[0], params[1], params[2], params[3]))
    quantile_val = np.round(np.quantile(fit_val_max,quantile))
    
    
    if verbose==True:
        print(f' The maximum predicted number of variants is {np.round(np.max(fit_val_max))} \n the {int(quantile*100)}th percentile value is {quantile_val}')
        fig, ax = plt.subplots(ncols=1, nrows=1)
        ax.hist(fit_val_max)
        ax.axvline(quantile_val, color='red')
        plt.show()
    return  quantile_val, fit_val_max, np.max(fit_val_max)