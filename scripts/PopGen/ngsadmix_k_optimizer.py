#!/usr/bin/env python3
"""
Script to determine the optimal K value from NGSadmix outputs (K=1 to K=6)
This script analyzes log-likelihood values and calculates delta K to find the best K.
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import glob

def detect_k_values(base_path=".", file_pattern="run_K*_rep*.log"):
    """
    Automatically detect K values from NGSadmix log files.
    
    Args:
        base_path (str): Base directory containing NGSadmix outputs
        file_pattern (str): Pattern to match log files
        
    Returns:
        tuple: (min_k, max_k) range of K values found
    """
    base_path = Path(base_path)
    k_values = set()
    
    # Find all files matching the pattern
    pattern_path = base_path / file_pattern
    log_files = list(base_path.glob(file_pattern))
    
    if not log_files:
        print(f"Warning: No files found matching pattern '{file_pattern}' in '{base_path}'")
        return None, None
    
    print(f"Found {len(log_files)} log files matching pattern")
    
    # Extract K values from filenames
    k_pattern = re.compile(r'[kK](\d+)', re.IGNORECASE)
    
    for log_file in log_files:
        filename = log_file.name
        match = k_pattern.search(filename)
        if match:
            k_val = int(match.group(1))
            k_values.add(k_val)
            print(f"  Found K={k_val} in file: {filename}")
    
    if not k_values:
        print("Error: No K values could be extracted from filenames")
        return None, None
    
    min_k = min(k_values)
    max_k = max(k_values)
    
    print(f"\nDetected K range: {min_k} to {max_k}")
    print(f"K values found: {sorted(k_values)}")
    
    return min_k, max_k

def parse_ngsadmix_log(log_file):
    """
    Parse NGSadmix log file to extract log-likelihood value.
    
    Args:
        log_file (str): Path to NGSadmix log file
        
    Returns:
        float: Log-likelihood value, or None if not found
    """
    try:
        with open(log_file, 'r') as f:
            content = f.read()
            # Look for log-likelihood pattern in NGSadmix output
            # Pattern may vary, adjust regex as needed for your output format
            patterns = [
                r'best like=([+-]?\d+\.?\d*)',
                r'loglike:\s*([+-]?\d+\.?\d*)',
                r'Log likelihood:\s*([+-]?\d+\.?\d*)',
                r'like=([+-]?\d+\.?\d*)'
            ]
            
            for pattern in patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    return float(match.group(1))
                    
        print(f"Warning: Could not find log-likelihood in {log_file}")
        return None
        
    except FileNotFoundError:
        print(f"Error: Log file {log_file} not found")
        return None
    except Exception as e:
        print(f"Error parsing {log_file}: {e}")
        return None

def calculate_delta_k(log_likelihoods):
    """
    Calculate delta K for Structure-like analysis.
    Delta K = mean(|L''(K)|) / sd(L(K))
    
    Args:
        log_likelihoods (dict): Dictionary with K values as keys and log-likelihood lists as values
        
    Returns:
        dict: Delta K values for each K
    """
    k_values = sorted(log_likelihoods.keys())
    delta_k = {}
    
    # Calculate mean log-likelihood for each K
    mean_lnpk = {}
    for k in k_values:
        if log_likelihoods[k]:
            mean_lnpk[k] = np.mean(log_likelihoods[k])
        else:
            mean_lnpk[k] = None
    
    # Calculate delta K
    for i, k in enumerate(k_values):
        if k == k_values[0] or k == k_values[-1]:
            delta_k[k] = None  # Cannot calculate for first and last K
            continue
            
        if (mean_lnpk[k-1] is not None and 
            mean_lnpk[k] is not None and 
            mean_lnpk[k+1] is not None):
            
            # Second derivative approximation
            l_double_prime = abs(mean_lnpk[k+1] - 2*mean_lnpk[k] + mean_lnpk[k-1])
            
            # Standard deviation of L(K)
            if len(log_likelihoods[k]) > 1:
                sd_lnpk = np.std(log_likelihoods[k])
                if sd_lnpk > 0:
                    delta_k[k] = l_double_prime / sd_lnpk
                else:
                    delta_k[k] = None
            else:
                delta_k[k] = None
        else:
            delta_k[k] = None
    
    return delta_k

def analyze_ngsadmix_results(base_path=".", file_pattern="run_K*_rep*.log", k_range=None):
    """
    Analyze NGSadmix results for different K values.
    
    Args:
        base_path (str): Base directory containing NGSadmix outputs
        file_pattern (str): Pattern to match log files
        k_range (tuple): Range of K values to analyze (min, max). If None, auto-detect from files
        
    Returns:
        dict: Analysis results
    """
    results = {
        'k_values': [],
        'log_likelihoods': {},
        'mean_log_likelihood': {},
        'best_log_likelihood': {},
        'delta_k': {},
        'best_k_likelihood': None,
        'best_k_delta': None
    }
    
    base_path = Path(base_path)
    
    # Auto-detect K range if not provided
    if k_range is None:
        print("Auto-detecting K values from filenames...")
        min_k, max_k = detect_k_values(base_path, file_pattern)
        if min_k is None or max_k is None:
            print("Error: Could not detect K values from files")
            return results
        k_range = (min_k, max_k)
    else:
        print(f"Using provided K range: {k_range[0]} to {k_range[1]}")
    
    # Collect log-likelihood values for each K
    for k in range(k_range[0], k_range[1] + 1):
        results['k_values'].append(k)
        results['log_likelihoods'][k] = []
        
        # Look for log files with the pattern run_K{k}_rep{rep}.log
        possible_files = []
        
        # Check for multiple replicates (common practice)
        for rep in range(1, 21):  # Check up to 20 replicates
            possible_files.extend([
                base_path / f"run_K{k}_rep{rep}.log",
                base_path / f"run_k{k}_rep{rep}.log",
            ])
        
        # Also use glob to find files matching the general pattern for this K
        k_specific_patterns = [
            f"run_K{k}_rep*.log",
            f"run_k{k}_rep*.log",
            f"*K{k}_*.log",
            f"*k{k}_*.log"
        ]
        
        for pattern in k_specific_patterns:
            possible_files.extend(list(base_path.glob(pattern)))
        
        # Remove duplicates while preserving order
        seen = set()
        unique_files = []
        for f in possible_files:
            if f not in seen:
                seen.add(f)
                unique_files.append(f)
        
        log_files_found = [f for f in unique_files if f.exists()]
        
        if not log_files_found:
            print(f"Warning: No log files found for K={k}")
            continue
        
        # Parse all found log files for this K
        for log_file in log_files_found:
            ll = parse_ngsadmix_log(log_file)
            if ll is not None:
                results['log_likelihoods'][k].append(ll)
        
        if results['log_likelihoods'][k]:
            results['mean_log_likelihood'][k] = np.mean(results['log_likelihoods'][k])
            results['best_log_likelihood'][k] = max(results['log_likelihoods'][k])
            print(f"K={k}: Found {len(results['log_likelihoods'][k])} replicates, "
                  f"Best LL = {results['best_log_likelihood'][k]:.2f}, "
                  f"Mean LL = {results['mean_log_likelihood'][k]:.2f}")
        else:
            print(f"K={k}: No valid log-likelihood values found")
    
    # Calculate delta K
    results['delta_k'] = calculate_delta_k(results['log_likelihoods'])
    
    # Determine best K based on highest log-likelihood
    valid_k_ll = {k: ll for k, ll in results['best_log_likelihood'].items() if ll is not None}
    if valid_k_ll:
        results['best_k_likelihood'] = max(valid_k_ll, key=valid_k_ll.get)
    
    # Determine best K based on delta K (highest delta K)
    valid_delta_k = {k: dk for k, dk in results['delta_k'].items() if dk is not None}
    if valid_delta_k:
        results['best_k_delta'] = max(valid_delta_k, key=valid_delta_k.get)
    
    return results

def plot_results(results, output_file="ngsadmix_k_analysis.png"):
    """
    Plot log-likelihood and delta K results.
    
    Args:
        results (dict): Results from analyze_ngsadmix_results()
        output_file (str): Output file for the plot
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # Plot 1: Log-likelihood vs K
    k_vals = []
    mean_ll = []
    best_ll = []
    
    for k in results['k_values']:
        if k in results['mean_log_likelihood'] and results['mean_log_likelihood'][k] is not None:
            k_vals.append(k)
            mean_ll.append(results['mean_log_likelihood'][k])
            best_ll.append(results['best_log_likelihood'][k])
    
    if k_vals:
        ax1.plot(k_vals, mean_ll, 'b-o', label='Mean Log-Likelihood', markersize=6)
        ax1.plot(k_vals, best_ll, 'r-s', label='Best Log-Likelihood', markersize=6)
        ax1.set_xlabel('K')
        ax1.set_ylabel('Log-Likelihood')
        ax1.set_title('Log-Likelihood vs K')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_xticks(k_vals)
    
    # Plot 2: Delta K vs K
    delta_k_vals = []
    delta_k_values = []
    
    for k in results['k_values']:
        if k in results['delta_k'] and results['delta_k'][k] is not None:
            delta_k_vals.append(k)
            delta_k_values.append(results['delta_k'][k])
    
    if delta_k_vals:
        ax2.plot(delta_k_vals, delta_k_values, 'g-^', label='Delta K', markersize=6)
        ax2.set_xlabel('K')
        ax2.set_ylabel('Delta K')
        ax2.set_title('Delta K vs K (Evanno Method)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.set_xticks(delta_k_vals)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved as {output_file}")
    plt.show()

def print_summary(results):
    """
    Print summary of the analysis.
    
    Args:
        results (dict): Results from analyze_ngsadmix_results()
    """
    print("\n" + "="*60)
    print("NGSadmix K Value Analysis Summary")
    print("="*60)
    
    print("\nLog-Likelihood Results:")
    print("-" * 40)
    for k in results['k_values']:
        if k in results['best_log_likelihood'] and results['best_log_likelihood'][k] is not None:
            n_reps = len(results['log_likelihoods'][k])
            best_ll = results['best_log_likelihood'][k]
            mean_ll = results['mean_log_likelihood'][k]
            print(f"K={k}: {n_reps} replicates, Best LL = {best_ll:.2f}, Mean LL = {mean_ll:.2f}")
    
    print("\nDelta K Results:")
    print("-" * 40)
    for k in results['k_values']:
        if k in results['delta_k'] and results['delta_k'][k] is not None:
            print(f"K={k}: Delta K = {results['delta_k'][k]:.4f}")
        else:
            print(f"K={k}: Delta K = N/A")
    
    print("\nRecommendations:")
    print("-" * 40)
    if results['best_k_likelihood']:
        print(f"Best K based on highest log-likelihood: K = {results['best_k_likelihood']}")
    if results['best_k_delta']:
        print(f"Best K based on highest Delta K: K = {results['best_k_delta']}")
    
    if results['best_k_likelihood'] and results['best_k_delta']:
        if results['best_k_likelihood'] == results['best_k_delta']:
            print(f"\nBoth methods suggest K = {results['best_k_likelihood']}")
        else:
            print(f"\nMethods disagree: Consider both K = {results['best_k_likelihood']} "
                  f"and K = {results['best_k_delta']}")
    print("="*60)

def main():
    """
    Main function to run the analysis.
    Modify the parameters below to match your file structure.
    """
    # Configuration - modify these parameters as needed
    BASE_PATH = "."  # Directory containing NGSadmix outputs
    FILE_PATTERN = "run_K*_rep*.log"  # Pattern to match log files (uses glob pattern)
    K_RANGE = None  # Set to None for auto-detection, or specify tuple like (1, 6)
    
    print("Starting NGSadmix K value analysis...")
    print(f"Base path: {BASE_PATH}")
    print(f"File pattern: {FILE_PATTERN}")
    
    if K_RANGE is None:
        print("K range: Auto-detect from filenames")
    else:
        print(f"K range: {K_RANGE[0]} to {K_RANGE[1]} (manually specified)")
    
    # Run analysis
    results = analyze_ngsadmix_results(
        base_path=BASE_PATH,
        file_pattern=FILE_PATTERN,
        k_range=K_RANGE
    )
    
    # Check if any results were found
    if not any(results['log_likelihoods'].values()):
        print("\nError: No valid log-likelihood values found for any K value.")
        print("Please check:")
        print("1. File paths and naming pattern")
        print("2. Log file format and content")
        print("3. Base directory location")
        return
    
    # Print summary
    print_summary(results)
    
    # Create plots
    plot_results(results)
    
    # Save results to CSV
    df_data = []
    for k in results['k_values']:
        row = {'K': k}
        if k in results['best_log_likelihood']:
            row['Best_LogLikelihood'] = results['best_log_likelihood'].get(k)
            row['Mean_LogLikelihood'] = results['mean_log_likelihood'].get(k)
            row['N_Replicates'] = len(results['log_likelihoods'][k]) if results['log_likelihoods'][k] else 0
        row['Delta_K'] = results['delta_k'].get(k)
        df_data.append(row)
    
    df = pd.DataFrame(df_data)
    df.to_csv('ngsadmix_k_analysis.csv', index=False)
    print("\nResults saved to 'ngsadmix_k_analysis.csv'")

if __name__ == "__main__":
    main()