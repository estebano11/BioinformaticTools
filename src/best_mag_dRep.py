#!/usr/bin/env python3
"""
best_mag_dRep.py - Integrated tool for clustering genomes and selecting best MAGs

Combines cluster generation and MAG selection in one workflow:
1. Generates genome clusters from CBD.csv input
2. Selects best MAG from each cluster based on quality metrics
"""

import pandas as pd
import argparse
from pathlib import Path

def generate_clusters(cbd_csv):
    """
    Generate genome clusters from CBD.csv file
    
    Args:
        cbd_csv (str): Path to CBD.csv file (format: genome,cluster)
        
    Returns:
        dict: {cluster: [genome1, genome2,...]}
    """
    clusters = {}
    with open(cbd_csv, "r") as f:
        # Skip header and process each line
        next(f)  
        for line in f:
            genome, cluster = line.strip().split(",")
            if cluster not in clusters:
                clusters[cluster] = []
            clusters[cluster].append(genome)
    return clusters

def best_mag(df, completeness_thresh=90, contam_thresh=5):
    """
    Select best MAG based on quality metrics with thresholds
    
    Args:
        df (DataFrame): MAG quality data
        completeness_thresh (int): Minimum completeness percentage
        contam_thresh (int): Maximum contamination percentage
        
    Returns:
        Series: Best MAG information
    """
    # Sort primarily by completeness, then contamination
    df = df.sort_values(
        by=["completeness", "contamination"],
        ascending=[False, True]
    )
    
    # Apply quality thresholds
    high_quality = df[
        (df["completeness"] >= completeness_thresh) &
        (df["contamination"] <= contam_thresh)
    ]
    
    return high_quality.iloc[0] if not high_quality.empty else df.iloc[0]

def main():
    parser = argparse.ArgumentParser(
        description="Integrated genome clustering and MAG selection tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input files
    parser.add_argument('cbd_csv', help="Input CBD.csv file (genome,cluster)")
    parser.add_argument('genome_info', help="Genome quality CSV: Genome,Completeness,Contamination,Length,N50")
    parser.add_argument('score_file', help="Genome scores CSV (Sdb.csv)")
    
    # Parameters
    parser.add_argument('--key', required=True, 
                       help="Comma-separated identifiers for your MAGs")
    parser.add_argument('-o', '--output', default="best_mags.tsv",
                      help="Output TSV filename")
    parser.add_argument('--comp_thresh', type=int, default=90,
                      help="Completeness threshold (%)")
    parser.add_argument('--contam_thresh', type=int, default=5,
                      help="Contamination threshold (%)")
    
    args = parser.parse_args()
    
    # Load quality data
    print("Loading genome quality data...")
    quality_df = pd.read_csv(args.genome_info, sep=",", index_col=0)
    score_df = pd.read_csv(args.score_file, sep=",", index_col=0)
    key_identifiers = args.key.split(",")
    
    # Generate clusters
    print("Generating genome clusters...")
    clusters = generate_clusters(args.cbd_csv)
    
    # Prepare output dataframe
    output_columns = [
        "cluster", "genomes", "count",
        "best_mag", "best_comp", "best_contam", "best_length", "best_n50", "best_score",
        "key_mag", "key_comp", "key_contam", "key_length", "key_n50", "key_score"
    ]
    results = []
    
    # Process each cluster
    print("Evaluating clusters...")
    for cluster, genomes in clusters.items():
        cluster_mags = []
        key_mags = []
        
        for mag in genomes:
            try:
                mag_data = {
                    "mag": mag,
                    "completeness": quality_df.loc[mag]["completeness"],
                    "contamination": quality_df.loc[mag]["contamination"],
                    "length": quality_df.loc[mag]["length"],
                    "n50": quality_df.loc[mag]["N50"],
                    "score": float(score_df.loc[mag]["score"])
                }
                cluster_mags.append(mag_data)
                
                if any(key in mag for key in key_identifiers):
                    key_mags.append(mag_data)
            except KeyError:
                continue
        
        if not cluster_mags:
            continue
            
        cluster_df = pd.DataFrame(cluster_mags)
        best_mag_info = best_mag(cluster_df, args.comp_thresh, args.contam_thresh)
        
        # Get best key mag if present
        key_mag_info = (
            best_mag(pd.DataFrame(key_mags), args.comp_thresh, args.contam_thresh)
            if key_mags else None
        )
        
        results.append({
            "cluster": cluster,
            "genomes": ",".join(genomes),
            "count": len(genomes),
            "best_mag": best_mag_info["mag"],
            "best_comp": best_mag_info["completeness"],
            "best_contam": best_mag_info["contamination"],
            "best_length": best_mag_info["length"],
            "best_n50": best_mag_info["n50"],
            "best_score": best_mag_info["score"],
            "key_mag": key_mag_info["mag"] if key_mags else "NA",
            "key_comp": key_mag_info["completeness"] if key_mags else "NA",
            "key_contam": key_mag_info["contamination"] if key_mags else "NA",
            "key_length": key_mag_info["length"] if key_mags else "NA",
            "key_n50": key_mag_info["n50"] if key_mags else "NA",
            "key_score": key_mag_info["score"] if key_mags else "NA"
        })
    
    # Save results
    print(f"Saving results to {args.output}...")
    pd.DataFrame(results).to_csv(args.output, sep="\t", index=False)
    print("Done!")

if __name__ == '__main__':
    main()
