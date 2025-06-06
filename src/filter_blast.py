#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

def filter_blast(input_data, identity_threshold=99.0, coverage_threshold=97.0, length_difference=750.0):
    """
    Filter BLAST results based on quality thresholds and remove reciprocal matches.
    
    Args:
        input_data: BLAST results as a list of strings or file-like object
        identity_threshold: Minimum identity percentage (default: 99.0)
        coverage_threshold: Minimum coverage percentage (default: 97.0)
        length_difference: Maximum allowed length difference (default: 750.0)
    
    Returns:
        Tuple of (filtered results, duplicates)
    """
    # Read input if it's a file-like object
    if hasattr(input_data, 'read'):
        input_data = input_data.readlines()
    
    # Initial filtering
    filtered = []
    for line in input_data:
        parts = line.strip().split('\t')
        if len(parts) < 7:  # Ensure we have enough columns
            continue
            
        # Skip self-matches
        if parts[0] == parts[1]:
            continue
            
        # Apply quality filters
        if (float(parts[2]) >= identity_threshold and 
            float(parts[3]) >= coverage_threshold and 
            abs(float(parts[5]) - float(parts[6])) <= abs(length_difference)):
            filtered.append(line)
    
    # Remove exact duplicates
    unique_results = set(filtered)
    
    # Find reciprocal matches (A->B and B->A)
    pairs = set()
    duplicates = set()
    
    for line in unique_results:
        parts = line.strip().split('\t')
        pair = (parts[0], parts[1])
        reverse_pair = (parts[1], parts[0])
        
        if reverse_pair in pairs:
            duplicates.add(line)
            duplicates.add('\t'.join([reverse_pair[0], reverse_pair[1]] + parts[2:]))
        pairs.add(pair)
    
    # Find all unique nodes
    all_nodes = set()
    duplicate_nodes = set()
    
    for line in unique_results:
        parts = line.strip().split('\t')
        for node in (parts[0], parts[1]):
            if node in all_nodes:
                duplicate_nodes.add(node)
            all_nodes.add(node)
    
    return (unique_results - duplicates, duplicates)

def main():
    # Read from stdin or a file if provided
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            filtered, duplicates = filter_blast(f)
    else:
        filtered, duplicates = filter_blast(sys.stdin)
    
    # Write outputs
    with open('filtered_blast.tab', 'w') as f:
        f.writelines(filtered)
    
    with open('duplicates.tab', 'w') as f:
        f.write('\n'.join(duplicates))
    
    print(f"Filtered results saved to filtered_blast.tab")
    print(f"Duplicate matches saved to duplicates.tab")
    print(f"Total filtered results: {len(filtered)}")
    print(f"Total duplicates found: {len(duplicates)}")

if __name__ == '__main__':
    main()
