#!/usr/bin/env python3
"""
search_pairs_of_otus.py - Find pairs of OTUs that exist in both network edges and OTU list

Input:
1. Network edge file (TSV format: OTU1\tOTU2)
2. OTU list file (one OTU per line)

Output:
All edges where both OTUs exist in the OTU list (TSV format: OTU1\tOTU2)
"""

import sys

def load_otus(otu_file):
    """Load OTU identifiers from file into a set for fast lookup"""
    with open(otu_file, "r") as f:
        return {line.strip() for line in f if line.strip()}

def load_edges(edge_file):
    """Load network edges from TSV file"""
    edges = []
    with open(edge_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                edges.append((parts[0], parts[1]))
    return edges

def find_matching_pairs(edges, otus):
    """Filter edges to only those where both nodes are in the OTU set"""
    return [(src, dst) for src, dst in edges if src in otus and dst in otus]

def main():
    if len(sys.argv) != 3:
        print("Usage: python search_pairs_of_otus.py <network_edges.tsv> <otu_list.txt>")
        sys.exit(1)

    edge_file = sys.argv[1]
    otu_file = sys.argv[2]

    try:
        # Load data
        otus = load_otus(otu_file)
        edges = load_edges(edge_file)
        
        # Find and output matching pairs
        matching_pairs = find_matching_pairs(edges, otus)
        for src, dst in matching_pairs:
            print(f"{src}\t{dst}")
            
        # Print summary statistics
        print(f"\nFound {len(matching_pairs)}/{len(edges)} edges with both OTUs in the list", 
              file=sys.stderr)
    
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
