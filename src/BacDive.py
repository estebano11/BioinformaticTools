#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BacDive Taxonomy Query Tool
Fetches taxonomic information (phylum and class) for bacterial species from BacDive API
"""

import sys
import requests
from requests.auth import HTTPBasicAuth

class BacDiveTaxonomyClient:
    """Client for querying BacDive taxonomic information"""
    
    def __init__(self, username, password):
        self.headers = {'Accept': 'application/json'}
        self.credentials = HTTPBasicAuth(username, password)
    
    def get_taxon_info(self, search_term):
        """
        Query BacDive API for taxonomic information
        
        Args:
            search_term (str): Genus or species name (min 4 chars)
            
        Returns:
            dict: API response containing taxonomic data or None if request fails
        """
        response = requests.get(
            f'http://bacdive.dsmz.de/api/pnu/taxon/{search_term}/',
            headers=self.headers,
            auth=self.credentials
        )
        
        return response.json() if response.status_code == 200 else None

    def process_species_name(self, raw_name):
        """
        Extract genus name from formatted string
        
        Args:
            raw_name (str): Input in format "gi:123456_Genus_species"
            
        Returns:
            str: Extracted genus name
        """
        parts = raw_name.split('_')
        if len(parts) < 2:
            return None
            
        # Handle special cases
        if parts[1] in ["Candidatus", "MULTISPECIES"] and len(parts) > 2:
            return parts[2]
        return parts[1]

def main():
    if len(sys.argv) != 2:
        print("Usage: python BacDivePNU.py input_file.txt")
        sys.exit(1)
    
    # Load credentials from environment or config in production
    client = BacDiveTaxonomyClient(
        username="your_email@example.com",
        password="your_password"
    )
    
    with open(sys.argv[1], "r") as f:
        for line in f:
            species_name = line.strip()
            genus = client.process_species_name(species_name)
            
            if not genus:
                print(f"{species_name}\tINVALID_FORMAT\tINVALID_FORMAT")
                continue
                
            taxon_data = client.get_taxon_info(genus)
            
            if taxon_data and taxon_data.get('count', 0) > 0:
                result = taxon_data['results'][0]
                print(f"{species_name}\t{result.get('phylum', 'NA')}\t{result.get('classis', 'NA')}")
            else:
                print(f"{species_name}\tNOT_FOUND\tNOT_FOUND")

if __name__ == '__main__':
    main()
