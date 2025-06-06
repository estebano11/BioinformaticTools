#!/usr/bin/env python3
"""
cazy_analyzer.py - Comprehensive CAZyme analysis pipeline

Performs:
1. Visualization of extracellular CAZyme distribution across bins
2. Clustering analysis of CAZyme abundance profiles
3. Substrate-product relationship mapping via KEGG
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from Bio.KEGG import REST, Enzyme
import re
from pathlib import Path
import argparse

# Constants
COLORS = ((1.0, 1.0, 1.0, 1.0), (0.0, 0.8, 0.0, 1.0))
SAMPLES = ["R1-008","R1-008b","R1-020","R1-030","R1-034","R1-050","R1-069",
           "R1-083","R1-114","R1-125","R1-128","R1-135","R1-140","R1-161",
           "R2-020","R2-030","R2-034","R2-050","R2-069","R2-083","R2-114",
           "R2-125","R2-128","R2-135","R2-140","R2-161","R3-020","R3-030",
           "R3-034","R3-050","R3-069","R3-083","R3-114","R3-125","R3-128",
           "R3-135","R3-140","R3-161","R4-020","R4-030","R4-034","R4-050",
           "R4-083","R4-114","R4-125","R4-128","R4-135","R4-140","R4-161"]

def load_data(base_path):
    """Load CAZyme annotation data"""
    gram_files = {
        'gram_neg': "extracellular_gram_negativas.caz.bin.tab",
        'gram_pos': "extracellular_gram_positivas.caz.bin.tab"
    }
    
    data = {}
    for key, file in gram_files.items():
        try:
            data[key] = pd.read_csv(
                Path(base_path) / file, 
                sep="\t", 
                header=None, 
                names=["CAZ", "BIN"]
            )
        except FileNotFoundError:
            print(f"Warning: {file} not found")
            data[key] = pd.DataFrame(columns=["CAZ", "BIN"])
    
    return pd.concat([data['gram_neg'], data['gram_pos']], ignore_index=True)

def create_heatmap(data, output_file):
    """Generate presence/absence heatmap of CAZymes across bins"""
    data_dummies = pd.get_dummies(data['CAZ'])
    df = pd.concat([data["BIN"], data_dummies], axis=1)
    result = df.groupby('BIN').sum()
    
    cmap = LinearSegmentedColormap.from_list('Custom', COLORS, len(COLORS))
    plt.figure(figsize=(16,5))
    ax = sns.heatmap(
        result, 
        cmap=cmap, 
        linewidths=.5, 
        linecolor='lightgray',
        xticklabels=True, 
        yticklabels=True
    )
    plt.yticks(rotation=0)
    
    # Customize colorbar
    colorbar = ax.collections[0].colorbar
    colorbar.set_ticks([0.25, 0.75])
    colorbar.set_ticklabels(['NO', 'YES'])
    
    plt.savefig(output_file, dpi=1200)
    plt.close()

def analyze_abundance(base_path, samples, output_file):
    """Analyze CAZyme abundance across samples"""
    caz_abundance = Path(base_path) / "uniq.extracellular.caz.abundance.tab"
    caz = pd.read_csv(caz_abundance, sep="\t", header=None, names=samples)
    
    heatmap = sns.clustermap(
        caz,
        metric="euclidean",
        standard_scale=0,
        method="ward",
        cmap="mako",
        col_cluster=False
    )
    heatmap.savefig(output_file)
    plt.close()

class KEGGAnalyzer:
    """Handles substrate-product relationship mapping using KEGG"""
    
    def __init__(self):
        self.compounds = self._load_kegg_compounds()
        self.trash = self._create_exclusion_list()
    
    def _load_kegg_compounds(self):
        """Load KEGG compound and glycan databases"""
        c = REST.kegg_list("compound")
        cpds = [row.replace("cpd:","").replace("\t", ";").replace("\n","").replace("  "," ") 
                for row in list(c)]
        cp_df = pd.DataFrame(sub.split(";") for sub in cpds)
        
        glycans = REST.kegg_list("glycan")
        gly = [row.replace("gl:","").replace("\t", ";").replace("\n","").replace("  "," ") 
               for row in list(glycans)]
        gly_df = pd.DataFrame(sub.split(";") for sub in gly)
        
        return pd.concat([cp_df, gly_df])
    
    def _create_exclusion_list(self):
        """Create list of compounds to exclude"""
        exclude_indices = list(range(10)) + [11,12,14,16,19,22,26,27,29]
        return self.compounds.loc[exclude_indices]
    
    def _clean_compound_name(self, name):
        """Clean compound names from KEGG"""
        name = re.sub(r'.*C', 'C', name)
        return name.replace("(n)","").replace("(n+m)","").replace("(m)","")
    
    def get_compound_names(self, compound_list):
        """Get proper names for compounds from KEGG"""
        valid_compounds = []
        names = []
        
        for cpd in compound_list:
            clean_cpd = self._clean_compound_name(cpd)
            
            # Skip excluded compounds
            if len(self.trash[self.trash[0] == clean_cpd]) > 0:
                continue
                
            # Look up in compound database
            matches = self.compounds[self.compounds[0] == clean_cpd]
            if len(matches) > 0:
                index = int(str(matches[1]).split(" ")[0])
                names.append(matches[1][index])
                valid_compounds.append(clean_cpd)
        
        return valid_compounds, names
    
    def map_reactions(self, ec_file, output_file):
        """Map EC numbers to reactions and substrates/products"""
        # Load CAZyme-EC mappings
        caz_ec = pd.read_csv(ec_file, sep="\t", header=None)
        
        # Expand multiple EC numbers per GH
        ec_data = []
        for _, row in caz_ec.iterrows():
            ecs = row[1].split(" ") if " " in row[1] else [row[1]]
            for ec in ecs:
                ec_data.append({"GH": row[0], "EC": ec})
        
        ec_df = pd.DataFrame(ec_data)
        unique_ecs = ec_df["EC"].unique()
        
        # Process each EC number
        results = []
        for ec in unique_ecs:
            print(f"Processing {ec}")
            ec_data = self._process_ec_number(ec)
            results.extend(ec_data)
        
        # Save results
        result_df = pd.DataFrame(results, columns=[
            "EC", "Reaction_ID", "Reaction", 
            "Substrate_IDs", "Product_IDs",
            "Substrate_Names", "Product_Names"
        ])
        result_df.to_csv(output_file, sep='\t', index=False)
    
    def _process_ec_number(self, ec):
        """Process a single EC number to get reaction data"""
        ec_results = []
        kegg_ec = f"ec:{ec}"
        
        try:
            # Get linked reactions from KEGG
            reactions = list(REST.kegg_link("reaction", kegg_ec))
            
            if not reactions:
                # If no linked reactions, try direct EC query
                ec_info = self._get_ec_info(ec)
                if ec_info:
                    ec_results.append(self._format_ec_result(ec, ec_info))
            else:
                # Process each reaction
                for reaction in reactions:
                    if "ec" not in reaction[:2]:  # Skip malformed entries
                        rx_data = self._process_reaction(ec, reaction)
                        if rx_data:
                            ec_results.append(rx_data)
        
        except Exception as e:
            print(f"Error processing {ec}: {e}")
        
        return ec_results
    
    def _get_ec_info(self, ec):
        """Get reaction info directly from EC number"""
        try:
            request = REST.kegg_get(ec).read()
            return request.replace(";\n", "")
        except:
            return None
    
    def _format_ec_result(self, ec, ec_info):
        """Format results from direct EC query"""
        lines = ec_info.split("\n")
        data = {"REACTION": "", "SUBSTRATE": "", "PRODUCT": ""}
        
        for line in lines:
            parts = line.split(maxsplit=1)
            if len(parts) > 1 and parts[0] in data:
                data[parts[0]] = parts[1]
        
        return {
            "EC": ec,
            "Reaction_ID": "",
            "Reaction": data["REACTION"],
            "Substrate_IDs": "",
            "Product_IDs": "",
            "Substrate_Names": data["SUBSTRATE"],
            "Product_Names": data["PRODUCT"]
        }
    
    def _process_reaction(self, ec, reaction):
        """Process a single reaction entry"""
        reaction = reaction.strip().replace("\n","").replace("rn:","").split("\t")[-1]
        try:
            reaction_data = REST.kegg_get(reaction).read()
            return self._parse_reaction_data(ec, reaction, reaction_data)
        except Exception as e:
            print(f"Error processing reaction {reaction}: {e}")
            return None
    
    def _parse_reaction_data(self, ec, rx_id, rx_data):
        """Parse reaction data from KEGG response"""
        rx_dict = {}
        current_key = None
        
        for line in rx_data.split("\n"):
            if line.startswith(" "):
                if current_key:
                    rx_dict[current_key] += " " + line.strip()
            else:
                parts = line.split(maxsplit=1)
                if parts:
                    current_key = parts[0]
                    rx_dict[current_key] = parts[1] if len(parts) > 1 else ""
        
        if "EQUATION" not in rx_dict:
            return None
            
        equation = rx_dict["EQUATION"].split("<=>")
        substrates, sub_names = self.get_compound_names(
            equation[0].replace(" + ", ",").replace(" ", "").split(",")
        )
        products, prod_names = self.get_compound_names(
            equation[1].replace(" + ", ",").replace(" ", "").split(",")
        )
        
        return {
            "EC": ec,
            "Reaction_ID": rx_id,
            "Reaction": rx_dict.get("NAME", ""),
            "Substrate_IDs": ",".join(substrates),
            "Product_IDs": ",".join(products),
            "Substrate_Names": ",".join(sub_names),
            "Product_Names": ",".join(prod_names)
        }

def main():
    parser = argparse.ArgumentParser(
        description="CAZyme analysis pipeline for biogas research"
    )
    parser.add_argument('base_path', help="Path to directory containing input files")
    parser.add_argument('--ec_file', default="extracelulares-ec.tab",
                      help="File mapping GH families to EC numbers")
    parser.add_argument('--output_dir', default="results",
                      help="Directory to save output files")
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output_dir).mkdir(exist_ok=True)
    
    print("Loading and processing CAZyme data...")
    caz_data = load_data(args.base_path)
    
    print("Generating heatmap...")
    create_heatmap(
        caz_data,
        Path(args.output_dir) / "cazyme_distribution_heatmap.svg"
    )
    
    print("Analyzing abundance profiles...")
    analyze_abundance(
        args.base_path,
        SAMPLES,
        Path(args.output_dir) / "cazyme_abundance_clustermap.svg"
    )
    
    print("Mapping substrate-product relationships...")
    kegg = KEGGAnalyzer()
    kegg.map_reactions(
        Path(args.base_path) / args.ec_file,
        Path(args.output_dir) / "substrate_product_mappings.tsv"
    )
    
    print("Analysis complete!")

if __name__ == "__main__":
    main()
