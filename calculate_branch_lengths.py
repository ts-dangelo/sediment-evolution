from ete3 import Tree
import os
import pandas as pd

def total_branch_length(treefile):
    """ traverse tree and add up branch distances """
    tree = Tree(treefile, format=1)
    return sum(node.dist for node in tree.traverse())

def process_gene_clusters(tree_dir):
  
    results = []
    
    tree_files = {f.split(".")[0]: os.path.join(tree_dir, f) 
                  for f in os.listdir(tree_dir) if f.endswith(".contree")}
    
    for gene_cluster, treefile in tree_files.items():
        total_length = total_branch_length(treefile)
        results.append({
            "gene_cluster": gene_cluster,
            "total_branch_length": total_length
        })
    
    return pd.DataFrame(results)

tree_dir = "../prodigal_mmseqs//aero_misag_gene_aln_sc_iqtree/" # directory with iqtree output files for gene clusters
proc_gen_clust = process_gene_clusters(tree_dir)
length_annots_df = pd.merge(proc_gen_clust, annotation_df, on='gene_cluster', how='inner') #add tree lengths to annotation dataframe 
