#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import glob
import shutil
import pandas as pd
from Bio import SeqIO


def merge_and_minimize_ps(df_self, df_swiss):
    df_self.columns = ['GeneID', 'Age']
    df_swiss.columns = ['GeneID', 'Age']
    combined_df = pd.concat([df_self, df_swiss], ignore_index=True)
    combined_df['Age'] = pd.to_numeric(combined_df['Age'], errors='coerce')
    valid_df = combined_df.dropna(subset=['Age'])
    result_df = valid_df.sort_values('Age').drop_duplicates('GeneID', keep='first')
    result_df = result_df.sort_values('GeneID')
    return result_df


def process_species_files(input_self_build_database, input_swiss_dir, output_dir, species_abbrev):
    self_file_path = os.path.join(input_self_build_database, f"{species_abbrev}_swissprot_200_add_PS_no_PS11_self_build_database.txt")
    try:
        df_self = pd.read_csv(self_file_path, sep='\t', header=0)
    except:
        df_self = pd.read_csv(self_file_path, sep='\t', header=None)

    swiss_pattern = os.path.join(input_swiss_dir, f"{species_abbrev}*")
    print(swiss_pattern)
    swiss_files = glob.glob(swiss_pattern)
    
    if not swiss_files:
        print(f"警告: 未找到 {species_abbrev} 的SwissProt文件")
        return None

    swiss_file_path = swiss_files[0]
    try:
        df_swiss = pd.read_csv(swiss_file_path, sep='\t', header=0)
    except:
        df_swiss = pd.read_csv(swiss_file_path, sep='\t', header=None)

    result_df = merge_and_minimize_ps(df_self, df_swiss)
    output_path = os.path.join(output_dir, f"{species_abbrev}_merged_PS.txt")
    result_df.to_csv(output_path, sep='\t', index=False)
    return output_path


def update_gene_ages(input_ps_dir, input_pep_dir, output_dir, fill_value):
    for ps_file in glob.glob(os.path.join(input_ps_dir, "*")):
        sp_name = os.path.basename(ps_file).split("_")[0]
        pep_file = os.path.join(input_pep_dir, f"{sp_name}.pep")
        output_file = os.path.join(output_dir, f"{sp_name}_geneid_PS.txt")
        with open(ps_file) as f:
            ps_dict = {line.split('\t')[0]: line.split('\t')[1].strip() for line in f if line.strip()}
        all_genes = {record.id for record in SeqIO.parse(pep_file, "fasta")}
        for gene in all_genes - set(ps_dict):
            ps_dict[gene] = fill_value
        with open(output_file, 'w') as f_out:
            f_out.writelines(f"{gene_id}\t{ps_value}\n" for gene_id, ps_value in ps_dict.items())


def count_phylostratum(input_dir, output_dir):
    for file_path in glob.glob(os.path.join(input_dir, "*.txt")):
        filename = os.path.basename(file_path)
        sp_name = filename.split("_")[0]
        df = pd.read_csv(file_path, sep='\t')
        counts = df['Age'].value_counts().sort_index()
        count_df = pd.DataFrame({'phylostratum': counts.index, 'Count': counts.values})
        output_file = os.path.join(output_dir, f"{sp_name}_count_PS_number.txt")
        count_df.to_csv(output_file, sep='\t', index=False)
        print(f"Completed: {filename} -> {output_file}")


if __name__ == "__main__":
    input_pep = "input_pep"
    fill_value = "11"
    input_self_build_database = "input_self_build_database"
    input_swiss_dir = "input_swissport"
    output_1_cat_self_swissport = "output_1_cat_self_swissport2"
    output_2_add_ps_last = "output_2_add_ps_last"
    output_3_PS_count_PS_num = "output_3_PS_count_PS_num"

    shutil.rmtree(output_1_cat_self_swissport, ignore_errors=True)
    os.makedirs(output_1_cat_self_swissport)
    shutil.rmtree(output_2_add_ps_last, ignore_errors=True)
    os.makedirs(output_2_add_ps_last)
    shutil.rmtree(output_3_PS_count_PS_num, ignore_errors=True)
    os.makedirs(output_3_PS_count_PS_num)

    self_files = glob.glob(os.path.join(input_self_build_database, "*.txt"))
    for self_file in self_files:
        filename = os.path.basename(self_file)
        species_abbrev = filename.split("_")[0]
        print(f"\n处理物种: {species_abbrev}")
        output_path = process_species_files(input_self_build_database, input_swiss_dir, output_1_cat_self_swissport, species_abbrev)
    print("\n处理完成!")
    update_gene_ages(input_ps_dir=output_1_cat_self_swissport, input_pep_dir=input_pep, output_dir=output_2_add_ps_last, fill_value=fill_value)
    count_phylostratum(output_2_add_ps_last, output_3_PS_count_PS_num)
