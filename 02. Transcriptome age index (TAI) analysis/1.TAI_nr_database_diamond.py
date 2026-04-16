#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import os
import glob
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor


def run_diamond_blastp(input_pep, output_dir, db, evalue, diamond_threads):
    base_name = os.path.basename(input_pep)
    output_file = os.path.join(output_dir, f"{base_name}.diamond")
    cmd = ["diamond", "blastp", "--db", db, "--query", input_pep, "--out", output_file, "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",  "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", "--evalue", str(evalue), "--sensitive", "--threads", str(diamond_threads)]
    try:
        subprocess.run(cmd, check=True)
        print(f"Success: {base_name}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {base_name}: {e}")


def batch_run_diamond(input_dir, output_dir, db, evalue, diamond_threads, max_workers):
    pep_files = glob.glob(os.path.join(input_dir, '*'))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for pep_file in pep_files:
            futures.append(executor.submit(run_diamond_blastp, pep_file, output_dir, db, evalue, diamond_threads))
        for future in futures:
            future.result()


def process_diamond_with_lineage(diamond_file, lineage_file, output_dir):
    base_name = os.path.basename(diamond_file)
    diamond_df = pd.read_csv(diamond_file, sep='\t', header=None)
    lineage_df = pd.read_csv(lineage_file, sep='\t', header=None)
    if len(diamond_df) != len(lineage_df):
        msg = f"Warning: The number of lines in the file does not match! {base_name}: diamond={len(diamond_df)}行, lineage={len(lineage_df)}行"
        print(msg)
    combined_df = diamond_df.copy()
    if lineage_df.shape[1] >= 2:
        combined_df['lineage'] = lineage_df.iloc[:, 1]
    else:
        combined_df['lineage'] = ""
    output_file = os.path.join(output_dir, f"{base_name.replace('.diamond', '_add_lineage.txt')}")
    combined_df.to_csv(output_file, sep='\t', header=False, index=False)


def batch_add_lineage(diamond_dir, lineage_dir, output_dir, max_workers):
    diamond_files = glob.glob(os.path.join(diamond_dir, '*.diamond'))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for diamond_file in diamond_files:
            future = executor.submit(process_diamond_with_lineage, diamond_file, lineage_dir, output_dir)
            futures.append(future)


def compute_per_line_ps(species_levels, input_dir, output_dir):
    for blast_file in glob.glob(os.path.join(input_dir, "*.txt")):
        filename = os.path.basename(blast_file)
        sp_name = filename.split('_')[0]
        if sp_name in species_levels:
            PHYLO_LEVELS = species_levels[sp_name]
        else:
            print(f"Warning: No PHYLO_LEVELS found for {sp_name}, skipping...")
            continue
        
        output_filename = f"{os.path.basename(blast_file).split('.')[0]}_PS.txt"
        output_file = os.path.join(output_dir, output_filename)
        
        with open(blast_file) as f_in, open(output_file, 'w') as f_out:
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    f_out.write(line)
                    continue
                lineage_str = parts[-1].strip()
                if not lineage_str:
                    f_out.write(f"{line.strip()}\t{len(PHYLO_LEVELS)}\n")
                    continue
                try:
                    ps_value = len(PHYLO_LEVELS)
                    for idx in range(len(PHYLO_LEVELS)-1, -1, -1):
                        if PHYLO_LEVELS[idx] in [lvl.strip() for lvl in lineage_str.split(';')]:
                            ps_value = idx + 1
                            break
                    f_out.write(f"{line.strip()}\t{ps_value}\n")
                except:
                    f_out.write(f"{line.strip()}\t{len(PHYLO_LEVELS)}\n")
        print(f"Processed: {blast_file} -> {output_file}")


def extract_min_phylostratum(input_dir, output_dir):
    for input_file in glob.glob(os.path.join(input_dir, "*")):
        # filename = os.path.basename(input_file)
        filename = f"{os.path.basename(input_file).split('.')[0]}_no_PS11.txt"
        output_file = os.path.join(output_dir, filename)
        gene_min_ps = {}
        with open(input_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                gene_id = parts[0]
                ps_value = int(parts[-1])
                if gene_id not in gene_min_ps or ps_value < gene_min_ps[gene_id]:
                    gene_min_ps[gene_id] = ps_value
        with open(output_file, 'w') as f_out:
            for gene_id, min_ps in gene_min_ps.items():
                f_out.write(f"{gene_id}\t{min_ps}\n")
        print(f"Processed: {input_file} -> {output_file} ({len(gene_min_ps)} genes)")


def update_gene_ages(input_ps_dir, input_pep_dir, output_dir, last_ps_value):
    for ps_file in glob.glob(os.path.join(input_ps_dir, "*")):
        try:
            sp_name = os.path.basename(ps_file).split("_")[0]
            pep_file = os.path.join(input_pep_dir, f"{sp_name}.pep")
            output_file = os.path.join(output_dir, f"{sp_name}_geneid_PS.txt")
            with open(ps_file) as f:
                ps_dict = {line.split('\t')[0]: line.split('\t')[1].strip() for line in f if line.strip()}
            all_genes = {record.id for record in SeqIO.parse(pep_file, "fasta")}
            for gene in all_genes - set(ps_dict):
                ps_dict[gene] = last_ps_value

            with open(output_file, 'w') as f_out:
                f_out.write("GeneID\tAge\n")
                f_out.writelines(f"{gene_id}\t{ps_value}\n" for gene_id, ps_value in ps_dict.items())

        except Exception as e:
            print(f"Error processing {ps_file}: {str(e)}")


def count_phylostratum(input_dir, output_dir):
    for file_path in glob.glob(os.path.join(input_dir, "*.txt")):
        try:
            filename = os.path.basename(file_path)
            sp_name = filename.split("_")[0]
            df = pd.read_csv(file_path, sep='\t')
            counts = df['Age'].value_counts().sort_index()
            count_df = pd.DataFrame({
                'phylostratum': counts.index,
                'Count': counts.values
            })
            output_file = os.path.join(output_dir, f"{sp_name}_count_PS_number.txt")
            count_df.to_csv(output_file, sep='\t', index=False)
            print(f"Completed: {filename} -> {output_file}")
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")


if __name__ == "__main__":
    input_pep_files = "input_pep_files"
    nr_database_path = "/gpfshddpool/home/lianshanshan/yuzijian/tools/nr_db/nr_taxdump.dmnd"
    lineage_file = "lineage.txt"
    sp_level_file = "input_sp_level.txt"
    last_ps_value = "11"
    max_workers_add_lineage = 6

    species_levels = {}
    with open(sp_level_file, 'r') as f:
        lines = f.readlines()

    for line in lines[1:]:
        parts = line.strip().split('\t')
        if not parts:
            continue
        species = parts[0]
        levels = [level.strip() for level in parts[1:] if level.strip()]
        species_levels[species] = levels
    
    print(f"Loaded species levels for {len(species_levels)} species:")
    for species, levels in species_levels.items():
        print(f"  {species}: {len(levels)} levels")

    output_1_diamond_result = "output_1_diamond_result"
    output_2_add_lineage = "output_2_add_lineage"
    output_3_add_lineage_num = "output_3_add_lineage_num"
    output_4_geneid_phylostranum = "output_4_geneid_phylostranum"
    output_5_add_ps_last = "output_5_add_ps_last"
    output_6_PS_count_PS_num = "output_6_PS_count_PS_num"

    shutil.rmtree(output_1_diamond_result, ignore_errors=True)
    os.makedirs(output_1_diamond_result)
    shutil.rmtree(output_2_add_lineage, ignore_errors=True)
    os.makedirs(output_2_add_lineage, exist_ok=True)
    shutil.rmtree(output_3_add_lineage_num, ignore_errors=True)
    os.makedirs(output_3_add_lineage_num)
    shutil.rmtree(output_4_geneid_phylostranum, ignore_errors=True)
    os.makedirs(output_4_geneid_phylostranum)
    shutil.rmtree(output_5_add_ps_last, ignore_errors=True)
    os.makedirs(output_5_add_ps_last)
    shutil.rmtree(output_6_PS_count_PS_num, ignore_errors=True)
    os.makedirs(output_6_PS_count_PS_num)

    batch_run_diamond(
        input_dir=input_pep_files,
        output_dir=output_1_diamond_result,
        db=nr_database_path,
        evalue=1e-3,
        diamond_threads=100,
        max_workers=2
    )

    batch_add_lineage(diamond_dir=output_1_diamond_result, lineage_dir=lineage_file, output_dir=output_2_add_lineage, max_workers=max_workers_add_lineage)

    compute_per_line_ps(species_levels, output_2_add_lineage, output_3_add_lineage_num)

    extract_min_phylostratum(output_3_add_lineage_num, output_4_geneid_phylostranum)

    update_gene_ages(input_ps_dir=output_4_geneid_phylostranum, input_pep_dir=input_pep_files, output_dir=output_5_add_ps_last, last_ps_value=last_ps_value)

    count_phylostratum(output_5_add_ps_last, output_6_PS_count_PS_num)

