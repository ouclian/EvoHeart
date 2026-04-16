#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @FileName      : 2.TAI_self_built_database.py
# @Time          : 2025-06-08 09:29:58
# @description   : 使用自建的层级库，获得PS结果
"""
import os
import glob
import shutil
import subprocess
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import concurrent.futures
from functools import partial


def build_blast_db(input_dir, output_dir):
    merged_file = os.path.join(output_dir, "self_blast_db.fasta")
    with open(merged_file, 'w') as outfile:
        for file in glob.glob(os.path.join(input_dir, '*')):
            outfile.write(open(file).read() + '\n')
    subprocess.run(["diamond", "makedb", "--in", merged_file, "--db", os.path.join(output_dir, "self_blast_db") ], check=True)


def run_single_blast(pep_file, db_path, output_dir, evalue, threads):
    base_name = os.path.splitext(os.path.basename(pep_file))[0]
    output_file = os.path.join(output_dir, f"{base_name}_swissprot.blast")
    cmd = ["diamond", "blastp", "--db", db_path, "--query", pep_file, "--out", output_file, "--outfmt", "6", "--evalue", str(evalue), "--sensitive", "--threads", str(threads)]
    subprocess.run(cmd, check=True)
    print(f"Completed: {pep_file} -> {output_file}")


def run_single_blast_star(args):
    return run_single_blast(*args)
def run_diamond_on_pepfiles(input_dir, db_path, output_dir, evalue, threads, max_workers):
    pep_files = glob.glob(os.path.join(input_dir, "*.pep"))
    tasks = [(pep_file, db_path, output_dir, evalue, threads) for pep_file in pep_files]
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(run_single_blast_star, tasks)
        success_count = sum(1 for result in results if result)
        print(f"完成处理: {success_count}/{len(pep_files)} 文件成功处理")


def build_gene_to_species_map(pep_folder):
    gene_to_species = {}
    extensions = [".pep", ".fa", ".fasta", ".faa", ".fna", ".fas"]
    for ext in extensions:
        for pep_file in glob.glob(os.path.join(pep_folder, f"*{ext}")):
            species = os.path.splitext(os.path.basename(pep_file))[0]
            try:
                for record in SeqIO.parse(pep_file, "fasta"):
                    if record.id not in gene_to_species:
                        gene_to_species[record.id] = species
            except Exception as e:
                print(f"警告：解析文件 {pep_file} 时出错: {e}")
    print(f"成功映射 {len(gene_to_species)} 个基因ID到物种")
    return gene_to_species


def add_species_label_to_file(blast_file, output_file, gene_to_species, sp_label, fill_value):
    blast_df = pd.read_csv(blast_file, sep='\t', header=None)
    labels = []
    for _, row in blast_df.iterrows():
        subject_id = row[1]
        species = gene_to_species.get(subject_id)
        label = sp_label.get(species, "0") if species else fill_value
        labels.append(label)
    blast_df[len(blast_df.columns)] = labels
    blast_df.to_csv(output_file, sep='\t', index=False, header=False)
    labeled_count = sum(1 for label in labels if label != "0")
    total_count = len(labels)
    print(f"{os.path.basename(blast_file)}: 标记记录 {labeled_count}/{total_count} ({labeled_count/total_count*100:.1f}%)")
    return labeled_count


def batch_process_blast_labeling(pep_folder, blast_dir, output_dir, sp_label, fill_value):
    gene_to_species = build_gene_to_species_map(pep_folder)
    blast_files = glob.glob(f"{blast_dir}/*")
    for blast_file in blast_files:
        filename = f"{os.path.basename(blast_file).split('.')[0]}_add_PS.blast"
        output_file = os.path.join(output_dir, filename)
        add_species_label_to_file(blast_file, output_file, gene_to_species, sp_label, fill_value)
    print(f"处理完成: 共处理 {len(blast_files)} 个BLAST文件")


def extract_min_phylostratum(input_dir, output_dir):
    for input_file in glob.glob(os.path.join(input_dir, "*")):
        filename = f"{os.path.basename(input_file).split('.')[0]}_no_PS11_self_build_database.txt"
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
            f_out.write("GeneID\tAge\n")
            for gene_id, min_ps in gene_min_ps.items():
                f_out.write(f"{gene_id}\t{min_ps}\n")
        print(f"Processed: {input_file} -> {output_file} ({len(gene_min_ps)} genes)")


def count_phylostratum(input_dir, output_dir):
    for file_path in glob.glob(os.path.join(input_dir, "*.txt")):
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


if __name__ == "__main__":
    input_pep = "input_pep"
    input_self_built_database = "input_self_built_database"
    output_0_self_blast_db = "output_0_self_blast_db"
    database_path = "output_0_self_blast_db/self_blast_db.dmnd"

    output_1_diamond_result = "output_1_diamond_result"
    output_2_blast_PS = "output_2_blast_PS"
    output_3_geneid_phylostranum = "output_3_geneid_phylostranum"
    output_4_PS_count_PS_num = "output_4_PS_count_PS_num"
    fill_value = "11"

    sp_label = {
    "Esco": "1", "Havo": "1", "Bama": "2", "Sace": "2", 
    "Mobr": "3", "Drme": "4", "Neve": "5", "Cael": "6", "Cate": "7", 
    "Hero": "8", "Ocbi": "9", "Magi": "10"
    }

    shutil.rmtree(output_0_self_blast_db, ignore_errors=True)
    os.makedirs(output_0_self_blast_db)
    shutil.rmtree(output_1_diamond_result, ignore_errors=True)
    os.makedirs(output_1_diamond_result)
    shutil.rmtree(output_2_blast_PS, ignore_errors=True)
    os.makedirs(output_2_blast_PS)
    shutil.rmtree(output_3_geneid_phylostranum, ignore_errors=True)
    os.makedirs(output_3_geneid_phylostranum)
    shutil.rmtree(output_4_PS_count_PS_num, ignore_errors=True)
    os.makedirs(output_4_PS_count_PS_num)

    build_blast_db(input_self_built_database, output_0_self_blast_db)

    run_diamond_on_pepfiles(
        input_dir=input_pep,
        db_path=database_path,
        evalue=1e-3,
        output_dir=output_1_diamond_result,
        threads=40,
        max_workers=1
    )
    batch_process_blast_labeling(input_self_built_database, output_1_diamond_result, output_2_blast_PS, sp_label, fill_value)
    extract_min_phylostratum(output_2_blast_PS, output_3_geneid_phylostranum)
    count_phylostratum(output_3_geneid_phylostranum, output_4_PS_count_PS_num)
