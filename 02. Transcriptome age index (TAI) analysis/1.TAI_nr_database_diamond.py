#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os, sys, glob, time, shutil, datetime, subprocess, collections
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import pandas as pd
from Bio import SeqIO

TAXID_COL = 12
MAX_TARGET_SEQS = 0
CMD_LOG = ".log"
TAB = b'\t'
FAILURES = []


def record_failure(stage, item, err):
    FAILURES.append((stage, str(item), str(err)))
    print(f" {stage} | {item} | {err}", flush=True)


def reset_dir(d, retries=12, wait=0.5):
    for _ in range(retries):
        shutil.rmtree(d, ignore_errors=True)
        try:
            os.makedirs(d, exist_ok=True); return
        except PermissionError:
            time.sleep(wait)
    os.makedirs(d, exist_ok=True)


def log_cmd(cmd):
    with open(CMD_LOG, "a") as f:
        f.write(f"[{datetime.datetime.now():%Y-%m-%d %H:%M:%S}] {' '.join(cmd)}\n")


def read_sp_level(sp_level_file):
    species_levels = {}
    with open(sp_level_file, 'r') as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.rstrip('\r\n').split('\t')
        if len(parts) < 2:
            continue
        species, levels = parts[0].strip(), [lv.strip() for lv in parts[1:] if lv.strip()]
        if not species or not levels:
            continue
        if len(levels) != len(parts) - 1:
            record_failure("read_sp_level", species, f"{len(parts) - 1 - len(levels)} ")
        species_levels[species] = levels
    return species_levels


def run_diamond_blastp(input_pep, output_dir, db, evalue, diamond_threads):
    base_name = os.path.basename(input_pep)
    output_file = os.path.join(output_dir, f"{base_name}.diamond")
    cmd = ["diamond", "blastp", "--db", db, "--query", input_pep, "--out", output_file,
           "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
           "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids",   
           "--evalue", str(evalue), "--sensitive",                                
           "--max-target-seqs", str(MAX_TARGET_SEQS),                            
           "--threads", str(diamond_threads)]                                    
    log_cmd(cmd); print(f"CMD: {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True)                                 
    n, sz = (sum(1 for _ in open(output_file)), os.path.getsize(output_file)) if os.path.exists(output_file) else (0, 0)
    print(f"Success: {base_name}  ({n} , {sz / 1024 ** 3:.2f} GB)", flush=True)


def batch_run_diamond(input_dir, output_dir, db, evalue, diamond_threads, max_workers):
    if not os.path.exists(db):
        raise FileNotFoundError(f"nr no: ./{db}")
    pep_files = sorted(glob.glob(os.path.join(input_dir, '*')))
    print(f"[diamond] 待比对 {len(pep_files)} 个 pep: {[os.path.basename(p) for p in pep_files]}")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:  
        futures = {executor.submit(run_diamond_blastp, p, output_dir, db, evalue, diamond_threads): p for p in pep_files}
        for fu, p in futures.items():
            try:
                fu.result()
            except Exception as e:
                record_failure("diamond", os.path.basename(p), e)


def build_taxid2lineage(lineage_file):
    t2l = {}
    with open(lineage_file) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 2 and p[0].strip() and p[1].strip():
                t2l.setdefault(p[0].strip(), p[1].strip())
    print(f"[{os.path.basename(lineage_file)}] 唯一 taxid = {len(t2l)}")
    return t2l


def resolve_lineage(staxids_str, t2l):
    ids = [t.strip() for t in str(staxids_str or "").split(';') if t.strip()]
    for t in ids:
        if t2l.get(t):
            return t2l[t], t
    return "", (ids[0] if ids else "")


def process_diamond_with_lineage(diamond_file, t2l, output_dir):
    base_name = os.path.basename(diamond_file)
    output_file = os.path.join(output_dir, base_name.replace('.diamond', '_add_lineage.txt'))
    n = nm = 0
    unmatched = collections.Counter()
    with open(diamond_file, newline='') as fin, open(output_file, 'w', newline='') as fout:
        for line in fin:
            raw = line.rstrip('\n').rstrip('\r')
            if not raw:
                continue
            parts = raw.split('\t')
            lineage, used = resolve_lineage(parts[TAXID_COL] if len(parts) > TAXID_COL else "", t2l)
            fout.write(raw + '\t' + lineage + '\n')
            n += 1
            if lineage:
                nm += 1
            else:
                unmatched[used] += 1
    print(f"--- {base_name}\n    {n} | lineage  {nm} ({nm / n * 100:.2f}%) |  {n - nm} ( {len(unmatched)}  taxid)")
    return unmatched


def batch_add_lineage(diamond_dir, lineage_file, output_dir, max_workers):
    t2l = build_taxid2lineage(lineage_file)
    diamond_files = sorted(glob.glob(os.path.join(diamond_dir, '*.diamond')))
    print(f"[lineage] 待处理 {len(diamond_files)} 个 diamond")
    all_unmatched = collections.Counter()
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_diamond_with_lineage, f, t2l, output_dir): f for f in diamond_files}
        for fu, f in futures.items():
            try:
                all_unmatched.update(fu.result())
            except Exception as e:
                record_failure("add_lineage", os.path.basename(f), e)
    if all_unmatched:
        rep = os.path.join(output_dir, "unmatched_taxid_count.tsv")
        with open(rep, "w") as f:
            f.write("taxid\trows\n" + "".join(f"{t}\t{c}\n" for t, c in all_unmatched.most_common()))
        print(f"[lineage] no taxid: {rep} ( {len(all_unmatched)} )")



def compute_per_line_ps(species_levels, input_dir, output_dir):
    for blast_file in sorted(glob.glob(os.path.join(input_dir, "*.txt"))):
        filename = os.path.basename(blast_file)
        sp_name = filename.split('.')[0]
        if sp_name not in species_levels:
            record_failure("compute_ps", filename, f"input_sp_level.txt no {sp_name}, skip")
            continue
        levels = species_levels[sp_name]
        output_file = os.path.join(output_dir, f"{sp_name}_PS.txt")
        n = blank = 0
        try:
            with open(blast_file) as f_in, open(output_file, 'w') as f_out:
                for line in f_in:
                    parts = line.rstrip('\r\n').split('\t')
                    if len(parts) < 2:
                        f_out.write(line); continue
                    n += 1
                    lineage_str = parts[-1].strip()
                    ps_value = len(levels)
                    if lineage_str:
                        marks = [lv.strip() for lv in lineage_str.split(';')]
                        for idx in range(len(levels) - 1, -1, -1):
                            if levels[idx] in marks:
                                ps_value = idx + 1; break
                    else:
                        blank += 1
                    f_out.write(f"{line.strip()}\t{ps_value}\n")
            print(f"Processed: {blast_file} -> {output_file} ({n}, lineage {blank} -> PS{len(levels)})")
        except Exception as e:
            record_failure("compute_ps", filename, e)


def _ps_min_worker(job):
    idx, path, start, end, levels_b, nlev = job
    gene_min, n, blank = {}, 0, 0
    ps_cache = {}
    with open(path, 'rb', buffering=1 << 22) as f:
        if start:
            f.seek(start - 1)
            if f.read(1) != b'\n':
                f.readline()
        first_pos = pos = f.tell()
        while pos < end:
            line = f.readline()
            if not line:
                break
            pos += len(line)
            i = line.find(TAB)
            if i < 0:
                continue
            lin = line[line.rfind(TAB) + 1:].strip()
            if not lin:
                blank += 1
            ps = ps_cache.get(lin)
            if ps is None:
                ps = nlev
                if lin:
                    marks = {m.strip() for m in lin.split(b';')}
                    for k in range(nlev - 1, -1, -1):
                        if levels_b[k] in marks:
                            ps = k + 1; break
                ps_cache[lin] = ps
            g = line[:i]
            if ps < gene_min.get(g, 10 ** 9):
                gene_min[g] = ps
            n += 1
    return idx, gene_min, n, blank, first_pos, pos


def compute_ps_and_min_ps(input_dir, output_dir, species_levels, nprocs, chunk_gb):
    chunk_bytes = int(chunk_gb * 1024 ** 3)
    results = []
    for input_file in sorted(glob.glob(os.path.join(input_dir, "*.txt"))):
        filename = os.path.basename(input_file)
        sp_name = filename.split('.')[0]
        if sp_name not in species_levels:
            record_failure("ps_min", filename, f"input_sp_level.txt"); continue
        levels = species_levels[sp_name]; nlev = len(levels)
        size = os.path.getsize(input_file)
        jobs = [(k, input_file, s, min(s + chunk_bytes, size), [lv.encode() for lv in levels], nlev)
                for k, s in enumerate(range(0, size, chunk_bytes))]          
        gene_min, n_all, blank_all, covered, t0 = {}, 0, 0, 0, time.time()
        print(f"[ps+min] {filename}  {size / 1024 ** 3:.2f} GB ->  {len(jobs)} , {nprocs}", flush=True)
        try:
            with mp.Pool(processes=nprocs) as pool:
                for done, (idx, part, n, blank, p0, p1) in enumerate(pool.imap(_ps_min_worker, jobs), 1):
                    for g, ps in part.items():
                        if ps < gene_min.get(g, 10 ** 9):
                            gene_min[g] = ps
                    n_all += n; blank_all += blank; covered += p1 - p0
                    print(f"    [{done}/{len(jobs)}] #{idx}  {n} |  {n_all} , gene {len(gene_min)}, "
                          f"time {(time.time() - t0) / 60:.1f} min", flush=True)
        except Exception as e:
            record_failure("ps_min", filename, f": {e}"); continue
        bytes_ok = (covered == size)
        if not bytes_ok:
            record_failure("ps_min", filename, f" {covered} !=  {size},")
        output_file = os.path.join(output_dir, f"{sp_name}_PS_no_PS{nlev}.txt")   
        with open(output_file, 'w') as f_out:
            f_out.write("".join(f"{g.decode()}\t{ps}\n" for g, ps in gene_min.items()))
        print(f"[ps+min] {filename} -> {output_file}  (coloum {n_all}, no lineage {blank_all} -> PS{nlev}, "
              f"gene {len(gene_min)},  {'OK' if bytes_ok else 'f'},  {(time.time() - t0) / 60:.1f} min)", flush=True)
        results.append({"file": filename, "lines": n_all, "blank": blank_all,
                        "genes": len(gene_min), "output": output_file, "bytes_ok": bytes_ok})
    return results


def extract_min_phylostratum(input_dir, output_dir, species_levels):
    for input_file in sorted(glob.glob(os.path.join(input_dir, "*"))):
        filename = os.path.basename(input_file)
        sp_name = filename.split('_')[0]
        if sp_name not in species_levels:
            record_failure("min_PS", filename, f"no {sp_name} "); continue
        output_file = os.path.join(output_dir, f"{filename.split('.')[0]}_no_PS{len(species_levels[sp_name])}.txt")  
        gene_min_ps = {}
        try:
            with open(input_file, 'r') as f:
                for line in f:
                    parts = line.rstrip('\r\n').split('\t')
                    if len(parts) >= 2:
                        ps = int(parts[-1])
                        if ps < gene_min_ps.get(parts[0], 10 ** 9):
                            gene_min_ps[parts[0]] = ps
        except Exception as e:
            record_failure("min_PS", filename, e); continue
        with open(output_file, 'w') as f_out:
            f_out.write("".join(f"{g}\t{p}\n" for g, p in gene_min_ps.items()))
        print(f"Processed: {input_file} -> {output_file} ({len(gene_min_ps)} genes)")


def update_gene_ages(input_ps_dir, input_pep_dir, output_dir, species_levels):
    for ps_file in sorted(glob.glob(os.path.join(input_ps_dir, "*"))):
        filename = os.path.basename(ps_file)
        sp_name = filename.split("_")[0]
        if sp_name not in species_levels:
            record_failure("fill_missing", filename, f"no {sp_name}"); continue
        last_ps_value = str(len(species_levels[sp_name]))
        pep_file = os.path.join(input_pep_dir, f"{sp_name}.pep")
        output_file = os.path.join(output_dir, f"{sp_name}_geneid_PS.txt")
        if not os.path.exists(pep_file):
            record_failure("fill_missing", filename, f"pep no: {pep_file}"); continue
        try:
            ps_dict, bad = {}, 0
            with open(ps_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line:
                        continue
                    p = line.split('\t')
                    if len(p) < 2:
                        bad += 1; continue
                    ps_dict[p[0]] = p[1].strip()
            if bad:
                record_failure("fill_missing", filename, f" {bad} ")
            all_genes = {r.id for r in SeqIO.parse(pep_file, "fasta")}
            n_fill = 0
            for gene in all_genes - set(ps_dict):
                ps_dict[gene] = last_ps_value; n_fill += 1
            with open(output_file, 'w') as f_out:
                f_out.write("GeneID\tAge\n" + "".join(f"{g}\t{ps_dict[g]}\n" for g in sorted(ps_dict)))
            print(f"Processed: {ps_file} -> {output_file}  (pep gene {len(all_genes)},  PS {len(all_genes) - n_fill}, fill {last_ps_value}  {n_fill})")
        except Exception as e:
            record_failure("fill_missing", filename, e)


def count_phylostratum(input_dir, output_dir):
    for file_path in sorted(glob.glob(os.path.join(input_dir, "*.txt"))):
        filename = os.path.basename(file_path)
        sp_name = filename.split("_")[0]
        try:
            df = pd.read_csv(file_path, sep='\t')
            if 'Age' not in df.columns:
                raise KeyError(f"no Age: {list(df.columns)}")
            counts = df['Age'].value_counts().sort_index()
            output_file = os.path.join(output_dir, f"{sp_name}_count_PS_number.txt")
            df6 = pd.DataFrame({'phylostratum': counts.index, 'Count': counts.values})
            df6['Percent'] = (df6['Count'] / df6['Count'].sum() * 100).round(4)
            df6.to_csv(output_file, sep='\t', index=False)
            print(f"Completed: {filename} -> {output_file}")
        except Exception as e:
            record_failure("count_ps", filename, e)


if __name__ == "__main__":
    input_pep_files = "input_pep"        
    nr_database_path = "./nr_db/.dmnd"   
    lineage_file = "lineage_v2.txt"            
    sp_level_file = "input_sp_level.txt"      
    
    evalue = 1e-5                              
    diamond_threads = 45                       
    diamond_max_workers = 1                   
    max_workers_add_lineage = 6                
    ps_nprocs = 24                             
    ps_chunk_gb = 10                           
    
    output_1_diamond_result = "output_1_diamond_result"              
    output_2_add_lineage = "output_2_add_lineage"                    
    output_3_add_lineage_num = "output_3_add_lineage_num"            
    output_4_geneid_phylostranum = "output_4_geneid_phylostranum"    
    output_5_add_ps_last = "output_5_add_ps_last"                    
    output_6_PS_count_PS_num = "output_6_PS_count_PS_num"            

    species_levels = read_sp_level(sp_level_file)   
    print(f"Loaded species levels for {len(species_levels)} species:")
    for sp, lv in species_levels.items():
        print(f"  {sp}: {len(lv)} levels -> {lv}")

    for d in (output_4_geneid_phylostranum, output_5_add_ps_last, output_6_PS_count_PS_num):
        reset_dir(d)

    
    batch_run_diamond(input_pep_files, output_1_diamond_result, nr_database_path, evalue, diamond_threads, diamond_max_workers)
    batch_add_lineage(output_1_diamond_result, lineage_file, output_2_add_lineage, max_workers_add_lineage)
    compute_per_line_ps(species_levels, output_2_add_lineage, output_3_add_lineage_num)
    compute_ps_and_min_ps(output_2_add_lineage, output_4_geneid_phylostranum, species_levels, ps_nprocs, ps_chunk_gb)   # ③+④ 合并
    update_gene_ages(output_4_geneid_phylostranum, input_pep_files, output_5_add_ps_last, species_levels)
    count_phylostratum(output_5_add_ps_last, output_6_PS_count_PS_num)

    print("\n" + "=" * 70)
    if FAILURES:                               
        print(f" {len(FAILURES)} :")
        for stage, item, err in FAILURES:
            print(f"   [{stage}] {item}: {err}")
        sys.exit(1)
    print("success")
