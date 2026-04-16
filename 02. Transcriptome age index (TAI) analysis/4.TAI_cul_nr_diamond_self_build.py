#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# @Author        : yuzijian
# @Email         : yuzijian1010@163.com
# @FileName      : 4.TAI_cul_nr_diamond_self_build.py
# @Time          : 2025-07-15 11:05:58
# @description   : 计算TAI值并绘图
"""
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import matplotlib.pyplot as plt


def calculate_tai(df, ancestral_level=None):
    phylostratum_col = df.columns[0]
    df[phylostratum_col] = pd.to_numeric(df[phylostratum_col], errors='coerce').fillna(11).astype(int)
    if ancestral_level is not None:
        df = df[df[phylostratum_col] <= ancestral_level].copy()
        print(f"计算祖先状态TAI (PS≤{ancestral_level})，保留{len(df)}个基因")
    phylostratum = df[phylostratum_col].values
    expression_cols = df.columns[2:]
    for col in expression_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    
    expression_data = df[expression_cols].values
    tai_values = {}
    for i, col_name in enumerate(expression_cols):
        expr_values = expression_data[:, i]
        total_expression = np.sum(expr_values)
        if total_expression > 0:
            # 使用公式计算：TAI = Σ(ps_i * (e_is / E_s))
            weighted_sum = 0
            for j in range(len(expr_values)):
                gene_contribution = expr_values[j] / total_expression
                weighted_sum += phylostratum[j] * gene_contribution
            tai = weighted_sum
        else:
            tai = 0
            print(f"警告: 发育时期 {col_name} 所有基因表达量为0，TAI设为0")
        tai_values[col_name] = tai
    return tai_values, phylostratum, expression_cols


def flat_line_test(df, tai_values, n_permutations=1000):
    stages = list(tai_values.keys())
    observed_tai = np.array([tai_values[stage] for stage in stages])
    observed_variance = np.var(observed_tai)
    phylostratum_col = df.columns[0]  # 进化层级列名
    expression_cols = df.columns[2:]  # 表达量列名
    for col in expression_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
    expression_data = df[expression_cols].values
    df[phylostratum_col] = pd.to_numeric(df[phylostratum_col], errors='coerce').fillna(11).astype(int)
    phylostratum = df[phylostratum_col].values
    permuted_variances = []
    for _ in tqdm(range(n_permutations), desc="Execute FlatLineTest"):
        permuted_data = np.copy(expression_data)
        for i in range(permuted_data.shape[0]):
            np.random.shuffle(permuted_data[i, :])
        permuted_tai = []
        total_expression = np.sum(permuted_data, axis=0)
        for j in range(permuted_data.shape[1]):
            expr_values = permuted_data[:, j]
            if total_expression[j] > 0:
                weighted_sum = 0
                for k in range(len(expr_values)):
                    gene_contribution = expr_values[k] / total_expression[j]
                    weighted_sum += phylostratum[k] * gene_contribution
                permuted_tai.append(weighted_sum)
            else:
                permuted_tai.append(0)
        permuted_variances.append(np.var(permuted_tai))
    count = sum(1 for var in permuted_variances if var >= observed_variance)
    p_value = (count + 1) / (n_permutations + 1)  # 使用保守估计避免p=0
    return p_value, permuted_variances, observed_variance


def plot_tai_curve(tai_values, title="TAI variation curve", save_path=None):
    plt.figure(figsize=(12, 6))
    stages = list(tai_values.keys())
    tai_vals = [tai_values[stage] for stage in stages]
    plt.plot(stages, tai_vals, 'o-', color='dodgerblue', linewidth=2, markersize=8)
    plt.xlabel('Developmental period', fontsize=12)
    plt.ylabel('TAI', fontsize=12)
    plt.title(title, fontsize=14)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"图表已保存至: {save_path}")


def plot_histogram(permuted_variances, observed_variance, p_value, save_path=None):
    plt.figure(figsize=(10, 6))
    sns.histplot(permuted_variances, bins=30, color='skyblue', edgecolor='white', kde=True)
    plt.axvline(observed_variance, color='red', linestyle='--', linewidth=2, 
                label=f'Observation variance = {observed_variance:.4f}')
    plt.title(f'FlatLineTest (p-value = {p_value:.4f})', fontsize=14)
    plt.xlabel('Variance', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"图表已保存至: {save_path}")


def plot_ancestral_tai_curves(df, ancestral_levels, expression_cols, save_path=None):
    plt.figure(figsize=(14, 8))
    for level in ancestral_levels:
        tai_values, _, _ = calculate_tai(df.copy(), ancestral_level=level)
        tai_vals = [tai_values[col] for col in expression_cols]
        plt.plot(expression_cols, tai_vals, 'o-', linewidth=2, 
                 label=f'PS≤{level} (Ancestral state)')
    global_tai, _, _ = calculate_tai(df.copy())
    global_tai_vals = [global_tai[col] for col in expression_cols]
    plt.plot(expression_cols, global_tai_vals, 'o--', color='black', linewidth=2, 
             label='Global TAI')
    plt.xlabel('Time', fontsize=12)
    plt.ylabel('TAI', fontsize=12)
    plt.xticks(rotation=45, ha='right', fontsize=10)
    plt.yticks(fontsize=10)
    plt.legend(fontsize=10, loc='best')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"图表已保存至: {save_path}")

def save_tai_results(tai_values, output_file):
    with open(output_file, 'w') as f:
        f.write("Developmental period\tTAI\n")
        for stage, tai in tai_values.items():
            f.write(f"{stage}\t{tai:.6f}\n")
    print(f"TAI结果已保存至: {output_file}")


if __name__ == "__main__":
    # 设置文件路径
    input_file = "new_Paye_4_PS_add_last.txt"
    global_output = "global_tai_results.txt"
    ancestral_output = "ancestral_tai_ps10_results.txt"
    global_plot = "global_tai_curve.png"
    ancestral_plot = "ancestral_tai_ps10_curve.png"
    flatline_hist = "flatline_test_histogram.png"
    ancestral_curves_plot = "ancestral_tai_curves.png"

    print(f"读取输入文件: {input_file}")
    df = pd.read_csv(input_file, sep='\t')

    print("\n计算全局TAI (包含所有基因):")
    global_tai, _, expression_cols = calculate_tai(df.copy())
    save_tai_results(global_tai, global_output)

    plot_tai_curve(global_tai, title="Global TAI", 
                   save_path=global_plot)

    print("\n执行FlatLineTest显著性分析...")
    p_value, permuted_variances, observed_variance = flat_line_test(df.copy(), global_tai, n_permutations=1000)
    print(f"FlatLineTest结果: p值 = {p_value:.4f}")

    plot_histogram(permuted_variances, observed_variance, p_value, save_path=flatline_hist)

    print("\n计算祖先状态TAI (双壳纲祖先 PS10):")
    ancestral_tai, _, _ = calculate_tai(df.copy(), ancestral_level=10)
    save_tai_results(ancestral_tai, ancestral_output)

    plot_tai_curve(ancestral_tai, title="ancestral  TAI (PS≤10)", 
                   save_path=ancestral_plot)

    print("\n绘制不同进化层级的TAI曲线比较...")
    plot_ancestral_tai_curves(df.copy(), ancestral_levels=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], 
                             expression_cols=expression_cols, 
                             save_path=ancestral_curves_plot)
    
    print("\n分析完成!")


