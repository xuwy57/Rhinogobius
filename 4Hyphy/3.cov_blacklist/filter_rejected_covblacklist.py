import os
import pandas as pd

# 输入文件
filtered_output_file = "../exportAln_genespace_rejected1229/filtered_genespace.orthogroups.txt"  # 包含 Group_ID 和基因信息
status_file = "../exportAln_genespace_1229/final_output_table.txt"       # 包含 Group_ID 和需要剔除的物种
rform_dir = "Rform"                          # 包含物种目录和 out_Pvalues.txt 文件
output_dir = "for_rejected1229"              # 输出目录

# 创建输出目录
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# 读取 filtered_output.txt 文件
filtered_df = pd.read_csv(filtered_output_file, sep="\t", dtype=str)

# 提取 Group_ID 和中间数字 ID
filtered_df['Rhinogobius_formosanus_ID'] = filtered_df['Rhinogobius_formosanus'].str.extract(r'\|(\d+)\|')

# 读取 final_output_table.txt 文件
status_df = pd.read_csv(status_file, sep="\t", dtype=str)
pass_ids = status_df[status_df['Status'] == 'PASS']

# 将所有物种的文件一次性加载到内存中
species_data = {}
for species_dir in os.listdir(rform_dir):
    species_path = os.path.join(rform_dir, species_dir)
    if os.path.isdir(species_path):
        input_file = os.path.join(species_path, "out_Pvalues.txt")
        if os.path.exists(input_file):
            species_data[species_dir] = pd.read_csv(input_file, sep="\t", header=None, dtype=str)

# 批量处理所有 PASS 的 OrthoID
for _, row in pass_ids.iterrows():
    group_id = row['ID']
    exclude_species = row['RemovedTips'].split(',') if pd.notna(row['RemovedTips']) else []

    if group_id in filtered_df['OrthoID'].values:
        target_id = filtered_df.loc[filtered_df['OrthoID'] == group_id, 'Rhinogobius_formosanus_ID'].values[0]
        if pd.notna(target_id):  # 确保 target_id 不是 NaN
            target_id = str(target_id).zfill(6)  # 使用 zfill 方法填充前导零

            # 修改所有相关物种的数据
            for species in exclude_species:
                if species in species_data:
                    pvalues_df = species_data[species]
                    target_rows = pvalues_df[pvalues_df[0] == target_id]
                    if not target_rows.empty:
                        pvalues_df.loc[target_rows.index, pvalues_df.columns[-2]] = "0.001"
                        pvalues_df.loc[target_rows.index, pvalues_df.columns[-1]] = "Exclude"

# 保存修改后的文件
for species_dir, pvalues_df in species_data.items():
    output_species_dir = os.path.join(output_dir, species_dir)
    if not os.path.exists(output_species_dir):
        os.makedirs(output_species_dir)

    output_file = os.path.join(output_species_dir, "out_Pvalues.txt")
    pvalues_df.to_csv(output_file, sep="\t", header=False, index=False)

print(f"Processing completed. Modified files are saved in {output_dir}")
