import pandas as pd

# 输入文件
status_file = "../exportAln_genespace_1229/final_output_table.txt"  # 包含 OrthoID、status 和需要排除的物种
orthogroup_file = "../exportAln_genespace_1229/genespace.orthogroups.txt"  # 包含 OrthoID 和物种详细信息
output_file = "filtered_genespace.orthogroups.txt"  # 输出文件

# 读取 status 文件
status_df = pd.read_csv(status_file, sep="\t")

# 筛选 status 为 PASS 的行
pass_ids = status_df[status_df['Status'] == 'PASS'][['ID', 'RemovedTips']]

# 读取 orthogroup 文件
orthogroup_df = pd.read_csv(orthogroup_file, sep="\t")

# **修改点：将 "ID" 替换为实际的列名，如 "OrthoID"**
pass_ids = pass_ids.rename(columns={"ID": "OrthoID"})
if "OrthoID" not in orthogroup_df.columns:
        raise KeyError("Column 'OrthoID' not found in orthogroup_df. Please check the file format.")

# 创建一个结果列表来存储处理后的数据
results = []

# 遍历 PASS 的 OrthoID
for _, row in pass_ids.iterrows():
        ortho_id = row['OrthoID']
        exclude_species = row['RemovedTips'].split(',') if pd.notna(row['RemovedTips']) else []

        # 找到对应的行
        if ortho_id in orthogroup_df['OrthoID'].values:
            ortho_row = orthogroup_df[orthogroup_df['OrthoID'] == ortho_id].iloc[0]
            processed_row = ortho_row.copy()
            # 对需要排除的物种列设置为 NA
            for species in exclude_species:
                for col in orthogroup_df.columns[3:]:  # 假设物种列从第 4 列开始
                    if "Rhinogobius_formosanus" in col:  # 如果是 Rhinogobius_formosanus 列，跳过
                        continue
                    cell_value = processed_row[col]
                    # 确保 cell_value 是字符串或 NaN，避免 float 报错
                    if pd.isna(cell_value):
                        processed_row[col] = "NA"
                    elif isinstance(cell_value, str) and species in cell_value:
                        processed_row[col] = "NA"
            # 添加到结果列表
            results.append(processed_row)

# 转换结果为 DataFrame
filtered_df = pd.DataFrame(results)
# 保存到输出文件
filtered_df.to_csv(output_file, sep="\t", index=False)
print(f"Processing completed. Filtered data saved to {output_file}")
