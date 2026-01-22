import re

# 输入文件
fasta_file = "ret_cleanDNA_F52.fasta"
status_file = "../exportAln_genespace_1229/final_output_table.txt"
output_fasta = "filtered_fasta.fasta"

# 读取 final_output_table.txt 文件
group_removal = {}
with open(status_file, "r") as status:
    next(status)  # 跳过标题行
    for line in status:
        columns = line.strip().split("\t")
        group_id = columns[0]
        removed_tips = columns[1].split(",") if columns[1] else []
        group_removal[group_id] = set(removed_tips)

# 处理 Fasta 文件
with open(fasta_file, "r") as fasta, open(output_fasta, "w") as output:
    write_flag = False  # 控制是否将当前序列写入输出
    current_group = None
    current_species = None

    for line in fasta:
        if line.startswith(">"):  # 这是一个序列名
            # 从序列名中提取 GroupID 和物种信息
            current_group = line.strip().split(';')[-1].split(':')[1]
            current_species = line.strip().split(';')[1].split(':')[1]

            # 检查是否需要剔除该物种
            if current_group in group_removal and current_species in group_removal[current_group]:
                print(current_species)
                write_flag = False  # 跳过该物种的序列
            else:
                write_flag = True  # 保留该物种的序列
                output.write(line)  # 写入序列头
        else:
            if write_flag:
                output.write(line)  # 写入序列内容

print(f"Processing completed. Filtered Fasta file saved to {output_fasta}")

