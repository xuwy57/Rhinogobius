import re
import sys
import os

def extract_by_groupid(input_file, output_file, target_group):
    current_entry = []
    found = False

    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                # 检测到新条目时处理已收集的条目
                if found:
                    f_out.write(''.join(current_entry))
                    found = False
                current_entry = [line]

                # 使用正则表达式提取GroupID
                group_match = re.search(r'GroupId:([^\s]+)', line)
                if group_match and group_match.group(1) == target_group:
                    found = True
            else:
                if found and line.strip():  # 仅保留有效序列行
                    current_entry.append(line)

        # 处理文件末尾的最后一个条目
        if found:
            f_out.write(''.join(current_entry))

def batch_extract(fasta_file, groupid_file, output_dir):
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    # 读取所有GroupID
    with open(groupid_file, 'r') as f:
        groupids = [line.strip() for line in f if line.strip()]
    
    # 为每个GroupID提取序列
    for groupid in groupids:
        output_file = os.path.join(output_dir, f"{groupid}.fasta")
        extract_by_groupid(fasta_file, output_file, groupid)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extractAA.py <input_fasta> <groupid_file> <output_dir>")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    groupid_file = sys.argv[2]
    output_dir = sys.argv[3]
    
    batch_extract(input_fasta, groupid_file, output_dir)
