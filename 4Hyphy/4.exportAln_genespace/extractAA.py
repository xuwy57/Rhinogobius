import re

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

# 使用示例
extract_by_groupid('ret_AA_F52.fasta', 'Group_19484_1.fasta', 'Group_19484_1')
