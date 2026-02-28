import os
import re
import pandas as pd
import matplotlib.pyplot as plt

# 你要叠加的三个文件
files = [
    "Data/new/Q200/BCC.Et.00.md3",
    "Data/new/Q1000/dt005/endt10000/BCC.Et.00.md3",
    "Data/new/Q5000/BCC.Et.00.md3",
]

# 可选：如果你不想从文件名解析label，就在这里手动写
labels = ["Q=200", "Q=1000", "Q=5000"]
#labels = None

def read_md3(path: str):
    # 第一行是说明
    with open(path, "r") as f:
        title = f.readline().strip()

    # 读取数据（空格分隔），header 会自动用文件里那行列名
    df = pd.read_csv(path, sep=r"\s+")
    return title, df

def label_from_filename(path: str):
    # 尝试从文件名里抓 Qxxx 作为标签（按你实际命名规则改）
    base = os.path.basename(path)
    m = re.search(r"Q([0-9]+(?:\.[0-9]+)?)", base)  # 支持 Q1 / Q1.5
    if m:
        return f"Q={m.group(1)}"
    return base  # 找不到就用文件名兜底

# 读入所有数据
datasets = []
for i, fp in enumerate(files):
    title, df = read_md3(fp)
    lab = labels[i] if labels is not None else label_from_filename(fp)
    datasets.append((lab, title, df))

# 假设三个文件结构完全一致：列名一致
ref_df = datasets[0][2]
x_col = ref_df.columns[0]
y_cols = ref_df.columns[1:]

# 针对每个 y 列画一张图，但叠加三个文件
for col in y_cols:
    plt.figure(figsize=(10, 4))
    for lab, _, df in datasets:
        plt.plot(df[x_col], df[col], label=lab)

    plt.xlabel(x_col)
    plt.ylabel(col)
    plt.title(col)
    plt.legend()

    plt.ticklabel_format(style="plain", axis="y")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)

    plt.tight_layout()
    plt.savefig(f"Data/BCC.Et.00.{col}.png")
    plt.close()
