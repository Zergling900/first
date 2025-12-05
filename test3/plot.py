import pandas as pd
import matplotlib.pyplot as plt

filename = "Data/BCC.Et.00.md3"
#filename = "Data/Diamond.Et.00.md3"
#filename = "Data/FCC.Et.00.md3"

# 第一行是列名说明，不作为数据
with open(filename, "r") as f:
    dataset_title = f.readline().strip()

# 读取实际数据（空格分隔）
df = pd.read_csv(filename, sep=r"\s+")

# 第一列做横轴
x_col = df.columns[0]

# plot each column except the first
for col in df.columns[1:]:
    plt.figure()
    plt.plot(df[x_col], df[col])

    plt.xlabel(x_col)   # x 轴使用第一列的列名
    plt.ylabel(col)     # y 轴使用该列列名
    plt.title(col)      # 图标题使用该列列名

    plt.ticklabel_format(style='plain', axis='y')   # 纵轴用普通数字
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    
    plt.tight_layout()
    plt.savefig(f"Data/BCC.Et.00.{col}.png")
    #plt.savefig(f"Data/Diamond.Et.00.{col}.png")
    #plt.savefig(f"Data/FCC.Et.00.{col}.png")
    plt.close()
