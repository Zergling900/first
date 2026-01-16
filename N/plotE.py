import pandas as pd
import matplotlib.pyplot as plt

file1 = "Data/BCC.Et.04.md3"
file2 = "Data/BCC.Et.05.md3"

# 第一行是列名，不需要 skiprows
df1 = pd.read_csv(file1, sep=r"\s+")
df2 = pd.read_csv(file2, sep=r"\s+")

# 第一列是横轴
time1 = df1.iloc[:, 0]
time2 = df2.iloc[:, 0]

# 第三列（索引 2）是能量 E
E1 = df1.iloc[:, 2]
E2 = df2.iloc[:, 2]

# 绘图
plt.figure()

plt.plot(time1, E1, label="File 1 E")
plt.plot(time2, E2, label="File 2 E")

plt.xlabel(df1.columns[0])  # 自动使用 "time"
plt.ylabel(df1.columns[2])  # 自动使用 "E"
plt.title("Energy")
plt.legend()
plt.tight_layout()

plt.savefig("E1.png")
