import glob
import io
import re

import matplotlib.pyplot as plt
import pandas as pd


def file_index(path: str) -> int:
    match = re.search(r"\.(\d+)\.md3$", path)
    if not match:
        return -1
    return int(match.group(1))


def load_et_file(path: str) -> pd.DataFrame:
    # Read as bytes first so we can safely strip unexpected NUL bytes.
    with open(path, "rb") as f:
        raw = f.read()
    nul_count = raw.count(b"\x00")
    if nul_count:
        print(f"[warn] {path}: removed {nul_count} NUL bytes")
        raw = raw.replace(b"\x00", b"")

    text = raw.decode("utf-8", errors="ignore")
    return pd.read_csv(io.StringIO(text), sep=r"\s+", engine="python")


files = sorted(glob.glob("Data/BCC.Et.*.md3"), key=file_index)
if not files:
    raise FileNotFoundError("No files matched: Data/BCC.Et.*.md3")

frames = [load_et_file(path) for path in files]
df = pd.concat(frames, ignore_index=True)

# Ensure full trajectory order even if file order or write order had glitches.
x_col = df.columns[0]
df = df.sort_values(by=x_col).reset_index(drop=True)

print(f"Loaded {len(files)} files, total rows = {len(df)}")

# Plot each column against time using the merged data.
for col in df.columns[1:]:
    plt.figure(figsize=(10, 5))
    plt.plot(df[x_col], df[col], linewidth=0.8)
    plt.xlabel(x_col)
    plt.ylabel(col)
    plt.title(f"BCC all files: {col}")
    plt.ticklabel_format(style="plain", axis="y")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.tight_layout()
    plt.savefig(f"Data/BCC.Et.ALL.{col}.png", dpi=160)
    plt.close()

# Force + Momentum in one figure.
force_momentum_cols = ["F_all_x", "F_all_y", "F_all_z", "P_all_x", "P_all_y", "P_all_z"]
missing = [c for c in force_momentum_cols if c not in df.columns]
if missing:
    print(f"[warn] Missing columns for combined plot: {missing}")
else:
    plt.figure(figsize=(12, 6))
    for col in force_momentum_cols:
        plt.plot(df[x_col], df[col], linewidth=0.9, label=col)
    plt.xlabel(x_col)
    plt.ylabel("value")
    plt.title("BCC all files: F_all_* and P_all_*")
    plt.ticklabel_format(style="plain", axis="y")
    plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
    plt.legend()
    plt.tight_layout()
    plt.savefig("Data/BCC.Et.ALL.FP.png", dpi=180)
    plt.close()
