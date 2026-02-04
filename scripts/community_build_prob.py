import micom
from micom import Community
from micom.workflows import workflow
from micom import load_pickle
from os import makedirs
from os.path import join, isfile
import pandas as pd
import numpy as np
import traceback
import datetime
import sys
import logging
import time

# -------------------------
# Global configuration
# -------------------------
try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 16  # 手动运行时的默认并行线程数

logger = micom.logger
logger.setLevel(logging.ERROR)

# -------------------------
# Inputs
# -------------------------
TAXONOMY_FILE = "/home/worker021/GSMM/mt_project/data/taxonomy_species_probiotic_treat-0.3_1.csv"
MEDIA_FILE = "/home/worker021/GSMM/mt_project/data/SVM_gapfilled_species_new.csv"
OUT_DIR = "/mnt/nfs/qnapNFSxin/mengya/ca_invade_new/models_0.1cal-0.3prob-treat1_+preb1"

PATH_AGORA = "/mnt/nfs/qnapNFSxin/mengya/database/AGORA201/"
PATH_FUNGI = "/mnt/nfs/qnapNFSxin/mengya/database/carvefungi/mapped_models/"

PREFIX1 = "0.1cal"
PREFIX2 = "0.3prob"

# -------------------------
# Load data
# -------------------------
makedirs(OUT_DIR, exist_ok=True)
taxonomy = pd.read_csv(TAXONOMY_FILE).query("relative > 5e-5")

def make_file_paths(row):
    """根据 Kingdom 选择相应模型库"""
    try:
        if row.Kingdom == "Bacteria":
            return [join(PATH_AGORA, i) for i in row.file.split("|")]
        elif row.Kingdom == "Fungi":
            return [join(PATH_FUNGI, j) for j in row.file.split("|")]
    except Exception:
        return np.nan

taxonomy["file"] = taxonomy.apply(make_file_paths, axis=1)

# 检查模型路径是否存在
missing_paths = taxonomy[taxonomy["file"].isna()]
if not missing_paths.empty:
    print(f"[WARN] {len(missing_paths)} species missing model file paths!")
taxonomy = taxonomy.dropna(subset=["file"])

# -------------------------
# Load media
# -------------------------
media = pd.read_csv(MEDIA_FILE)
media.index = media.reaction
media = media.flux

# -------------------------
# Model build function
# -------------------------
def build_and_save(args):
    try:
        s, tax, prefix1, prefix2, out_dir = args
        fname = join(out_dir, f"{s}_{prefix1}-{prefix2}.pickle")
        # 如果结果文件已存在 → 直接跳过
        if isfile(fname):
            print(f"[SKIP] {s}: file already exists.")
            return
        # 模型构建
        t0 = datetime.datetime.now()
        print(f" Building sample: {s}")

        com = Community(tax, id=s, progress=False)

        # 培养基，仅保留在模型反应列表中的部分
        ex_ids = [r.id for r in com.exchanges]
        found = media.index.isin(ex_ids).sum()
        logger.info("%d/%d import reactions found in model.", found, len(media))

        com.media = media[media.index.isin(ex_ids)]
        com.to_pickle(fname)

        t1 = datetime.datetime.now()
        dt = (t1 - t0).total_seconds()
        print(f"[OK] {s}: saved → {fname} ({dt:.1f}s)")

    except Exception as e:
        err = traceback.format_exc(limit=3)
        print(f"[ERROR] {s}: {type(e).__name__} - {e}\n{err}")

# -------------------------
# Run workflow
# -------------------------
if __name__ == "__main__":
   start = time.time()
   samples = taxonomy.sample_id.unique()
   args = [(s, taxonomy[taxonomy.sample_id == s], PREFIX1, PREFIX2, OUT_DIR) for s in samples]

   print(f"=== MICOM model builder started ===")
   print(f"Total samples: {len(samples)} | Using {max_procs} threads")

   workflow(build_and_save, args, max_procs)

   elapsed = time.time() - start
   print("\n=== All tasks finished. ===")
   print(f" Completed all tasks in {elapsed:.2f} seconds.")

# # -------------------------
# # workflow test (10 random samples)
# # -------------------------
# if __name__ == "__main__":
#     all_samples = taxonomy.sample_id.unique()
#     np.random.seed(42)  # 确保随机选择可复现
#     test_samples = np.random.choice(all_samples, size=min(10, len(all_samples)), replace=False)
#     print(f" Selected {len(test_samples)} samples for test run:")
#     print(test_samples)
#
#     args = [(s, taxonomy[taxonomy.sample_id == s], PREFIX1, PREFIX2, OUT_DIR) for s in test_samples]
#
#     print("\n=== MICOM model builder started (TEST MODE, 10 samples) ===")
#     print(f"Total samples: {len(test_samples)} | Using {max_procs} threads\n")
#
#     workflow(build_and_save, args, max_procs)
#
#     print("\n=== Test build finished. ===")
