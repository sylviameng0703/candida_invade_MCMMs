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
logger = micom.logger
logger.setLevel(logging.ERROR)

try:
    max_procs = snakemake.threads
except NameError:
    max_procs = 16

# -------------------------
# Inputs
# -------------------------
TAXONOMY_FILE = "/home/worker021/GSMM/mt_project/data/taxonomy_invaded0.1_species_5e-5.csv"
MEDIA_FILE = "/home/worker021/GSMM/mt_project/data/SVM_gapfilled_species_new.csv"
OUT_DIR = "/mnt/nfs/qnapNFSxin/mengya/ca_invade_new/models_invade0.1_species_medianew/"

PATH_AGORA = "/mnt/nfs/qnapNFSxin/mengya/database/AGORA201/"
PATH_FUNGI = "/mnt/nfs/qnapNFSxin/mengya/database/carvefungi/mapped_models"

PREFIX = "0.1cal"

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


# 添加模型路径列
taxonomy["file"] = taxonomy.apply(make_file_paths, axis=1)

# 检查模型路径是否存在,删除路径缺失项
missing_paths = taxonomy[taxonomy["file"].isna()]
if not missing_paths.empty:
    print(f"[WARN] {len(missing_paths)} species missing model file paths!")
taxonomy = taxonomy.dropna(subset=["file"])

# -------------------------
# Load media
# -------------------------
media = pd.read_csv(MEDIA_FILE)
media.index = media.reaction
media_flux = media.flux

# -------------------------
# Model build function
# -------------------------
def build_and_save(args):
    try:
        s, tax, prefix, out_dir = args
        fname = join(out_dir, f"{s}_{prefix}-invaded.pickle")
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
        found = media_flux.index.isin(ex_ids).sum()
        logger.info("%d/%d import reactions found in model.", found, len(media))

        com.media = media_flux[media_flux.index.isin(ex_ids)]
        com.to_pickle(fname)

        t1 = datetime.datetime.now()
        dt = (t1 - t0).total_seconds()
        print(f"[OK] {s}: saved → {fname} ({dt:.1f}s)")

    except Exception as e:
        err = traceback.format_exc(limit=3)
        print(f"[ERROR] {s}: {type(e).__name__} - {e}\n{err}")


# -------------------------
# run workflow
# -------------------------
if __name__ == '__main__':

    samples = taxonomy.sample_id.unique()
    args = [(s, taxonomy[taxonomy.sample_id == s], PREFIX, OUT_DIR) for s in samples]

    print("=== MICOM model builder started ===")
    print(f'taxonomy行数: {len(taxonomy)}')
    print(f"Total samples: {len(samples)} | Using {max_procs} threads")

    workflow(build_and_save, args, max_procs)
    print("\n=== All tasks finished. ===")


## -------------------------
## workflow test (8 random samples)
## -------------------------
#if __name__ == "__main__":
#    all_samples = taxonomy.sample_id.unique()
#    np.random.seed(36)  # 确保随机选择可复现
#    test_samples = np.random.choice(all_samples, size=min(10, len(all_samples)), replace=False)
#    print(f" Selected {len(test_samples)} samples for test run:")
#    print(test_samples)
#
#    args = [(s, taxonomy[taxonomy.sample_id == s], PREFIX, OUT_DIR) for s in test_samples]
#
#    print("\n=== MICOM model builder started (TEST MODE, 10 samples) ===")
#    print(f'taxonomy行数: {len(taxonomy)}')
#    print(f"Total samples: {len(test_samples)} | Using {max_procs} threads\n")
#
#    workflow(build_and_save, args, max_procs)
#
#    print("\n=== Test build finished. ===")
