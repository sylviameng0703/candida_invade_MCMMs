import warnings
from os import listdir
from os.path import join
from time import time
from sys import argv
import numpy as np
import pandas as pd
import micom
from micom import load_pickle
from micom.workflows import workflow
from pandas.errors import SettingWithCopyWarning
from tqdm import tqdm
from pathlib import Path

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

start = time()
logger = micom.logger

def calc_rates(args):
    f, frac, folder_path = args
    print(args)
    try:
        smp = f.split('.pickle')[0].split('_')[0]
        cond = f.split('.pickle')[0].split('_')[1]
        filepath = Path(folder_path) / f
        logger.info(f"Loading community {smp}-{cond} with tradeoff fraction {frac}...")
        com = load_pickle(filepath)
        abund = com.abundances.copy()
        com.taxa = list(abund.index.values)
        sol = com.cooperative_tradeoff(fraction=frac, fluxes=True, pfba=True, rtol=1e-4, atol=1e-5)

        rates = sol.members[["growth_rate",'abundance']].copy()
        rates = rates.drop('medium', errors='ignore')
        rates.loc["community", :] = [sol.growth_rate, rates['abundance'].sum()]
        rates['sample'] = smp
        rates['condition'] = cond
        rates['tradeoff'] = frac

        fluxes = sol.fluxes.copy()
        fluxes['sample'] = smp
        fluxes['condition'] = cond
        fluxes['tradeoff'] = frac
        return {"rates": rates, "fluxes": fluxes}
    except Exception as e:
        logger.error(f"Issue with sample {f}: {e}")
        return None


def main():
    try:
        max_procs = snakemake.threads
    except NameError:
        max_procs = 16

    # 输入模型 pickle 存放文件夹
    folder_path = "/mnt/nfs/qnapNFSxin/mengya/ca_invade_new/models_invade0.1_species/"
    # folder_path = "/mnt/nfs/qnapNFSxin/mengya/ca_invade/models_nocal_species"
    # folder_path = "/mnt/nfs/qnapNFSxin/mengya/ca_invade/models_novirus_species"

    # 结果输出文件夹
    results_dir = "/mnt/nfs/qnapNFSxin/mengya/ca_invade_new/results/"

    from os import makedirs
    from os.path import exists

    if not exists(results_dir):
        try:
            makedirs(results_dir, exist_ok=True)
            logger.info(f"Created results directory: {results_dir}")
        except Exception as e:
            logger.error(f"Could not create results directory {results_dir}: {e}")
            raise

    files = [f for f in listdir(folder_path) if f.endswith(".pickle")]
    # tradeoffs = [round(x, 2) for x in np.arange(0.1, 1.0, 0.1)]
    tradeoffs = [0.8]

    args = [(f, t, folder_path) for f in files for t in tradeoffs]
    logger.info(f"Max processes: {max_procs}, Total tasks: {len(args)}")

    res = workflow(calc_rates, args, max_procs)

    all_rates = []
    all_fluxes = []

    for i, result in enumerate(tqdm(res, total=len(args), desc="Processing samples")):
        if result is not None:
            all_rates.append(result['rates'])
            all_fluxes.append(result['fluxes'])
        else:
            logger.warning(f"Issue with sample {args[i][0]}")

    today_str = f"{pd.Timestamp('today'):%Y%m%d}"

    if all_rates:
        rates_df = pd.concat(all_rates, sort=False)
        rates_path = join(results_dir, f"{today_str}_invade0.1_growth_rates_tradeoff0.8_pfba.pkl")
        rates_df.to_pickle(rates_path)
        logger.info(f"Saved growth rates to {rates_path}")
    else:
        logger.warning("No rates data collected.")

    if all_fluxes:
        fluxes_df = pd.concat(all_fluxes, sort=False)
        fluxes_path = join(results_dir, f"{today_str}_invade0.1_fluxes_tradeoff0.8_pfba.pkl")
        fluxes_df.to_pickle(fluxes_path)
        logger.info(f"Saved fluxes to {fluxes_path}")
    else:
        logger.warning("No fluxes data collected.")

    elapsed = time() - start
    logger.info(f"Elapsed time: {elapsed:.2f} seconds.")


if __name__ == "__main__":
    main()
