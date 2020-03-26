import pandas as pd
import numpy as np
import limix
import glob
import re
import csv


# get file paths and other info
phen_paths = glob.glob('data/hpc/heritability/*_prewas_output_with_ancestral_reconstruction_given_tree_phenotype.csv')
phen_paths.sort()

shared_ref_multi_paths = glob.glob('data/hpc/heritability/*shared_ref_multi.csv')
shared_ref_multi_paths.sort()
shared_maj_multi_paths = glob.glob('data/hpc/heritability/*shared_maj_multi.csv')
shared_maj_multi_paths.sort()
shared_anc_multi_paths = glob.glob('data/hpc/heritability/*shared_anc_multi.csv')
shared_anc_multi_paths.sort()

shared_ref_nomulti_paths = glob.glob('data/hpc/heritability/*shared_ref_nomulti.csv')
shared_ref_nomulti_paths.sort()
shared_maj_nomulti_paths = glob.glob('data/hpc/heritability/*shared_maj_nomulti.csv')
shared_maj_nomulti_paths.sort()
shared_anc_nomulti_paths = glob.glob('data/hpc/heritability/*shared_anc_nomulti.csv')
shared_anc_nomulti_paths.sort()

datasets = [re.sub('data/hpc/heritability/|_prewas.*csv','',x) for x in phen_paths]
basename = [re.sub('data/hpc/heritability/|_phenotype.*csv','',x) for x in phen_paths]

# initialize lists
basenames = []
phen_type = []
ref_type = []
multi_type = []
her_est = []

# calculate heritability estimates
for i in range(len(phen_paths)):
    phenotypes = pd.read_csv(phen_paths[i])
    
    shared_ref_multi = pd.read_csv(shared_ref_multi_paths[i],index_col=0).to_numpy()
    shared_maj_multi = pd.read_csv(shared_maj_multi_paths[i],index_col=0).to_numpy()
    shared_anc_multi = pd.read_csv(shared_anc_multi_paths[i],index_col=0).to_numpy()
    
    shared_ref_nomulti = pd.read_csv(shared_ref_nomulti_paths[i],index_col=0).to_numpy()
    shared_maj_nomulti = pd.read_csv(shared_maj_nomulti_paths[i],index_col=0).to_numpy()
    shared_anc_nomulti = pd.read_csv(shared_anc_nomulti_paths[i],index_col=0).to_numpy()
    
    bm_phen = phenotypes[['BM_pheno_1']].to_numpy()
    wn_phen = phenotypes[['WN_pheno_1']].to_numpy()
    
    bm_ref_multi = limix.her.estimate(bm_phen,'normal',shared_ref_multi,verbose=False)
    bm_maj_multi = limix.her.estimate(bm_phen,'normal',shared_maj_multi,verbose=False)
    bm_anc_multi = limix.her.estimate(bm_phen,'normal',shared_anc_multi,verbose=False)
    wn_ref_multi = limix.her.estimate(wn_phen,'normal',shared_ref_multi,verbose=False)
    wn_maj_multi = limix.her.estimate(wn_phen,'normal',shared_maj_multi,verbose=False)
    wn_anc_multi = limix.her.estimate(wn_phen,'normal',shared_anc_multi,verbose=False)
    
    bm_ref_nomulti = limix.her.estimate(bm_phen,'normal',shared_ref_nomulti,verbose=False)
    bm_maj_nomulti = limix.her.estimate(bm_phen,'normal',shared_maj_nomulti,verbose=False)
    bm_anc_nomulti = limix.her.estimate(bm_phen,'normal',shared_anc_nomulti,verbose=False)
    wn_ref_nomulti = limix.her.estimate(wn_phen,'normal',shared_ref_nomulti,verbose=False)
    wn_maj_nomulti = limix.her.estimate(wn_phen,'normal',shared_maj_nomulti,verbose=False)
    wn_anc_nomulti = limix.her.estimate(wn_phen,'normal',shared_anc_nomulti,verbose=False)
    
    [basenames.append(x) for x in [basename[i]]*12]
    [phen_type.append(x) for x in ['bm']*6]
    [phen_type.append(x) for x in ['wn']*6]
    [ref_type.append(x) for x in ['ref','ref','maj','maj','anc','anc']*2]
    [multi_type.append(x) for x in ['multi','nomulti']*6]
    [her_est.append(x) for x in [bm_ref_multi,bm_ref_nomulti, 
                                 bm_maj_multi,bm_maj_nomulti,
                                 bm_anc_multi,bm_anc_nomulti,
                                 wn_ref_multi,wn_ref_nomulti,
                                 wn_maj_multi,wn_maj_nomulti,
                                 wn_anc_multi,wn_anc_nomulti]]

rows = zip(basenames,phen_type,ref_type,multi_type,her_est)

# save heritability estimates
with open('data/local/heritability/her_ests_shared.csv', "w") as f:
    writer = csv.writer(f)
    for row in rows:
        writer.writerow(row)
