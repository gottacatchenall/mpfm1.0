#! /usr/bin/env python

import csv, argparse, os, pandas


# =============================================
#  Globals, etc.
# =============================================

loci_path = 'loci.csv'

local_ld_path = 'linkage.csv'
global_ld_path = 'global_linkage.csv'
fst_path = 'f_st.csv'
allele_freq_path = 'allele_freq.csv'

att_mig_path = 'attempted_migration.csv'
eff_mig_path = 'eff_migration.csv'
patches_path = 'patches.csv'

params_path = 'params.ini'

dynamics_path = 'dynamics.csv'
demography_path = 'demography.csv'
metadata_path = 'metadata.csv'

demography_csv_header = 'run_id, generation, fst_mean, fst_sigma, global_ld_mean, global_ld_sigma, local_ld_mean, local_ld_sigma, n_loci_fixed, neutral, fitness\n'
dynamics_csv_header = 'run_id, generation, eff_mig_mean, eff_mig_sd, att_mig_mean, att_mig_sd, prop_patches_extinct, prop_of_k_mean, prop_of_k_sd\n'

# =============================================
#  Source
# =============================================

def get_folder():
    parser = argparse.ArgumentParser(description='pass the folder which contains the runs you want to aggregate')
    parser.add_argument('folder', nargs=1)
    args = parser.parse_args()
    folder = args.folder[0]
    if not os.path.exists(folder):
        print('invalid path!')
        exit(-1)
    return folder

def init_files(folder):
    dem = folder + '/' + demography_path
    dyn = folder + '/' + dynamics_path
    meta = folder + '/' + metadata_path

    with open(dem, 'a') as file:
        file.write(demography_csv_header)
    with open(dyn, 'a') as file:
        file.write(demography_csv_header)
    with open(meta, 'a') as file:
        file.write('run_id,')

        p = os.path.abspath(folder) + '/' + os.listdir(folder)[0] + '/params.ini'

        print 'okay'
        df = pandas.read_csv(p, sep=' ', error_bad_lines=False, header=None)
        for row in df.itertuples():
            file.write(str(row[1]) + ',')
        file.write('\n')

def write_metadata(run_id, source_dir_path, target_dir_path):
    param_file_path = source_dir_path + '/' + params_path
    target_file_path = target_dir_path + '/' + metadata_path

    df = pandas.read_csv(param_file_path, error_bad_lines=False, sep=' ', header=None)
    with open(target_file_path, 'a') as file:
        file.write(str(run_id) + ',')

        for row in df.itertuples():
            file.write(str(row[2]) + ',')

        file.write('\n')

def get_fst(data):
    return data.F_st.mean(), data.F_st.std()

def get_ld(data):
    return data.D.mean(), data.D.std()

def get_att_mig(data):
    return data.prop_of_new_patch_migrants.mean(), data.prop_of_new_patch_migrants.std()

def get_eff_mig(data):
    return data.eff_migration.mean(), data.eff_migration.std()

def get_prop_extinct(data):
    nrow = data.shape[0]
    return float(data.query('n_indiv == 0').shape[0])/float(nrow)


def get_prop_k(data):
    return data.prop_of_k.mean(), data.prop_of_k.std()

def get_num_fixed(data):
    patches = data['patch_num'].unique()
    s = 0
    for patch in patches:
        s += data.query('patch_num == ' + str(patch) + ' & frequency == 1.0').shape[0]

    return float(s)/float(len(patches))


def write_demography(run_id, source_dir_path, target_dir_path):

    target_file_path = target_dir_path + '/' + demography_path

    gld = source_dir_path + '/' + global_ld_path
    lld = source_dir_path + '/' + local_ld_path
    fst = source_dir_path + '/' + fst_path
    alf = source_dir_path + '/' + allele_freq_path

    global_ld_df = pandas.read_csv(gld)
    local_ld_df = pandas.read_csv(lld)
    fst_df = pandas.read_csv(fst)
    allele_freq_df =  pandas.read_csv(alf)

    gens = allele_freq_df['generation'].unique()

    types = ['fitness', 'neutral']

    for gen in gens:
        for t in types:

            q = 'generation == ' + str(gen) + ' & type == \"' + str(t) + '\"'

            this_gen_alf = allele_freq_df.query(q)
            this_gen_fst = fst_df.query(q)
            this_gen_lld = local_ld_df.query(q)
            this_gen_gld = global_ld_df.query(q)


            fst_mean, fst_sigma = get_fst(this_gen_fst)
            lld_mean, lld_sigma = get_ld(this_gen_lld)
            gld_mean, gld_sigma = get_ld(this_gen_gld)
            n_fixed = get_num_fixed(this_gen_alf)

            fit = 0
            neut = 0
            if t == "fitness":
                fit = 1
            if t == "neutral":
                neut = 1

            write_line = "%d,%d,%f,%f,%f,%f,%f,%f,%f,%d,%d\n" % (run_id, gen, fst_mean, fst_sigma, gld_mean, gld_sigma, lld_mean, lld_sigma, n_fixed, neut, fit)

            with open(target_file_path, 'a') as file:
                file.write(write_line)

def write_dynamics(run_id, source_dir_path, target_dir_path):
    target_file_path = target_dir_path + '/' + dynamics_path

    att_mig = source_dir_path + '/' + att_mig_path
    eff_mig = source_dir_path + '/' + eff_mig_path
    pat = source_dir_path + '/' + patches_path

    att_mig_data = pandas.read_csv(att_mig)
    eff_mig_data = pandas.read_csv(eff_mig)
    patch_data = pandas.read_csv(pat)

    gens = att_mig_data['generation'].unique()

    for gen in gens:
        q = 'generation == ' + str(gen)

        this_gen_att = att_mig_data.query(q)
        this_gen_eff = eff_mig_data.query(q)
        this_gen_patches = patch_data.query(q)

        att_mig_mean, att_mig_sigma = get_att_mig(this_gen_att)
        eff_mig_mean, eff_mig_sigma = get_eff_mig(this_gen_eff)

        prop_extinct = get_prop_extinct(this_gen_patches)
        prop_k_mean, prop_k_sigma = get_prop_k(this_gen_patches)

        write_line = '%d,%d,%f,%f,%f,%f,%f,%f,%f\n' % (run_id, gen, eff_mig_mean, eff_mig_sigma, att_mig_mean, att_mig_sigma, prop_extinct, prop_k_mean, prop_k_sigma)

        with open(target_file_path, 'a') as file:
            file.write(write_line)

def main():
    folder = get_folder()
    init_files(folder)

    id_ct = 0

    for x in (os.listdir(folder)):
        path = os.path.abspath(folder) + '/' + x
        if os.path.isdir(path):
            write_metadata(id_ct, path, folder)
            write_demography(id_ct, path, folder)
            write_dynamics(id_ct, path, folder)
            id_ct += 1


if __name__ == '__main__':
    main()
