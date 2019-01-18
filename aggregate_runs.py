#! /usr/bin/env python

import csv, argparse, os, pandas
import networkx as nx

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
        file.write(dynamics_csv_header)
    with open(meta, 'a') as file:
        file.write('run_id,')

        p = os.path.abspath(folder) + '/' + os.listdir(folder)[0] + '/params.ini'
        print 'okay'
        df = pandas.read_csv(p, sep=' ', error_bad_lines=False, header=None)
        for row in df.itertuples():
            file.write(str(row[1]) + ',')
        file.write('\n')


metadata = pandas.DataFrame()

def add_metadata(run_id, source_dir_path, target_dir_path):
    global metadata
    param_file_path = source_dir_path + '/' + params_path
    target_file_path = target_dir_path + '/' + metadata_path

    df = pandas.read_csv(param_file_path, error_bad_lines=False, sep=',', header=None)
    df = df.transpose()
    df.columns = df.iloc[0]
    df = df.reindex(df.index.drop(0))
    df['run_id'] = run_id

    print df
    if len(metadata) == 0:
        metadata = df
    else:
        metadata = pandas.concat([metadata, df], sort=True)
    print metadata
    #with open(target_file_path, 'w') as f:
        #metadata.to_csv(f, header=True)


def write_metadata(target_dir_path):
    global metadata
    target_file_path = target_dir_path + '/' + metadata_path

    with open(target_file_path, 'a') as f:
        metadata.to_csv(f, header=True)

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

def get_graph_stats(data):
    gen0 = data.query('generation == 0')

    def f(x):
        return np.exp(-1*x^2)

    patches = []

    for index,row in gen0.iterrows():
        x,y = row['x'], row['y']
        patches.append((x,y))

    n_patches = len(patches)

    mat = np.zeros((n_patches, n_patches))

    for i, (x,y) in enumerate(patches):
        for j, (x2,y2) in enumerate(patches):
            mat[i,j] = np.sqrt((x2-x)*(x2-x) + (y2-y)*(y2-y))

    incidence_matrix = np.zeros((n_patches,n_patches))


    for i in range(n_patches):
        for j in range(n_patches):
            incidence_matrix[i,j] = f(mat[i,j])


    dispersal_probs = np.zeros((n_patches, n_patches))
    G = nx.Graph()


    for i in range(n_patches):
        row_sum = 0
        for j in range(n_patches):
            row_sum += incidence_matrix[i,j]
        for j in range(n_patches):
            dispersal_probs[i,j] = incidence_matrix[i,j] / row_sum
            G.add_edge(i,j, w=dispersal_probs[i,j])

    centrality =  nx.eigenvector_centrality(G, weight='w')


def write_demography(run_id, source_dir_path, target_dir_path):
    try:
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
    except:
        print('failed w/ ' + source_dir_path)
def write_dynamics(run_id, source_dir_path, target_dir_path):
    try:
        target_file_path = target_dir_path + '/' + dynamics_path

        att_mig = source_dir_path + '/' + att_mig_path
        eff_mig = source_dir_path + '/' + eff_mig_path
        pat = source_dir_path + '/' + patches_path

        att_mig_data = pandas.read_csv(att_mig)
        eff_mig_data = pandas.read_csv(eff_mig)
        patch_data = pandas.read_csv(pat)

        gens = eff_mig_data['generation'].unique()

        for gen in gens:
            q = 'generation == ' + str(gen)

            this_gen_att = att_mig_data.query(q)
            this_gen_eff = eff_mig_data.query(q)
            this_gen_patches = patch_data.query(q)

            #att_mig_mean, att_mig_sigma = get_att_mig(this_gen_att)
            att_mig_mean, att_mig_sigma = 0,0
            eff_mig_mean, eff_mig_sigma = get_eff_mig(this_gen_eff)

            prop_extinct = get_prop_extinct(this_gen_patches)
            prop_k_mean, prop_k_sigma = get_prop_k(this_gen_patches)

            write_line = '%d,%d,%f,%f,%f,%f,%f,%f,%f\n' % (run_id, gen, eff_mig_mean, eff_mig_sigma, att_mig_mean, att_mig_sigma, prop_extinct, prop_k_mean, prop_k_sigma)

            with open(target_file_path, 'a') as file:
                file.write(write_line)
    except:
        print('failed w/ ' + source_dir_path)

def main():
    folder = get_folder()
    init_files(folder)

    id_ct = 0

    for x in (os.listdir(folder)):
        path = os.path.abspath(folder) + '/' + x
        if os.path.isdir(path):
            print path
            add_metadata(id_ct, path, folder)
            write_demography(id_ct, path, folder)
            write_dynamics(id_ct, path, folder)
            id_ct += 1

    write_metadata(folder)

if __name__ == '__main__':
    main()
