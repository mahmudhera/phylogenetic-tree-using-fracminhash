# This script will compare the handful of ways we can construct a distance matrix/tree
# using ANI & containment

# This command imports all the functions in compare_distances so I don't need to copy over code
from compare_distances import *
import subprocess
import os
import numpy as np
import sourmash as sm

out_dir = 'comparing_matrices'
os.makedirs(out_dir, exist_ok=True)

###########
# Mahmudur code to generate the matrices. Note: this uses the construct_fmh_sketches.py
genome_list_filename = 'genome-list-bacteria'
sketch_directory = 'fmh_sketches_bacteria'
k = 21
scale_factor = 0.1
seed = 0
dist_matrix_filename = 'pairwise_dist_matrix'

list_genomes_sketches = []

# read sketches
genome_list = read_genome_list(genome_list_filename)
for (gname, gpath) in genome_list:
    sketch_filename = sketch_directory + '/fmh_sketch_k_' + str(k) + '_scale_f_' + str(scale_factor) + '_seed_' + str(seed) + '_genome_' + gname
    fmh = read_sourmash_sketch(sketch_filename, scale_factor)
    list_genomes_sketches.append( (gname, fmh) )

# construct the ANI, distance, and (average) containment matrices
mrh_ANI_matrix = []
mrh_distance_matrix = []
mrh_containment_matrix = []
for i in range(len(list_genomes_sketches)):
    distance_list = []
    containment_list = []
    ANI_list = []
    for j in range(len(list_genomes_sketches)):
        if i == j:
            distance = 0.0
            containment = 1.0
            ANI = 1.0
        else:
            fmh1 = list_genomes_sketches[i][1]
            fmh2 = list_genomes_sketches[j][1]
            d1 = containment_to_mutation_rate(fmh1.get_containment(fmh2), k)
            d2 = containment_to_mutation_rate(fmh2.get_containment(fmh1), k)
            c1 = fmh1.get_containment(fmh2)
            c2 = fmh2.get_containment(fmh1)
            containment = (c1 + c2) / 2.0
            distance = (d1 + d2) / 2.0
            ANI = 1 - distance
        distance_list.append(distance)
        containment_list.append(containment)
        ANI_list.append(ANI)
    mrh_distance_matrix.append(distance_list)
    mrh_containment_matrix.append(containment_list)
    mrh_ANI_matrix.append(ANI_list)
mrh_distance_matrix = np.array(mrh_distance_matrix)
mrh_containment_matrix = np.array(mrh_containment_matrix)
mrh_ANI_matrix = np.array(mrh_ANI_matrix)
np.set_printoptions(precision=4)
print(f"mrh ANI matrix: {mrh_ANI_matrix}")
# export to numpy
for mat_type in ["distance", "ANI", "containment"]:
    save_name = f"mrh_{mat_type}_matrix"
    np.save(os.path.join(out_dir, save_name), vars()['mrh_' + mat_type + '_matrix'])
    # save the file names
    with open(f"{save_name}.labels.txt", 'w') as fid:
        for name, _ in genome_list:
            fid.write(name)

####################################
# Use the built-in sourmash approach
# Remove the sketches if they are already there
sketch_dir = os.path.join(out_dir, 'sourmash_sketches_bacteria')
os.makedirs(sketch_dir, exist_ok=True)
# Build the sketches
sketch_command = f"sourmash compute -k 21 --scale 10 --output-dir {sketch_dir} "
for _, gpath in genome_list:
    sketch_command += gpath + ' '
subprocess.call(sketch_command, shell=True)
# Then run sourmash compare with ANI
save_prefix = 'sourmash_compare_ANI'
compare_ANI_command = f"sourmash compare -k 21 --ani -o {save_prefix}.cmp "
for _, gpath in genome_list:
    sig_path = os.path.join(sketch_dir, gpath.split('/')[-1]) + '.sig'
    compare_ANI_command += sig_path + ' '
subprocess.call(compare_ANI_command, shell=True)
sm_compare_ANI_matrix = np.load(f"{save_prefix}.cmp")
# Then run sourmash compare with containment
save_prefix = 'sourmash_compare_containment'
compare_containment_command = f"sourmash compare -k 21 --containment -o {save_prefix}.cmp "
for _, gpath in genome_list:
    sig_path = os.path.join(sketch_dir, gpath.split('/')[-1]) + '.sig'
    compare_containment_command += sig_path + ' '
subprocess.call(compare_containment_command, shell=True)
sm_compare_containment_matrix = np.load(f"{save_prefix}.cmp")


####################################
# Use the ANI estimate from https://github.com/sourmash-bio/sourmash/blob/851dc2b1dc215ef59f0b5a81d67900329690aa91/tests/test_signature.py#L431
tessa_dist_matrix = []
tessa_ANI_matrix = []
tessa_containment_matrix = []
for _, gpath1 in genome_list:
    ANI_list = []
    dist_list = []
    cont_list = []
    for _, gpath2 in genome_list:
        sig_path1 = os.path.join(sketch_dir, gpath1.split('/')[-1]) + '.sig'
        sig_path2 = os.path.join(sketch_dir, gpath2.split('/')[-1]) + '.sig'
        ss1 = sm.load_one_signature(sig_path1, ksize=21)
        ss2 = sm.load_one_signature(sig_path2, ksize=21)
        s1_cont_s2 = ss1.containment_ani(ss2, estimate_ci=True)
        s2_cont_s1 = ss2.containment_ani(ss1, estimate_ci=True)
        dist = (s1_cont_s2.dist + s2_cont_s1.dist) / 2.0
        ANI = 1 - dist
        ANI_list.append(ANI)
        dist_list.append(dist)
        cont_list.append(ss1.avg_containment(ss2))
    tessa_dist_matrix.append(dist_list)
    tessa_ANI_matrix.append(ANI_list)
    tessa_containment_matrix.append(cont_list)
tessa_dist_matrix = np.array(tessa_dist_matrix)
tessa_ANI_matrix = np.array(tessa_ANI_matrix)
tessa_containment_matrix = np.array(tessa_containment_matrix)
np.save(os.path.join(out_dir, "tessa_ANI_matrix"), tessa_ANI_matrix)
np.save(os.path.join(out_dir, "tessa_dist_matrix"), tessa_dist_matrix)
np.save(os.path.join(out_dir, "tessa_containment_matrix"), tessa_containment_matrix)

################
# Then compare everything
print(f"Mahmudur ANI: {mrh_ANI_matrix}")
print(f"Sourmash compare ANI: {sm_compare_ANI_matrix}")
print(f"Tessa ANI: {tessa_ANI_matrix}")
