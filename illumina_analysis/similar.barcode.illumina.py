from sys import argv
from fuzzysearch import find_near_matches
import pandas as pd
import regex as re

name_prefix, in_path, out_path = argv[1:]

fn = pd.read_csv("/mnt/data3/disk/guoziyan/Pacbio/V2max/similar/output230516/final.result/final_bar_geno_230517.txt",
                 names=['barcode','geno'],sep="\t")

pacbio = list(fn['barcode'])

todo = []
with open(in_path + name_prefix, 'r') as f:
    for line in f:
        bar = line.strip()
        todo.append(bar)

f.close()

with open(out_path + name_prefix + '.havesim.txt', 'w') as f1, open(out_path + name_prefix + '.morethanonesim.txt',
                                                                    'w') as f2, open(
    out_path + name_prefix + '.nosim.txt', 'w') as f3:
    co = 1
    for u in todo:
        pp = '(%s){s<=1}' % (u)
        main = u  ## for multiple mapping
        mat = []
        for i in pacbio:
            result = re.fullmatch(pp, i)
            if result is not None:
                mat.append(i)
        if len(set(mat)) == 1:  ## either u is in pacbio or u has one sim in pacbio
            f1.write('%s\t%s\n' % (u,mat[0]))
        elif len(set(mat)) == 0:
            f3.write(u + ':' + '\t'.join(mat) + '\n')
        else:
            f2.write('%s\t%s\n' % (u, mat[0]))
        co = co + 1
        if co == 500:
            f1.flush()
            f2.flush()
            f3.flush()
            co = 1

f1.close()
f2.close()
f3.close()
print('%s success' % (name_prefix))