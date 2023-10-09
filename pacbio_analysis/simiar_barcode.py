from sys import argv
from fuzzysearch import find_near_matches
import pandas as pd
import regex as re

name_prefix, in_path, out_path = argv[1:]

count = []
barcode = []
geno = []
with open("/mnt/data3/disk/guoziyan/Pacbio/V2max/filtered.bar_geno.txt", 'r') as f:
    for line in f:
        line_sp = line.strip().split('\t')
        geno.append(line_sp[1])
        line_sp = line_sp[0]
        line_sp = line_sp.strip().split(' ')
        count.append(line_sp[0])
        barcode.append(line_sp[1])
f.close()
fn = pd.DataFrame({'barcode': barcode, 'geno': geno, 'count': count})
bar_uni = list(set(barcode))
todo = []
with open(in_path + name_prefix, 'r') as f:
    for line in f:
        bar = line.strip()
        todo.append(bar)

f.close()

with open(out_path + name_prefix + '.sim_bar.txt', 'w') as f1, open(out_path + name_prefix + '.filtered_bar_geno.txt',
                                                                    'w') as f2, open(
    out_path + name_prefix + '.nothavemajority.txt', 'w') as f3, open(out_path + name_prefix + '.unsolved.txt',
                                                                      'w') as f4:
    for u in todo:
        pp = '(%s){s<=1}' % (u)
        main = u  ## for multiple mapping
        mat = []
        conti = False
        for i in bar_uni:
            result = re.fullmatch(pp, i)
            if result is not None:
                mat.append(i)
        if len(set(mat)) == 1:  ## if u is uniq
            if mat[0] == u:
                conti = True
            else:
                print('%s uniq wrong\n' %(u))
                f4.write('%s\n' % (u))
        elif len(set(mat)) > 1:  ## if u is multiple similar barcode
            f1.write(u + ':' + '\t'.join(mat) + '\n')  ## record similar barcode clusters
            ## define several cluster structures: only main barcode were continued
            cluster_dict = {}
            for eli in mat:
                mm = []
                pp = '(%s){s<=1}' % (eli)
                for i in bar_uni:
                    result = re.fullmatch(pp, i)
                    if result is not None:
                        mm.append(i)
                cluster_dict[eli] = mm
                if len(mm) > len(set(mat)):
                    main = eli
            if main == u:  ## which means that either the barcode is uniq that has no similar barcode, or u is the main
                ## to check if the main has no two-step further barcode
                total = []
                for eliment in cluster_dict:
                    total = total + cluster_dict[eliment]
                if set(total) > set(
                        mat):  ## if the main has two-step far barcode, not continue and write to unsolved.txt
                    conti = False
                    f4.write('%s\n' % (u))
                else:
                    conti = True
        if conti == True:
            sub = fn[fn['barcode'].isin(mat)]
            sub = sub.astype({'count': 'int64'})
            if sub.dtypes['count'] != 'int64':
                print('%s type wrong\n' % (name_prefix))
            ss = sub.groupby(['geno']).sum('count')
            ss = ss.sort_values('count', ascending=False)
            if ss['count'].iloc[0] / sum(ss['count']) > 0.5:
                sub = sub[sub['geno'] == ss.index[0]]
                sub = sub.sort_values('count',ascending = False)
                if len(sub) == 1:
                    f2.write('%s\t%s\t%s\n' % (sub['barcode'].iloc[0], sub['geno'].iloc[0], sub['count'].iloc[0]))
                elif len(sub) > 1:
                    if  sub['count'].iloc[0] > sub['count'].iloc[1]:
                        f2.write('%s\t%s\t%s\n' % (sub['barcode'].iloc[0], sub['geno'].iloc[0], sub['count'].iloc[0]))
                    else:
                        print('%s majority01 wrong\n' % (name_prefix))                                                
            else:
                if ss['count'].iloc[0] < ss['count'].iloc[1]:
                    print('%s majority02 wrong' % (name_prefix))
                else:
                    f3.write('%s\n' % (u))
print('%s success' % (name_prefix))
