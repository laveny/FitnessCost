from sys import argv

name_prefix, in_path, out_path = argv[1:]


print('%s start' % name_prefix)

todo = []
with open(in_path + name_prefix, 'r') as f:
    for line in f:
        bar = line.strip()
        todo.append(bar)

f.close()

with open(out_path + name_prefix + '.havesim.txt', 'w') as f1:
    cc = 1
    for a in todo:
        total = []
        for i in range(len(a)):
            for m in ['A', 'T', 'C', 'G']:
                st = a[:i] + m + a[i + 1:]
                total.append(st)
        total = list(set(total))
        total.sort()
        for g in total:
            f1.write('%s\t%s\n' % (g,a))
        cc = cc + 1
        if cc == 500:
            f1.flush()
            cc = 1
f1.close()

print('%s success' % name_prefix)
