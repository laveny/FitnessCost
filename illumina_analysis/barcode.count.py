import collections
import os
import json
from sys import argv

input, output_dir, output_base = argv[1:]

count_dict = {}

with open(input, "r") as f1:
    str1 = f1.read().split("\n")
    for k in str1:
        if k in count_dict:
            count_dict[k] += 1
        else:
            count_dict[k] = 1

output = output_dir + output_base + ".count.txt"


with open(output, "w") as out1:
    for i in count_dict:
        out1.write('%s\t%s\n' % (i, count_dict[i]))

out1.close()


output = output_dir + output_base + ".json"

json.dump(count_dict , open(output,"w"))