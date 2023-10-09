import numpy as np
import pandas as pd
import regex as re
import json
import os

input_dir = "./Native/barcode.count/"
a = os.listdir(input_dir)
file_list=[]
for k in a:
    if 'json' in k:
        k = os.path.join(input_dir,k)
        file_list.append(k)
barcode = []
for i in file_list:
    dic = json.load(open(i,'r'))
    barcode = barcode + list(dic.keys())
    barcode = list(set(barcode))
    print(i)

input_dir = "./Inducible/barcode.count/"
a = os.listdir(input_dir)
file_list=[]
for k in a:
    if 'json' in k:
        k = os.path.join(input_dir,k)
        file_list.append(k)
barcode_i = []
for i in file_list:
    dic = json.load(open(i,'r'))
    barcode_i = barcode_i + list(dic.keys())
    barcode_i = list(set(barcode_i))
    print(i)

todo = list(set(barcode+barcode_i))
todo = todo[1:]
with open('illumina.barcode.txt','w') as f1:
    for i in todo:
        f1.write('%s\n' % i )
