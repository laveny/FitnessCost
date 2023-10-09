from sys import argv
import os
import regex

fq1, fq2, output_dir, output_base = argv[1:]

output = output_dir + output_base + ".txt"
out1=open(output, "w")

f_re = "(?e)(?P<p1>CTGAGTAGGACAAATCCGCC){e<=1}(?P<umi_f1>[ATCG]{15})(?P<p2>ATATGAGGCTTATCGTGAAG){e<=1}(?P<umi_f2>[ATCG]{15})(?P<p2>ACGAGTGTTACTGCAGCTG){e<=1}"
r_re = "(?e)(?P<p3>CAGCTGCAGTAACACTCGT){e<=1}(?P<umi_r1>[ATCG]{15})(?P<p2>CTTCACGATAAGCCTCATAT){e<=1}(?P<umi_r2>[ATCG]{15})(?P<p2>GGCGGATTTGTCCTACTCAG){e<=1}"
#new_f_re = "(?e)(?P<p1>ACAGTAACGGGATCG){e<=0}(?P<umi_f5>[ATCG]{7,12})(?P<p2>ACGAGTGTTACTGCAGCTG){e<=0}"
#new_r_re = "(?e)(?P<p3>CAGTAACACTCGT){e<=0}(?P<umi_r5>[ATCG]{7,12})(?P<p2>CGATCCCGTTACTGT){e<=0}"
new_f_re = "(?e)(?P<p1>AGTAACGGGATCGACGAG){e<=1}(?P<umi_f5>[ATCG]{9,12})(?P<p2>TGTTACTGCAGCTG){e<=1}"
new_r_re = "(?e)(?P<p3>CAGCTGCAGTAACA){e<=1}(?P<umi_r5>[ATCG]{9,12})(?P<p2>CTCGTCGATCCCGTTACT){e<=1}"

nofunc_umi = ['CGATAGTTA','CCGTTCGTAAG','TTTCCATGAGCA']

with open(fq1, "r") as f1, open(fq2, "r") as f2:
    for j, element in enumerate(zip(f1, f2)):
        if j % 4 == 1:
            forward_seq = regex.search(f_re, element[0].strip(), flags=0)
            reverse_seq = regex.search(r_re, element[1].strip(), flags=0)
            if forward_seq is not None and reverse_seq is not None:
                forward_umi = forward_seq.groupdict()["umi_f1"] + forward_seq.groupdict()["umi_f2"]
                reverse_umi = reverse_seq.groupdict()["umi_r1"] + reverse_seq.groupdict()["umi_r2"]
                reverse_umi = reverse_umi[::-1].replace('A', 't').replace('T', 'a').replace('G', "c").replace("C",'g').upper()
                if forward_umi == reverse_umi:
                    out1.write(forward_umi + "\n")
            else:
                forward_kongzai = regex.search(new_f_re, element[0].strip(), flags=0)
                reverse_kongzai = regex.search(new_r_re, element[1].strip(), flags=0)
                if forward_kongzai is not None and reverse_kongzai is not None:
                    kongzai_forward_umi = forward_kongzai.groupdict()["umi_f5"]
                    kongzai_reverse_umi = reverse_kongzai.groupdict()["umi_r5"]
                    kongzai_reverse_umi = kongzai_reverse_umi[::-1].replace('A', 't').replace('T', 'a').replace('G',"c").replace("C", 'g').upper()
                    if kongzai_forward_umi == kongzai_reverse_umi:
                        out1.write(kongzai_forward_umi + "\n")

out1.close()
