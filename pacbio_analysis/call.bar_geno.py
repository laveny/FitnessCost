from sys import argv
import regex as re

name_prefix = argv[1]
mud_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
Reverse_complemrnt = lambda x: "".join([mud_dict[i] for i in x][::-1])

geno_pattern = re.compile("(GTTTTTTTGAGTGTTTCTC){e<=1}\w{1626}(GTAAAACGACGGCCAGTCAA){e<=1}")
umi_pattern = re.compile("CTGAGTAGGACAAATCCGCC(?P<Barcode>[ATCG]+ATATGAGGCTTATCGTGAAG[ATCG]+)ACGAGTGTTACTGCAGCTG")
outfile = open(name_prefix+".call.bar_geno.txt", "w")
with open(name_prefix + ".n.ccs.sam", "r") as inputr:  # nagetive strand ccs result
    for nega in inputr:
        if nega[0] == "m":
            nega_line_list = nega.strip().split()
            nega_seq = nega_line_list[9]  # the result sequence after ccs
            zmw_nega = nega_line_list[0].split("/")[1]  # the number of ZMW
            nega_reverse = Reverse_complemrnt(nega_seq)
            geno2 = re.search(geno_pattern, nega_reverse, flags=0)
            umi2 = re.search(umi_pattern, nega_reverse, flags=0)
            if geno2 is not None and umi2 is not None:
                geno3 = geno2.group()[19:1647]
                umi3 = umi2.groupdict()["Barcode"]
                outfile.write(">" + str(zmw_nega) + "n" + "\n" + "n_" + geno3 + " " + umi3 +  "\n")
with open(name_prefix + ".p.ccs.sam", "r") as inputf:  # positive strand ccs result
    for posi in inputf:
        if posi[0] == "m":
            posi_line_list = posi.strip().split()
            posi_seq = posi_line_list[9]
            zmw_posi = posi_line_list[0].split("/")[1]
            geno = re.search(geno_pattern, posi_seq, flags=0)
            umi = re.search(umi_pattern, posi_seq, flags=0)
            if geno is not None and umi is not None:
                geno1 = geno.group()[19:1647]
                umi1 = umi.groupdict()["Barcode"]
                outfile.write(">" + str(zmw_posi) + "p" + "\n" + "p_" + geno1 + " " + umi1 +  "\n")

outfile.close()
