import pysam
import sys

def grepUsingReads(in_bam, out_bam):
    bf = pysam.AlignmentFile(in_bam, "rb", check_sq = False)
    bf_head_dict = dict(bf.header)

    with pysam.AlignmentFile(out_bam, "wb", header = bf_head_dict) as outf:
        for r in bf:
            readType = r.get_tag("XF")
            if readType in "CODING" or readType in "INTRONIC":
                outf.write(r)
        outf.close()

if __name__ == "__main__":
    in_bam = sys.argv[1]
    out_bam = sys.argv[2]

    grepUsingReads(in_bam, out_bam)