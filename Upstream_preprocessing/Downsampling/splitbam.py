##### Code has not been tested on unsorted bam files, sort on barcode (CB):
##### samtools sort -t CB unsorted.bam  -o sorted_tags.bam
###
##### INPUT: .bam file to be sorted and output directory to place split BC
##### OUTPUT: .bam file for each unique barcode, best to make a new directory

### Python 3.6.8
import pysam
import sys
def splitbarcode(unsplit_file, out_dir):

# variable to hold barcode index
        CB_hold='unset'
        itr=0
# read in upsplit file and loop reads by line
        samfile=pysam.AlignmentFile( unsplit_file, "rb")
        for read in samfile.fetch( until_eof=True):
                        # barcode itr for current read
                        CB_itr = read.get_tag( 'XC')
# if change in barcode or first line; open new file
                        if( CB_itr!=CB_hold or itr==0):
                                        # close previous split file, only if not first read in file
                                        if( itr!=0):
                                                        split_file.close()
                                        CB_hold = CB_itr
                                        itr+=1
                                        split_file=pysam.AlignmentFile( out_dir + "XC_{}.bam".format( itr), "wb", template=samfile)

    # write read with same barcode to file
                        split_file.write( read)
        split_file.close()
        samfile.close()


if __name__ == "__main__":
        unsplit_file=sys.argv[1]
        out_dir=sys.argv[2]
        splitbarcode(unsplit_file, out_dir)
