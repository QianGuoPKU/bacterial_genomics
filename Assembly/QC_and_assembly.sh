#inputfile decompressed file
dir_fq_infile=$1
dir_prinseq_outfile=$2
for name in `cat $3`
do
echo "$name"
sed -i '/@/ s/$/\/1/' $dir_fq_infile/$name'_1.fq'
sed -i '/@/ s/$/\/2/' $dir_fq_infile/$name'_2.fq'
sed -i '/@/ s/ /#/g' $dir_fq_infile/$name'_1.fq'
sed -i '/@/ s/ /#/g' $dir_fq_infile/$name'_2.fq'
sed -i '/@/ s/#2/#1/g' $dir_fq_infile/$name'_2.fq'
#tick out bad sequences
/lustre1/hqzhu_pkuhpc/mli/bin/prinseq-lite -verbose -fastq $dir_fq_infile/$name'_1.fq' -ns_max_p 10 -min_qual_mean 25 -out_good $dir_prinseq_outfile/$name'_1_good' -out_bad $dir_prinseq_outfile/$name'_1_bad'
/lustre1/hqzhu_pkuhpc/mli/bin/prinseq-lite -verbose -fastq $dir_fq_infile/$name'_2.fq' -ns_max_p 10 -min_qual_mean 25 -out_good $dir_prinseq_outfile/$name'_2_good' -out_bad $dir_prinseq_outfile/$name'_2_bad'
manageread extractpair pairfastq=$2/$name'_1_good.fastq',$2/$name'_2_good.fastq' outprefix=$4/$name
spades.py -1 $4/$name'_1.fastq' -2 $4/$name'_2.fastq'  --phred-offset 33 -o $5/$name -t 30 -m 1024
done

