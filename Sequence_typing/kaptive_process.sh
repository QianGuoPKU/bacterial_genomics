for name in `cat $1`
do
echo $name
kaptive.py -a $2/$name/$name'.fna' -k /data2/qguo/software/Kaptive-master/reference_database/Acinetobacter_baumannii_k_locus_primary_reference.gbk -o $3/$name
kaptive.py -a $2/$name/$name'.fna' -k /data2/qguo/software/Kaptive-master/reference_database/Acinetobacter_baumannii_OC_locus_primary_reference.gbk -o $4/$name
done



