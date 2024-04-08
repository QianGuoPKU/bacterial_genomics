dir_genepredicted_infile=$2
dir_geneannotated_outfile=$3
for name in `cat $1`
do
echo $name
diamond blastp --db /data2/qguo/reference_databases/BacMet/bacmet_exp.dmnd --ungapped-score 95 --min-score 95 --query $dir_genepredicted_infile/$name/$name'.faa' --out $dir_geneannotated_outfile/$name'_best.tab' -p 40 --query-cover 80 --subject-cover 80 --outfmt 6  --max-target-seqs 1 --evalue 1e-3  --id 80
done
