dir_genepredicted_infile=$2
dir_geneannotated_outfile=$3
for name in `cat $1`
do
echo $name
#diamond blastp --db /data2/qguo/reference_database/CARD_diamonddb/CARD.dmnd --query $dir_genepredicted_infile/$name'.gene.csqa' --out $dir_geneannotated_outfile/$name'_best.tab'  --outfmt 6  --max-target-seqs 1 --evalue 1e-3  --id 30 --tmpdir ../tmp/
#diamond blastp --db /data2/qguo/reference_databases/card/diamond_db/protein_fasta_protein_homolog_model.dmnd --ungapped-score 95 --min-score 95 --query $dir_genepredicted_infile/$name/$name'.faa' --out $dir_geneannotated_outfile/$name'_best.tab' -p 40 --query-cover 60 --subject-cover 60 --outfmt 6  --max-target-seqs 1 --evalue 1e-3  --id 50
diamond blastp --db /data2/qguo/reference_databases/card/card.dmnd --ungapped-score 95 --min-score 95 --query $dir_genepredicted_infile/$name/$name'.faa' --out $dir_geneannotated_outfile/$name'_best.tab' -p 40 --query-cover 80 --subject-cover 80 --outfmt 6  --max-target-seqs 1 --evalue 1e-3  --id 80
done
