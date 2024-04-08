for name in `cat $1`
do
echo $name
prokka --compliant --kingdom Bacteria --genus Acinetobacter --species baumannii --cpus 40 --outdir $name --prefix $name $2/$name/scaffolds.fasta
done



