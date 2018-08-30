#! /bin/bash


input_file=$1

echo $input_file
thisdir=$(pwd)
output_dir=$(pwd)/out
prokka_dir=${output_dir}/prokka_output
prophet_dir=${output_dir}/ProphET_output
mobsuite_dir=${output_dir}/mobsuite_output
prefix="${RANDOM}${RANDOM}"
logfile=${output_dir}/mobilephone.log

mkdir $output_dir


echo $output_dir
echo $prokka_dir
echo $prophet_dir
echo $prefix

echo "test" >> $logfile
echo "Running Prokka"
prokka $input_file -o ${prokka_dir} --prefix $prefix 2>> $logfile

prokka_fna=${prokka_dir}/${prefix}.fna
prokka_gff=${prokka_dir}/${prefix}.gff
prokka_new_gff=${prokka_dir}/new.gff

echo "Reformatting gff"
cd submodules/ProphET/UTILS.dir/GFFLib/
perl ./gff_rewrite.pl --input  $prokka_gff --output $prokka_new_gff --add_missing_features 2>> $logfile
cd $thisdir
echo "Running ProphET"
./submodules/ProphET/ProphET_standalone.pl --fasta_in $prokka_fna --gff_in $prokka_new_gff --outdir ${prophet_dir} 2>> $logfile


echo "Running MOBSUITE to detect and type plasmids"
mob_recon --infile $input_file --outdir $mobsuite_output --run_typer
