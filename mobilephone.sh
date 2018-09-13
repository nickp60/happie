#! /bin/bash
set -euo pipefail
IFS=$'\n\t'
readonly LOG_FILE="/tmp/$(basename "$0").log"
info()    { echo "[INFO]    $*" | tee -a "$LOG_FILE" >&2 ; }
warning() { echo "[WARNING] $*" | tee -a "$LOG_FILE" >&2 ; }
error()   { echo "[ERROR]   $*" | tee -a "$LOG_FILE" >&2 ; }
fatal()   { echo "[FATAL]   $*" | tee -a "$LOG_FILE" >&2 ; exit 1 ; }
echo "Logging to $LOG_FILE"
#/ Usage: extract mobile regions from <fasta file> uing <cores> cores
usage() {
    grep '^#/' "$0" | cut -c4-
    exit 0
}
expr "$*" : ".*--help" > /dev/null && usage



input_file=$1
cpus=$2


echo $input_file
echo $cpus
main_output_dir=`pwd`/out

thisdir=$(pwd)
scripts_dir=$(pwd)/scripts
logfile=${main_output_dir}/mobilephone.log
mainprefix="${RANDOM}${RANDOM}"


mkdir $main_output_dir

main(){
    file=$1
    counter=$2
    
    prefix=$mainprefix$counter
    output_dir=${main_output_dir}/${counter}/
    
    prokka_dir=${output_dir}/prokka_output
    prophet_dir=${output_dir}/ProphET_output
    mobsuite_dir=${output_dir}/mobsuite_output
    mlplasmids_dir=${output_dir}/mlplasmids_output
    island_dir=${output_dir}/islandpath_dimob_output



    echo "test" >> $logfile
    info  "Running Prokka"
    prokka $file -o ${prokka_dir} --prefix $prefix --fast --cpus $cpus 2>> $logfile

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
    mob_recon --infile $file --outdir $mobsuite_dir --run_typer --keep_tmp 2>> $logfile


    echo "Running mlplasmids to detect and type plasmids"
    mkdir $mlplasmids_dir
    # this takes ages the first time, cause its got to get the model
    Rscript ${scripts_dir}/run_mlplasmids.R $file ${mlplasmids_dir}/mlplasmids_results.tab .8 "Escherichia coli" 2>> $logfile

    # echo "Running islandpath to detect genomic islands"
    # mkdir $island_dir
    # perl ./submodules/islandpath/Dimob.pl ${prokka_dir}/${prefix}.gbk ${island_dir}/output.txt

}


############################ Main ###################################

python $scripts_dir/splitfasta.py $input_file $main_output_dir/split_files/
counter=0
for path in $main_output_dir/split_files/* ; do
    echo $path
    main $path $counter
    counter=$((counter+1))
done


#main
#while read id start stop; do cat
