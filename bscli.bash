#!bash
PROJECTNAME=$1


mkdir -p $1/fastq

module load basespace-cli

bs download project --name $1 -o $1/fastq --extension fastq.gz


for f in $PROJECTNAME/fastq/*/*.gz; do
    folder_name=${f#*/fastq/}
    folder_name=${folder_name%%/*}
    s=${folder_name%%_*}
    mv $f $PROJECTNAME/fastq/$s.${f##*/}
done

cd $1

module load mothur/1.39.5
mothur "#make.file(inputdir=fastq, type=gz);"
