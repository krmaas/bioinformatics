



#!bash
PROJECTNAME=$1

module load basemount


basemount ../basemountpoint/basespace

chmod 777 /yOUReNTIREdIRECTORY

mkdir -p $1/fastq

        for f in ../basemountpoint/basespace/Projects/$PROJECTNAME/Samples/*/Files/*.gz;
        do s=${f##../basemountpoint/basespace/Projects/$PROJECTNAME/Samples/}; 
        s=${s%%/*};
        cp  "$f" "$PROJECTNAME"/fastq/"$s"."${f##*Files/}" ;
        done

cd $1

module load mothur
mothur "#make.file(inputdir=fastq, type=gz);"

basemount --unmount /home/CAM/kmaas/basemountpoint/basespace
chmod 700 /yOUReNTIREdIRECTORY

## I have reference files and mothur/slurm files in one folder, so I can easily populate the folder with exactly the files I'll use
#rsync ../../mothur.files/* .

mv fastq/stability.files $1".file"
