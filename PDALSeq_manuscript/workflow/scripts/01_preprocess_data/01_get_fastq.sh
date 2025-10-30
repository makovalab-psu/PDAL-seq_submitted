#Get P5 Illumina TruSeq reads from SRA
#Args $1 = the SRA extension
#Args $2 = The output file in fastq.gz format

echo $1
ACN=$(echo $1 | awk '{split($0,a,"_TruSeq_"); print a[2]}' | awk '{gsub("_forward","",$1); print}')
echo $ACN

if [ ! -f $2 ]; then
    echo "File not found!"

    fastq-dump --accession $ACN --split-files --split-spot --gzip --outdir temp/

    echo "temp/"$ACN"_1.fastq.gz" 
    echo "temp/"$ACN"*"
    mv "temp/"$ACN"_1.fastq.gz" $2
    rm temp/$ACN*

fi


