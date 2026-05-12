# download reprocessing code
wget https://github.com/cellgeni/reprocess_public_10x/archive/0aa3606d57003407f1080777d646977a0ae5e15e.tar.gz
tar -xzvf 0aa3606d57003407f1080777d646977a0ae5e15e.tar.gz
mv reprocess_public_10x-0aa3606d57003407f1080777d646977a0ae5e15e reprocess_public_10x

# run reprocessing
# fill datasets.txt with ids to be reprocessed
while read d
do
 echo $d
 nohup ./reprocess_public_10x/reprocess_public_10x.sh $d > ${d}.log 2>&1 &
done < datasets.txt


###########################################################################################
# some datasets required some manual adjustments, these are explained below ###############
# some submissions has ATAC libraries, these were skipped #################################
###########################################################################################

#############################
## E-MTAB-9536 ##############
#############################

# some samples (ERS23679623 ERS23679624 ERS23679625 ERS23679626 ERS23679627 ERS23679628 ERS23679629 ERS23679630 ERS23679631)
# have R1 of different length (25 and 28) that seems to be result of corresponding libraries being sequenced twice with different settings. 
# most reads (>99%) are 28nt long, so lets just remove all short
# names of reads with shorter R1 start with MS6_31352, so I'll just remove these reads

SERIES=E-MTAB-9536
cd $SERIES
# make list with samples that didn't work, keep full list
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list

cd fastqs
mkdir ../fastqs_orig
while read s
do
    echo $s
    #mv $s ../fastqs_orig/${s}_orig
    #mkdir $s
    for f in `ls ../fastqs_orig/${s}_orig`
    do
        echo "-${f}"
        bsub -n 3 -q normal \
            -R"span[hosts=1] select[mem>1000] rusage[mem=1000]" -M1000 \
            -J fq \
            -o %J.fq.log -e %J.fq.err \
            "zcat ../fastqs_orig/${s}_orig/${f} | awk '/^@MS6_31352/ {skip=3; next} skip>0 {skip--; next} {print}' | gzip > ${s}/${f}"
    done
done < ../${SERIES}.sample.list 

cd ..
# remove falied
while read s
do
 rm -rf $s
done < ${SERIES}.sample.list

# restart mapping
nohup ./run_starsolo.sh $SERIES  > pm19_solo.log 2>&1 &

./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list

#############################
## GSM4088785 ###############
#############################

# some samples (GSM4088774 GSM4088775 GSM4088776 GSM4088777 GSM4088778 GSM4088779 GSM4088780 GSM4088781 GSM4088782 GSM4088783 GSM4088785 GSM4088786 GSM4088787 GSM4088788) has failed due to issues with raw data download, lets do it manually
SERIES=GSE137804
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list

# get urls
while read srr
do
 ../../actions/hdca_v2_reprocessing/utils/query_sdl.sh $srr
done < ${SERIES}.sample.list > srr_urls.txt

# download
while read srr url
do
  echo $srr
  wget -b -O done_wget/$srr  -o ${srr}_wget.log  $url
done <  srr_urls.txt

# stopped here
# seems like the rest was downloaded as fq, so should be save just tun it....:
nohup ./convert_to_fastq.sh $SERIES > to_fq.log 2>&1 &
./reorganise_fastqs.sh $SERIES


while read s
do
 rm -rf $s
done < ${SERIES}.sample.list

# restart mapping
nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 &

./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list

#############################
## GSE142526 ################
#############################

# some samples (GSM4231316 GSM4231320 GSM4231321 GSM4231323 GSM4231324)
# didn't finish succesfully because UMI has unusuall length - 5nt
SERIES=GSE142526
cd $SERIES

# detect unsuccesfull samples:
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample.list



# change expected length of UMI to 5 in the script
sed -i 's/$L3 == 10/$L3 == 5/' bam_to_10x_fastq_gz.sh
./convert_to_fastq.sh  $SERIES
mv fastqs fastqs_old
./reorganise_fastqs.sh $SERIES

# change UMI setting in STARsolo script
sed -i 's/UMILEN=10/UMILEN=5/g' starsolo_10x_auto.sh
sed -i 's/$R1LEN < 24/$R1LEN < 14/g' starsolo_10x_auto.sh

while read s
do
 rm -rf $s
done < ${SERIES}.sample.list

# restart mapping
nohup ./run_starsolo.sh $SERIES  > pm19_solo.log 2>&1 &

./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list
# two samples (GSM4231317 and GSM4231318) seems to be technical replicates as the have almost identical list of barcodes, but on GEO they have individual count matrices (the also have high barcode overlap)
# these samples doesn't seems to be used in HDCA
