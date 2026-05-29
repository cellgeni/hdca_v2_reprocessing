# download reprocessing code
wget https://github.com/cellgeni/reprocess_public_10x/archive/0aa3606d57003407f1080777d646977a0ae5e15e.tar.gz
tar -xzvf 0aa3606d57003407f1080777d646977a0ae5e15e.tar.gz
mv reprocess_public_10x-0aa3606d57003407f1080777d646977a0ae5e15e reprocess_public_10x

# run reprocessing
# fill datasets.txt with ids to be reprocessed
while read d
do
 echo $d
 nohup ../actions/reprocess_public_10x/reprocess_public_10x.sh $d > ${d}.log 2>&1 &
done < datasets2.txt


###########################################################################################
# some datasets required some manual adjustments, these are explained below ###############
# some submissions has ATAC libraries, these libraries were skipped #######################
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
## GSE137804 ################
#############################
# some samples (GSM4088774 GSM4088775 GSM4088776 GSM4088777 GSM4088778 GSM4088779 GSM4088780 GSM4088781 GSM4088782 GSM4088783 GSM4088785 GSM4088786 GSM4088787 GSM4088788) 
# has failed due to issues with raw data download, lets do it manually
SERIES=GSE137804
cd $SERIES
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list

# get urls
while read gsm
do
 srr=`grep $gsm ${SERIES}.accessions.tsv | cut -f 4`
 ../../actions/hdca_v2_reprocessing/utils/query_sdl.sh $srr
done < ${SERIES}.sample.list | cut -d' ' -f1 > srr_urls.txt

# download
while read srr url
do
  echo $srr
  wget -b -O done_wget/$srr  -o ${srr}_wget.log  $url
done <  srr_urls.txt

  
# ./convert_to_fastq.sh  $SERIES # this one tries to rename sra files, but there is not need, so just use:
./bsub_sra2fastq.sh $SERIES sra_to_10x_fastq_gz.sh 

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
nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 &

./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list
# two samples (GSM4231317 and GSM4231318) seems to be technical replicates as the have almost identical list of barcodes, but on GEO they have individual count matrices (the also have high barcode overlap)
# these samples doesn't seem to be used in HDCA


#############################
## GSE155683 ################
#############################
# downloading failed for all samples
SERIES=GSE155683
cd $SERIES

# get urls
while read gsm
do
 srr=`grep $gsm ${SERIES}.accessions.tsv | cut -f 4`
 ../../actions/hdca_v2_reprocessing/utils/query_sdl.sh $srr
done < ${SERIES}.sample.list | cut -d' ' -f1 > srr_urls.txt

# download
while read srr url
do
  echo $srr
  wget -b -O done_wget/$srr  -o ${srr}_wget.log  $url
done <  srr_urls.txt

# ./convert_to_fastq.sh  $SERIES # this one tries to rename sra files, but there is no need, so just use:
./bsub_sra2fastq.sh $SERIES sra_to_10x_fastq_gz.sh  
rm -rf  GSM47104*
./reorganise_fastqs.sh $SERIES
nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 & 
./solo_QC.sh > $SERIES.solo_qc.tsv
# these are four atac samples, lets delete them (GSM4710478 GSM4710479 GSM4710480 GSM4710481)
for s in `grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample.list`
do
 echo $s
 rm -r $s
done

#############################
## PRJEB77091 ###############
#############################
# most samples are visium/TCR, visiums has failed, but TCR worked, lets delete them
rm -rf ERS20601045 ERS20601046
./solo_QC.sh > $SERIES.solo_qc.tsv
# delete failed (visiums) samples
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample.list | while read e; do rm -rf $e; done

#############################
## E-MTAB-8813 ##############
#############################
SERIES=E-MTAB-8813
cd $SERIES

# detect unsuccesfull samples:
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample.list
# raw data are unavailable for some of sampes


#############################
## E-MTAB-10552 #############
#############################
SERIES=E-MTAB-10552
cd $SERIES
# looks links from sdrf are broken
#ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10552/FCAImmP7862084_S1_L001_R1_001.fastq.gz # doesn't work
# but these seems to be working:
#https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/552/E-MTAB-10552/Files/FCAImmP7862084_S1_L001_R1_001.fastq.gz # works


mv ${SERIES}.urls.list ${SERIES}.urls_orig.list

sed  's|ftp://ftp\.ebi\.ac\.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10552/|https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/552/E-MTAB-10552/Files/|g' ${SERIES}.urls_orig.list > ${SERIES}.urls.list

nohup ./continuous_download.sh $SERIES > download.log 2>&1 & 

./reorganise_fastqs.sh $SERIES
nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 & 
./solo_QC.sh > $SERIES.solo_qc.tsv


#############################
## E-MTAB-8901 ##############
#############################
SERIES=E-MTAB-8901
cd $SERIES

mv ${SERIES}.sample.list ${SERIES}.sample_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list

# few samples (ERS4414927 ERS4414930 ERS4414932 ERS4414939 ERS4414946 ERS4414953)
# have strange layout - R2 is very short (26nt). lets try to accomodate them
sed -i 's/elif (( $R2LEN < 40 ))/elif (( $R2LEN < 20 ))/' starsolo_10x_auto.sh

for s in `cat ${SERIES}.sample.list`
do
 echo $s
 rm -r $s
done

nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 & 
./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list

#############################
## GSE171892 ################
#############################
SERIES=GSE171892
cd $SERIES
# GSM5236519 failed, likely due to sra->fastq conversion as fastq are corrupted
# delete and rerun
rm -rf GSM5236519 fastqs/GSM5236519
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
mv $SERIES.urls.list $SERIES.urls_all.list
grep -Fxv -f <(cut -f1 ${SERIES}.solo_qc.tsv | tail -n +2) ${SERIES}.sample_all.list > ${SERIES}.sample.list

grep -f <(grep GSM5236519 GSE171892.accessions.tsv | cut -f4 | sed 's/,/\n/g') $SERIES.urls_all.list > $SERIES.urls.list
nohup ./continuous_download.sh $SERIES > download.log 2>&1 & 
./reorganise_fastqs.sh $SERIES
nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 & 

./solo_QC.sh > $SERIES.solo_qc.tsv

mv ${SERIES}.sample_all.list ${SERIES}.sample.list


#############################
## GSE245310 ################
#############################
# in work
SERIES=GSE245310
cd $SERIES
mv ${SERIES}.sample.list ${SERIES}.sample_all.list
#mv $SERIES.urls.list $SERIES.urls_all.list
# so, they sumbitted lanes as individual samples... 
# here are real samples.
wget -q -O - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245310/suppl/GSE245310%5Fhuman%5FDRG%5Fmeta%2Etxt%2Egz | zcat | cut -f4 | tail -n +2 | sort | uniq -c


find . -type d -path "*/output" -empty | cut -d'/' -f 2 > ${SERIES}.sample.list
# 52 samples failed with segmentation fault... lets rerun them
for i in `cat ${SERIES}.sample.list`; do rm -rf $i; done


# there is a rumor that it is because of Velocyto, lets skip this stage
sed -i 's|Velocyto/raw Velocyto/filtered||' starsolo_10x_auto.sh
sed -i 's|Velocyto ||g' starsolo_10x_auto.sh

nohup ./run_starsolo.sh $SERIES  > solo.log 2>&1 & 
# looks like removal of velocyto helped
./solo_QC.sh > $SERIES.solo_qc.tsv
mv ${SERIES}.sample_all.list ${SERIES}.sample.list


#############################
# find incomplete starsolos #
#############################
find . -type d -path "*/*/output" -empty
