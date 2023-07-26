outdir=/ldfssz1/ST_BI/USER/zhaozijian/parsefq/test
barcode=/jdfssz1/ST_SUPERCELLS/PUB/snATAC/pipeline/03.scATAC_v2/config/C4scATAClib_seqT1_R1_70_R2_50.json
fastq1=/zfssz8/CNGB_DATA/BGISEQ01/DIPSEQ/DIPSEQT1/P21Z10200N0090_Temp/ALL-467-N-230223-ATAC-2/230320_SEQ102_DP8450009866BL_L01_SP2303130451/DP8450009866BL_L01_6_1.fq.gz
fastq2=/zfssz8/CNGB_DATA/BGISEQ01/DIPSEQ/DIPSEQT1/P21Z10200N0090_Temp/ALL-467-N-230223-ATAC-2/230320_SEQ102_DP8450009866BL_L01_SP2303130451/DP8450009866BL_L01_6_2.fq.gz
root=/jdfssz1/ST_SUPERCELLS/PUB/scRNA/pipeline/v3.1.5/common

mkdir -p ${outdir}/01.DataStat

#export LD_LIBRARY_PATH=/ldfssz1/ST_BI/USER/software/lib/:$LD_LIBRARY_PATH
#export PATH=/ldfssz1/ST_BI/USER/zhaofuxiang/lib/gcc-9.1.0/bin:$PATH
#export LD_LIBRARY_PATH="/ldfssz1/ST_BI/USER/zhaofuxiang/lib/gcc-9.1.0/lib:/ldfssz1/ST_BI/USER/zhaofuxiang/lib/gcc-9.1.0/lib64:$LD_LIBRARY_PATH"

echo "$fastq1" > ${outdir}/01.DataStat/in1.list
echo "$fastq2" > ${outdir}/01.DataStat/in2.list

bcPara=${outdir}/01.DataStat/bc.Para
echo "in1=${outdir}/01.DataStat/in1.list" > $bcPara
echo "in2=${outdir}/01.DataStat/in2.list" >> $bcPara
echo "config=${barcode}" >> $bcPara
echo "cbdis=${outdir}/01.DataStat/barcode_counts_raw.txt" >> $bcPara
echo "report=${outdir}/01.DataStat/sequencing_report.txt" >> $bcPara
echo "outFq1=${outdir}/01.DataStat/read_1.fq" >> $bcPara
echo "outFq2=${outdir}/01.DataStat/read_2.fq" >> $bcPara
echo "threads=10" >> $bcPara
echo "meanQual=4" >> $bcPara

${root}/bin/parseFqDev $bcPara
