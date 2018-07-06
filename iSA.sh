[ $# != 3 ] && { echo "Usage: iSA.sh <name> <t1> <t2>"
                 echo "<name> name of the bam file from Bismark (only the string before '.bam')"
                 echo "<t1> number of bases uniformly trimmed off from 5' end of mate 1 reads (0 if no trimming)"
                 echo "<t2> number of bases uniformly trimmed off from 5' end of mate 2 reads (0 if no trimming)"
                 exit 1
               }

echo Checking library type...
pbat=`samtools view -H $1.bam | grep ' --pbat' | wc -l | cut -f 1 -d ' '`
if [ $pbat -gt 0 ]
then echo Library appears to be made with PBAT method and not compatible with iSA.
     echo Exiting...
     exit
fi

echo Writing bam file for Watson strand...
samtools view -h $1.bam | awk '$3!="chrM" && $3!="ChrM" && $16!="XG:Z:GA"' | samtools view -bh -@ 8 - > $1.Wat.bam;
echo Writing bam file for Crick strand...
samtools view -h $1.bam | awk '$3!="chrM" && $3!="ChrM" && $16!="XG:Z:CT"' | samtools view -bh -@ 8 - > $1.Cri.bam;
echo Writing bed file for Watson strand...
bamToBed -bedpe -i $1.Wat.bam | awk -v t1=$2 -v t2=$3 '{if($2<t1){$9=0} else{$9=$2-t1}; print $1,$9,$6+t2,NR,$7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -u -S 64G > $1.Wat.bed;
echo Writing bed file for Crick strand...
bamToBed -bedpe -i $1.Cri.bam | awk -v t1=$2 -v t2=$3 '{if($2<t2){$9=0} else{$9=$2-t2}; print $1,$9,$6+t1,NR,$7}' OFS='\t' | sort -k1,1 -k2,2n -k3,3n -u -S 64G > $1.Cri.bed;
echo Heating...
echo Cooling down... Strands annealing...$'\n'
intersectBed -a $1.Wat.bed -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | cut -f 1-5,9,10 > $1.iSA.bed;
num=`wc -l $1.iSA.bed | cut -f 1 -d ' '`
echo You found $num pairs of alignments between Watson and Crick. | tee $1.iSA.report.txt
wa=`wc -l $1.Wat.bed | cut -f 1 -d ' '`
cr=`wc -l $1.Cri.bed | cut -f 1 -d ' '`
wat=$(echo "scale=2; $num*100/$wa" | bc)
cri=$(echo "scale=2; $num*100/$cr" | bc)
echo Pairing efficiency: $wat% for Watson, $cri% for Crick.$'\n' | tee -a $1.iSA.report.txt
echo Calculating fold enrichment over random pairing...
N30=`awk '{if($2>29) print $1,$2-30,$3-30}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $N30 pairs of alignments found with -30 bp shift, | tee -a $1.iSA.report.txt
N20=`awk '{if($2>19) print $1,$2-20,$3-20}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $N20 pairs of alignments found with -20 bp shift, | tee -a $1.iSA.report.txt
N10=`awk '{if($2>9) print $1,$2-10,$3-10}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $N10 pairs of alignments found with -10 bp shift, | tee -a $1.iSA.report.txt
P10=`awk '{print $1,$2+10,$3+10}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $P10 pairs of alignments found with +10 bp shift, | tee -a $1.iSA.report.txt
P20=`awk '{print $1,$2+20,$3+20}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $P20 pairs of alignments found with +20 bp shift, | tee -a $1.iSA.report.txt
P30=`awk '{print $1,$2+30,$3+30}' OFS='\t' $1.Wat.bed | intersectBed -a - -b $1.Cri.bed -wa -wb -f 1 -F 1 -sorted | wc -l`
echo $P30 pairs of alignments found with +30 bp shift, | tee -a $1.iSA.report.txt
fold=$(echo "scale=1; $num*6/($N30+$N20+$N10+$P10+$P20+$P30)" | bc)
echo You got $fold-fold enrichment over random pairing.$'\n' | tee -a $1.iSA.report.txt

read -p "Do you want to proceed (Y/N)?" answer
case $answer in
Y|y) echo Continue on... Extracting paired alignments...
     cut -f 4 $1.iSA.bed | awk '{print $1*2-1"\n"$1*2}' > $1.Wat.LN;
     cut -f 6 $1.iSA.bed | awk '{print $1*2-1"\n"$1*2}' > $1.Cri.LN;
     samtools view -H $1.Wat.bam > header;
     samtools view $1.Wat.bam | awk 'FNR==NR {h[$1];next} (FNR in h)' $1.Wat.LN - | cat header - | samtools view -bh -@ 48 - > $1.Wat.iSA.bam;
     samtools view $1.Cri.bam | awk 'FNR==NR {h[$1];next} (FNR in h)' $1.Cri.LN - | cat header - | samtools view -bh -@ 48 - > $1.Cri.iSA.bam;
     echo Extracting DNA methylation calls...
     bismark_methylation_extractor -p --multicore 8 --no_header --gzip $1.Wat.iSA.bam
     bismark_methylation_extractor -p --multicore 8 --no_header --gzip $1.Cri.iSA.bam
     echo Pairing DNA methylation calls between Watson and Crick...
     awk '{print $5"\t"NR}' $1.iSA.bed | sort -k1,1 > $1.Wat.ID;
     awk '{print $7"\t"NR}' $1.iSA.bed | sort -k1,1 > $1.Cri.ID;
     gzip -cd CpG_OT_$1.Wat.iSA.txt.gz | sort -k1,1 | join -j 1 - $1.Wat.ID | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4,$3,$4-1,$4+1,$5}' OFS='\t' | sort -k1,1 > $1.Wat.CG.me;
     gzip -cd CpG_OB_$1.Cri.iSA.txt.gz | sort -k1,1 | join -j 1 - $1.Cri.ID | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4-1,$3,$4-2,$4,$5}' OFS='\t' | sort -k1,1 > $1.Cri.CG.me;
     gzip -cd CHG_OT_$1.Wat.iSA.txt.gz | sort -k1,1 | join -j 1 - $1.Wat.ID | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4,$3,$4-1,$4+2,$5}' OFS='\t' | sort -k1,1 > $1.Wat.CHG.me;
     gzip -cd CHG_OB_$1.Cri.iSA.txt.gz | sort -k1,1 | join -j 1 - $1.Cri.ID | sed 's/ /\t/g' | awk '{print $6"_"$3"_"$4-2,$3,$4-3,$4,$5}' OFS='\t' | sort -k1,1 > $1.Cri.CHG.me;
     join -j 1 $1.Wat.CG.me $1.Cri.CG.me | sed 's/ /\t/g' | awk '{if($5=="z" && $9=="z"){print $2,$3,$4,1,0,0,0} else if($5=="z" && $9=="Z"){print $2,$3,$4,0,1,0,0} else if($5=="Z" && $9=="z"){print $2,$3,$4,0,0,1,0} else{print $2,$3,$4,0,0,0,1}}' OFS='\t' | sort -k1,1 -k2,2n | groupBy -g 1-3 -c 4,5,6,7 -o sum > $1.intraCpG.bed;
     join -j 1 $1.Wat.CHG.me $1.Cri.CHG.me | sed 's/ /\t/g' | awk '{if($5=="x" && $9=="x"){print $2,$3,$4,1,0,0,0} else if($5=="x" && $9=="X"){print $2,$3,$4,0,1,0,0} else if($5=="X" && $9=="x"){print $2,$3,$4,0,0,1,0} else{print $2,$3,$4,0,0,0,1}}' OFS='\t' | sort -k1,1 -k2,2n | groupBy -g 1-3 -c 4,5,6,7 -o sum > $1.intraCWG.bed;
     awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$5; HW+=$6; Me+=$7} END{print "You found "To" intraCpGs, of which:\n"Un" are unmethylated,\n"HW" are hemi-Watson,\n"HC" are hemi-Crick,\n"Me" are methylated.\n"}' $1.intraCpG.bed | tee -a $1.iSA.report.txt
     awk '{To+=$4+$5+$6+$7; Un+=$4; HC+=$5; HW+=$6; Me+=$7} END{print "You found "To" intraCWGs, of which:\n"Un" are unmethylated,\n"HW" are hemi-Watson,\n"HC" are hemi-Crick,\n"Me" are methylated.\n"}' $1.intraCWG.bed | tee -a $1.iSA.report.txt;;
N|n) echo Cleaning up...
     rm $1.Wat.bam $1.Cri.bam $1.Wat.bed $1.Cri.bed
     exit;;
esac

read -p "Remove temp files (Y/N)?" clean
case $clean in
Y|y) echo Cleaning up...
     rm $1.Wat.bam $1.Cri.bam $1.Wat.bed $1.Cri.bed $1.Wat.LN $1.Cri.LN header $1.Wat.iSA.M-bias_R?.png $1.Cri.iSA.M-bias_R?.png $1.Wat.ID $1.Cri.ID $1.Wat.CG.me $1.Cri.CG.me $1.Wat.CHG.me $1.Cri.CHG.me;;
N|n) exit;;
esac

echo The end of iSA.
