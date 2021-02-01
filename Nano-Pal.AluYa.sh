echo "Nano-Pal Version1.0 Jan 22nd 2021"

#PB_DATA='/home/arthurz/arthur_remillsscr/19.01.22.NA12878_PacBio_pbmm2/19.01.29.alignment'
lib='lib'

#####
echo "Enter the your output prefix: "  
read FC 
#FC='FAO08905'

#####
echo "Enter the your element (ALU): "  
read MEI 
MEI='ALU'


#####
echo "Enter the your input data: "  
read RAW
#RAW='/home/arthurz/arthur_remillsscr/19.10.26.Nanopore.validaton/workspace.12.29.20.rebatch/data/'${FC}


######
echo "Enter the your input reference GRCh38.fa file: "  
read REF
#REF='/nfs/turbo/dcmb-brainsom/technical/reference/GRCh38_bsm_reference_genome/GRCh38_BSM.fa'

######
echo "Enter the your reference version (GRCh38 only for now): "  
read REF_VER
#REF_VER='GRCh38'

######
echo "Enter the path of minimap2: "  
read MINIMAP
#MINIMAP='/home/arthurz/app/minimap2/minimap2'

######
echo "Enter the path of PALMER: "  
read PALMER
#PALMER='/home/arthurz/arthur_remillsscr/18.02.23.PALMER.upgrade/v1.7.2/PALMER'

#MEI reference sequence
L1Hs=${lib}'/L1.3'
AluYa5=${lib}'/AluYa5'
AluYb8=${lib}'/AluYb8'
SVA_E=${lib}'/SVA_E'
SVA_F=${lib}'/SVA_F'
#ref MEI information from RepeatMasker
S_refLINE=${lib}'/hg38.RM.L1.ref'
S_refALU=${lib}'/hg38.RM.ALU.ref'
S_refSVA=${lib}'/hg38.RM.SVA.ref'

#PALMER call for NA12878
S_origLINE=${lib}'/PALMER.NA12878.L1.txt'
S_origALU=${lib}'/PALMER.NA12878.ALU.txt'
S_origSVA=${lib}'/PALMER.NA12878.SVA.txt'


#P&P call for NA12878
S_PP_LINE=${lib}'/union.1214/L1.inter.fi'
S_PP_ALU=${lib}'/union.1214/ALU.inter.fi'
S_PP_SVA=${lib}'/union.1214/SVA.inter.fi'


mv ${FC}.workspace.${MEI}.ya ${FC}.workspace.old
mkdir ${FC}.workspace.${MEI}.ya
cd ${FC}.workspace.${MEI}.ya

##RAW fastq data & generate fasta
cat ${RAW}/*.fastq > batch.fastq
cat batch.fastq  | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > batch.fasta 

wait 

##alignment
${MINIMAP} -ax map-ont  ${REF} batch.fastq  >  Nanopore.sam 
samtools sort Nanopore.sam  > Nanopore.sorted.bam 
samtools index Nanopore.sorted.bam

wait 
mkdir palmer

##PALMER
##pull out L1Hs signals
while read -r chr
do
${PALMER} --input Nanopore.sorted.bam --workdir palmer/ --output ${chr} --ref_ver ${REF_VER} --ref_fa ${REF} --type ${MEI} --chr ${chr}
done < ${lib}/chr.list

wait

##pull out all reads having putative L1Hs signal reported by PALMER
fasta_i=batch.fasta

rm blastn_refine.all.txt
cat palmer/chr*_*_*/blastn_refine.txt >>  blastn_refine.all.txt 

###Tool: Cigar information parsing
cp ${lib}/cigar.parser.1106.o .

samtools view Nanopore.sorted.bam -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 | awk '{print $1,$3,$4,$6}' > mapped.info.txt

rm cigar_results.all.txt
while read -r a b c d
do
echo ${d} > input_cigar
./cigar.parser.1106.o
wait
cat cigar_results >> cigar_results.all.txt
done < mapped.info.txt

##Get the information of start and end for each read with L1Hs signal
awk '{print $4+$6+$10}' cigar_results.all.txt > cigar.ref.txt
paste mapped.info.txt cigar.ref.txt | awk '{print $1,$2,$3,$3+$5}' >  mapped.info.final.txt

rm	5.plus.txt 
rm	5.minus.txt 
rm 	3.plus.txt 
rm 	3.minus.txt 
rm read.length.txt
rm read.name.txt
mkdir blastn

#blastn -evalue 0.001 -task blastn -gapopen 4 -query ${L1Hs} -subject ${fasta_i} -outfmt "7 qacc sacc evalue qstart qend sstart send" > blastn/all.blastn.txt

line=$(wc -l $fasta_i| awk '{print $1}')

##Grep the L1Hs sequence details (in L1Hs reference sequence) using blastn  
for((i=0;i<${line}/2;i++))
do
let t=${line}-i*2
tail -n ${t} ${fasta_i} | head -n 2 > blastn/${t}.fasta

#######MEI
blastn -evalue 0.001 -task blastn -gapopen 4 -query ${AluYa5} -subject blastn/${t}.fasta -outfmt "7 qacc sacc evalue qstart qend sstart send" > blastn/${t}.blastn.txt

awk 'NR%2'	blastn/${t}.fasta | awk '{print $1}' | sed -e "s/>//g" >> read.name.txt
len=$(awk '!(NR%2)'	blastn/${t}.fasta | awk '{print length($1)}')
let len_fix=len-100
echo ${len} >> read.length.txt

#######MEI
cat blastn/${t}.blastn.txt | grep -v "#"  | awk 'BEGIN{flag=0}{if($5>225&&($6<100||$7<100)&&($6<$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' >> 5.plus.txt 
cat blastn/${t}.blastn.txt | grep -v "#"  | awk 'BEGIN{flag=0}{if($5>225&&($6<100||$7<100)&&($6>$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' >> 5.minus.txt 
cat blastn/${t}.blastn.txt | grep -v "#"  | awk 'BEGIN{flag=0}{if($5>225&&($6>'$len_fix'||$7>'$len_fix')&&($6<$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' >> 3.plus.txt 
cat blastn/${t}.blastn.txt | grep -v "#"  | awk 'BEGIN{flag=0}{if($5>225&&($6>'$len_fix'||$7>'$len_fix')&&($6>$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' >> 3.minus.txt 

wait
done


paste read.name.txt read.length.txt 5.plus.txt 5.minus.txt 3.plus.txt 3.minus.txt > read.all.txt

wait

rm	palmer.5.plus.txt
rm	palmer.5.minus.txt
rm	palmer.3.minus.txt
rm	palmer.3.plus.txt
rm	palmer.map.txt

##Grep the polymorphic L1Hs sequence details (in L1Hs reference sequence) from PALMER
while read -r a b c d e f 
do
let b_fix=b-100

#######MEI
awk 'BEGIN{flag=0}{if($2~"'$a'"&&$5>225&&($6<100||$7<100)&&($6<$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' blastn_refine.all.txt >> palmer.5.plus.txt
awk 'BEGIN{flag=0}{if($2~"'$a'"&&$5>225&&($6<100||$7<100)&&($6>$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' blastn_refine.all.txt >> palmer.5.minus.txt
awk 'BEGIN{flag=0}{if($2~"'$a'"&&$5>225&&($6>'$b_fix'||$7>'$b_fix')&&($6<$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' blastn_refine.all.txt >> palmer.3.plus.txt
awk 'BEGIN{flag=0}{if($2~"'$a'"&&$5>225&&($6>'$b_fix'||$7>'$b_fix')&&($6>$7)) flag+=1} END{if(flag>=1) print "1";else print "0"}' blastn_refine.all.txt >> palmer.3.minus.txt

awk 'BEGIN{flag=0}{if($1=="'$a'") {flag+=1;print $2,$3,$4}} END{if(flag==0) print "NON","0","0"}' mapped.info.final.txt >> palmer.map.txt
wait
done < read.all.txt

paste read.all.txt palmer.5.plus.txt palmer.5.minus.txt palmer.3.plus.txt palmer.3.minus.txt palmer.map.txt > read.all.palmer.final.txt
wait

#cp read.all.palmer.txt read.all.palmer.final.txt

mkdir inter

##Tool: intersect with reference L1NE data from RepeatMasker and PALMER calls
cp ${lib}/inter inter/

##Tool: exclude redundant intersection for each read
cp ${lib}/RM_collapse.0921.o  inter/


cd inter

awk '{print $11,$12,$12+1,$1}' ../read.all.palmer.final.txt  | grep -v "NON" > read.5.Q.txt
awk '{print $11,$13,$13+1,$1}' ../read.all.palmer.final.txt  | grep -v "NON" > read.3.Q.txt

cp read.5.Q.txt Q.txt 

#######MEI
cp ${S_refALU}  S.txt
./inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt 
./RM_collapse.0921.o 
mv output_rm_cluster.txt read.5.Q.ref.txt

#######MEI
cp ${S_origALU}  S.txt
./inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt| awk '{print $1,$2,$3,$5,$4}'| sort -n | uniq -f 4 |  awk '{print $1,$2,$3,$5,$4}'  > read.5.Q.P.txt

cp read.3.Q.txt Q.txt 

#######MEI
cp ${S_refALU}  S.txt
./inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt > input_rm_cluster.txt 
./RM_collapse.0921.o 
mv output_rm_cluster.txt read.3.Q.ref.txt

#######MEI
cp ${S_origALU}  S.txt
./inter
awk '{if($8=="") print $1,$2,$3,$4,"NON"; else print $1,$2,$3,$4,$8}' inter.txt| awk '{print $1,$2,$3,$5,$4}'| sort -n | uniq -f 4 |  awk '{print $1,$2,$3,$5,$4}'  > read.3.Q.P.txt

cd ..

wait

sort -k 4 inter/read.5.Q.P.txt > 5.P.inter.txt
sort -k 4 inter/read.3.Q.P.txt > 3.P.inter.txt
sort -k 4 inter/read.5.Q.ref.txt > 5.ref.inter.txt
sort -k 4 inter/read.3.Q.ref.txt > 3.ref.inter.txt

paste 5.P.inter.txt 3.P.inter.txt 5.ref.inter.txt 3.ref.inter.txt | awk '{print $4,$5,$10,$15,$20}' > read.RM.txt

rm read.RM.2.txt
while read -r a b c d e f g h i j k l m
do
awk 'BEGIN{score=0}{if($1=="'$a'") {score+=1;print $2,$3,$4,$5}} END{if(score==0) print "NA","NA","NA","NA"}' read.RM.txt >> read.RM.2.txt
wait
done < read.all.palmer.final.txt

##Summary table of L1Hs for each read
paste read.all.palmer.final.txt read.RM.2.txt > summary.final.txt



##Tool: interscetion
mkdir inter.2
cp ${lib}/inter.0118.o inter.2/


##different types of LINE
#######MEI
less ${S_refALU} | grep AluYa > ref.AluYa5
less ${S_refALU} | grep AluYb > ref.AluYb8
less ${S_refALU} | grep -v AluYb | grep -v AluYa | grep AluY > ref.AluY
less ${S_refALU} | grep -v AluYb | grep -v AluYa | grep -v AluY | grep Alu > ref.Alu

##Valid mapping reads
samtools view Nanopore.sorted.bam -q 10 -F 0x100 -F 0x200 -F 0x800 -F 0x400 -f 0x10 | awk '{print $1}' >  RC.all.list

rm RC.tag
while read -r aa b c d e f g h i j  k l m n o p q r
do
awk 'BEGIN{flag=0}{if($1=="'$aa'") flag+=1} END{if(flag>=1) print "1";else print "0"}' RC.all.list >> RC.tag

wait
done < summary.final.txt

paste summary.final.txt RC.tag  > summary.final.2.txt

SIGNAL1=$(less summary.final.2.txt |awk '{if((($3+$4)!=0&&($5+$6)==0)||(($3+$4)==0&&($5+$6)!=0)) print $0}' | wc -l)
SIGNAL2=$(less summary.final.2.txt |awk '{if(($3+$4)!=0&&($5+$6)!=0) print $0}' | wc -l)
FAIL=$(less summary.final.2.txt |awk '{if(($3+$4+$5+$6)==0) print $0}' | wc -l)



less summary.final.2.txt |awk '{if($11=="NON") print $0}' >  summary.final.unmap.read.txt

less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6)==0&&($7+$8+$9+$10)!=0) print $0}'>  summary.final.odd.read.txt

less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6+$7+$8+$9+$10)==0) print $0}'>  summary.final.no.read.txt

less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6)!=0) print $0}'>  summary.final.all.read.txt

#less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0) print $0}'>  summary.final.PALMER.L1.read.txt
less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6)!=0&&($7+$8+$9+$10)!=0) print $0}'>  summary.final.PALMER.read.txt


less summary.final.2.txt |awk '{if($11!="NON") print $0}' |awk '{if(($3+$4+$5+$6)!=0&&($7+$8+$9+$10)==0) print $0}'>  summary.final.ref.read.txt

#less summary.final.PALMER.L1.read.txt |awk '{if(($7+$8)>0&&($9+$10)==0) print $11,$12,$14,$18,"5"; else if (($7+$8)==0&&($9+$10)>0) print $11,$13,$15,$18,"3"; else if(($7+$8)>0&&($9+$10)>0) print $11,$12,$14,$18,"5""\n"$11,$13,$15,$18,"3"}' | sort -k 2 -n | sort -k 1 > capture.loci.palmer
less summary.final.PALMER.read.txt |awk '{if(($7+$8)>0&&($9+$10)==0&&$7>0) print $11,$12,$14,$18,"5","+";else if(($7+$8)>0&&($9+$10)==0&&$8>0) print $11,$12,$14,$18,"5","-";else if (($7+$8)==0&&($9+$10)>0&&$9>0) print $11,$13,$15,$18,"3","+";else if (($7+$8)==0&&($9+$10)>0&&$10>0) print $11,$13,$15,$18,"3","-"; else if(($7+$8)>0&&($9+$10)>0&&$7>0&&$9>0) print $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","+";else if(($7+$8)>0&&($9+$10)>0&&$7>0&&$10>0) print $11,$12,$14,$18,"5","+""\n"$11,$13,$15,$18,"3","-";else if(($7+$8)>0&&($9+$10)>0&&$8>0&&$9>0) print $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","+";else if(($7+$8)>0&&($9+$10)>0&&$8>0&&$10>0) print $11,$12,$14,$18,"5","-""\n"$11,$13,$15,$18,"3","-"}' | sort -k 2 -n | sort -k 1 > capture.loci.palmer

less summary.final.PALMER.read.txt |awk '{if(($7+$8)>0&&($9+$10)==0&&$18==0&&($5+$6)>0) print $11,$13,$17,$18,"3"; else if(($7+$8)>0&&($9+$10)==0&&$18==1&&($3+$4)>0) print $11,$13,$17,$18,"3"; else if (($7+$8)==0&&($9+$10)>0&&$18==0&&($3+$4)>0) print $11,$12,$16,$18,"5"; else if (($7+$8)==0&&($9+$10)>0&&$18==1&&($5+$6)>0) print $11,$12,$16,$18,"5"}' | sort -k 2 -n | sort -k 1 > capture.loci.ref.add

less summary.final.ref.read.txt |awk '{if(($3+$4)>0&&($5+$6)==0&&$18==0) print $11,$12,$16,$18,"5"; else if(($3+$4)>0&&($5+$6)==0&&$18==1) print $11,$13,$17,$18,"3"; else if (($3+$4)==0&&($5+$6)>0&&$18==0) print $11,$13,$17,$18,"3"; else if (($3+$4)==0&&($5+$6)>0&&$18==1) print $11,$12,$16,$18,"5"; else if(($3+$4)>0&&($5+$6)>0) print $11,$12,$16,$18,"5""\n"$11,$13,$17,$18,"3"}' | sort -k 2 -n | sort -k 1 > capture.loci.ref


less capture.loci.palmer | grep  cluster |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.palmer.process
#less capture.loci.palmer | grep  -v cluster |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}'  > capture.loci.potential.process
less capture.loci.palmer | grep  -v cluster |  awk '{print $1,$2,$2+1,$3,$4,$5,$6,"Nanopore"}'  > capture.loci.potential.process


cat capture.loci.ref capture.loci.ref.add > capture.loci.ref.all
less capture.loci.ref.all | grep AluYa |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYa5
less capture.loci.ref.all | grep -v AluYa | grep AluYb  |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYb8
less capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep AluY |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYother
less capture.loci.ref.all | grep -v AluYa | grep -v AluYb | grep -v AluY |  awk '{print $1,$2,$2+1,$3,$4,$5,"Nanopore"}' > capture.loci.r.AluYnon

cd inter.2

#######MEI
cp ${S_PP_ALU} Q.txt
cp ../capture.loci.palmer.process S.txt
./inter.0118.o
wait
mv inter.txt inter.palmer.txt

cp ../capture.loci.potential.process S.txt
./inter.0118.o
wait
mv inter.txt inter.palmer.add.txt

#######MEI
cp ../ref.AluYa5 Q.txt
cp ../capture.loci.r.AluYa5 S.txt
./inter.0118.o
wait
mv inter.txt inter.r.AluYa5.txt

#######MEI
cp ../ref.AluYb8 Q.txt
cp ../capture.loci.r.AluYb8 S.txt
./inter.0118.o
wait
mv inter.txt inter.r.AluYb8.txt

#######MEI
cp ../ref.AluY Q.txt
cp ../capture.loci.r.AluYother S.txt
./inter.0118.o
wait
mv inter.txt inter.r.AluY.txt

#######MEI
cp ../ref.Alu Q.txt
cp ../capture.loci.r.AluYnon S.txt
./inter.0118.o
wait
mv inter.txt inter.r.Alu.txt

##report the numbers of captured events 
#cat inter.palmer.txt inter.palmer.add.txt | grep Nano | awk '{print $4,$30 ,$17, $(NF-1)}' | sort | uniq -c > p.txt.fi.2
cat inter.palmer.txt inter.palmer.add.txt | grep Nano | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | sort | uniq -c > p.txt.fi

#######MEI
less inter.r.AluYa5.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.AluYa5.txt.fi
less inter.r.AluYb8.txt | grep Nano | awk '{print $1,$2,$3}' | sort | uniq -c > r.AluYb8.txt.fi
less inter.r.AluY.txt | grep Nano | awk '{print $1,$2,$3, $4}' | sort | uniq -c > r.AluY.txt.fi
less inter.r.Alu.txt | grep Nano | awk '{print $1,$2,$3, $4}' | sort | uniq -c > r.Alu.txt.fi


cp ../capture.loci.potential.process Q.txt

#######MEI
cp ${S_PP_ALU} S.txt
./inter.0118.o
wait
less inter.txt | grep -v cluster > inter.potential.txt

cd ..

#######MEI
cp ${lib}/cluster.1219.o .

cp inter.2/inter.potential.txt input_cluster.txt 
./cluster.1219.o
cp clustered.txt  potential.clustered.txt.fi


pwd

echo There are ${SIGNAL1} reads caputuring putative ${MEI} signals on one end.
echo There are ${SIGNAL2} reads caputuring putative ${MEI} signals on both ends.
echo There are ${FAIL} reads having no putative L1 signals.


less inter.2/p.txt.fi | awk '{sum+=$1} END {print "Non-reference AluY reads sum = ", sum}'
less inter.2/r.AluYa5.txt.fi | awk '{sum+=$1} END {print "Reference AluYa reads sum = ", sum}'
less inter.2/r.AluYb8.txt.fi | awk '{sum+=$1} END {print "Reference AluYb reads sum = ", sum}'
less inter.2/r.AluY.txt.fi | awk '{sum+=$1} END {print "Other reference AluY reads sum = ", sum}'
less capture.loci.r.AluYnon | grep Alu | wc -l | awk '{print "Other reference SVA reads sum = "$1}'
less potential.clustered.txt.fi | awk '{sum+=$4} END {print "Potential specific non-reference AluY reads sum = ", sum}'
echo "Number of non_Alu reads"
less capture.loci.r.AluYnon | grep -v Alu | wc -l | awk '{print $1}'

echo "Number of non-reference AluY and the file (number of suppoting reads + coordinate + overlap information)"
wc -l inter.2/p.txt.fi
echo "Number of reference AluYa and the file (number of suppoting reads + coordinate)"
wc -l inter.2/r.AluYa5.txt.fi
echo "Number of reference AluYb and the file (number of suppoting reads + coordinate)"
wc -l inter.2/r.AluYb8.txt.fi
echo "Number of other reference AluY and the file (number of suppoting reads + coordinate)"
wc -l inter.2/r.AluY.txt.fi
echo "Number of other reference Alu and the file (number of suppoting reads + coordinate)"
wc -l inter.2/r.Alu.txt.fi
echo "Number of potential specific non-reference AluY and the file (coordinate + number of suppoting reads + number of left suppoting reads + number of right suppoting reads + strand)"
wc -l potential.clustered.txt.fi



