
module load jellyfish/2.2.8


cat Alu.primate AluY AluYa5 AluYb8 L1.3 L1PA2 SVA_D SVA_E> SVA_F.outbound
jellyfish count -m 23 -C SVA_F.outbound -s 100M -o SVA_F.outbound.23mer

rm SVA_F.CC.23mer.count.1
rm SVA_F.CC.23mer.count.3

while read -r a
do

jellyfish query SVA_F.outbound.23mer ${a} >> SVA_F.CC.23mer.count.1
jellyfish query SVA_F.23mer ${a} >> SVA_F.CC.23mer.count.3

wait

done < SVA_F.CC.23mer.list
paste SVA_F.CC.23mer.count.1 SVA_F.CC.23mer.count.3 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > SVA_F.CC.23mer.list.candidate

rm SVA_F.GG.23mer.count.1
rm SVA_F.GG.23mer.count.3

while read -r a
do

jellyfish query SVA_F.outbound.23mer ${a} >> SVA_F.GG.23mer.count.1
jellyfish query SVA_F.23mer ${a} >> SVA_F.GG.23mer.count.3

wait

done < SVA_F.GG.23mer.list
paste SVA_F.GG.23mer.count.1 SVA_F.GG.23mer.count.3 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > SVA_F.GG.23mer.list.candidate




cat Alu.primate AluY AluYa5 AluYb8 L1.3 L1PA2 SVA_D SVA_F> SVA_E.outbound
jellyfish count -m 23 -C SVA_E.outbound -s 100M -o SVA_E.outbound.23mer

rm SVA_E.CC.23mer.count.1
rm SVA_E.CC.23mer.count.3

while read -r a
do

jellyfish query SVA_E.outbound.23mer ${a} >> SVA_E.CC.23mer.count.1
jellyfish query SVA_E.23mer ${a} >> SVA_E.CC.23mer.count.3

wait

done < SVA_E.CC.23mer.list
paste SVA_E.CC.23mer.count.1 SVA_E.CC.23mer.count.3 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > SVA_E.CC.23mer.list.candidate

rm SVA_E.GG.23mer.count.1
rm SVA_E.GG.23mer.count.3

while read -r a
do

jellyfish query SVA_E.outbound.23mer ${a} >> SVA_E.GG.23mer.count.1
jellyfish query SVA_E.23mer ${a} >> SVA_E.GG.23mer.count.3

wait

done < SVA_E.GG.23mer.list
paste SVA_E.GG.23mer.count.1 SVA_E.GG.23mer.count.3 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > SVA_E.GG.23mer.list.candidate



rm SVA_E.*.list.candidate.count

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> SVA_E.CC.list.candidate.count

wait

done < SVA_E.CC.23mer.list.candidate

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> SVA_E.GG.list.candidate.count

wait

done < SVA_E.GG.23mer.list.candidate


rm SVA_F.*.list.candidate.count

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> SVA_F.CC.list.candidate.count

wait

done < SVA_F.CC.23mer.list.candidate

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> SVA_F.GG.list.candidate.count

wait

done < SVA_F.GG.23mer.list.candidate

