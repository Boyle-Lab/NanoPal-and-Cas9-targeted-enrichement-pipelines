

module load jellyfish/2.2.8

cat AluY Alu.primate L1.3 SVA_E SVA_F AluYa5 > AluYb8.outbound
jellyfish count -m 23 -C AluYb8.outbound -s 100M -o AluYb8.outbound.23mer

rm AluYb8.CC.23mer.count.1
rm AluYb8.CC.23mer.count.2
while read -r a
do

jellyfish query AluYb8.outbound.23mer ${a} >> AluYb8.CC.23mer.count.1
jellyfish query AluYb8.23mer ${a} >> AluYb8.CC.23mer.count.2

wait

done < AluYb8.CC.23mer.list
paste AluYb8.CC.23mer.count.1 AluYb8.CC.23mer.count.2 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > AluYb8.CC.23mer.list.candidate


rm AluYb8.GG.23mer.count.1
rm AluYb8.GG.23mer.count.2

while read -r a
do

jellyfish query AluYb8.outbound.23mer ${a} >> AluYb8.GG.23mer.count.1
jellyfish query AluYb8.23mer ${a} >> AluYb8.GG.23mer.count.2

wait

done < AluYb8.GG.23mer.list
paste AluYb8.GG.23mer.count.1 AluYb8.GG.23mer.count.2 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > AluYb8.GG.23mer.list.candidate






cat AluY Alu.primate L1.3 SVA_E SVA_F AluYb8 > AluYa5.outbound
jellyfish count -m 23 -C AluYa5.outbound -s 100M -o AluYa5.outbound.23mer

rm AluYa5.CC.23mer.count.1
rm AluYa5.CC.23mer.count.2
while read -r a
do

jellyfish query AluYa5.outbound.23mer ${a} >> AluYa5.CC.23mer.count.1
jellyfish query AluYa5.23mer ${a} >> AluYa5.CC.23mer.count.2

wait

done < AluYa5.CC.23mer.list
paste AluYa5.CC.23mer.count.1 AluYa5.CC.23mer.count.2 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > AluYa5.CC.23mer.list.candidate


rm AluYa5.GG.23mer.count.1
rm AluYa5.GG.23mer.count.2

while read -r a
do

jellyfish query AluYa5.outbound.23mer ${a} >> AluYa5.GG.23mer.count.1
jellyfish query AluYa5.23mer ${a} >> AluYa5.GG.23mer.count.2

wait

done < AluYa5.GG.23mer.list
paste AluYa5.GG.23mer.count.1 AluYa5.GG.23mer.count.2 | awk '{if($2==0&&$4==1&&$1==$3) print $1}' > AluYa5.GG.23mer.list.candidate


rm AluYb8.*.list.candidate.count

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> AluYb8.CC.list.candidate.count

wait

done < AluYb8.CC.23mer.list.candidate

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> AluYb8.GG.list.candidate.count

wait

done < AluYb8.GG.23mer.list.candidate


rm AluYa5.*.list.candidate.count

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> AluYa5.CC.list.candidate.count

wait

done < AluYa5.CC.23mer.list.candidate

while read -r a
do

jellyfish query ref.jf.build/ref.23mer ${a} >> AluYa5.GG.list.candidate.count

wait

done < AluYa5.GG.23mer.list.candidate
