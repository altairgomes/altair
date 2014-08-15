arquivo="$1"
i=$(cat $arquivo | wc -l)
j=$(echo "$i-2" | bc)
k=$(echo "$j-43" | bc)
head -$j $arquivo > .aux1
tail -$k .aux1 > .aux2
while read -a linha;do
  echo "${linha[25]}  JD" >> jd_$arquivo
done < .aux2
rm .aux1 .aux2
