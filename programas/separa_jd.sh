while read -a linha
do
n1=${linha[13]}
echo "$n1 JD" >> julian_date
done < output
