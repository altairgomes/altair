entrada="Leda.eph"
saida1="Leda_cor.eph"
saida2="Leda_err.eph"
limitedist=800.0
limiteelev=15.0
while read -a linha
do
n1=${linha[20]}
n2=${linha[21]}
elev=${linha[23]}
dist=$(echo "sqrt($n1*$n1+$n2*$n2)" | bc)
if [ $(echo "$dist > $limitedist" | bc) -ne 0 ] && [ $(echo "$elev > $limiteelev" | bc) -ne 0 ];then
   grep "${linha[2]}" $entrada >> $saida1
 else
   grep "${linha[2]}" $entrada >> $saida2
fi
done < $entrada
