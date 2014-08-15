atual=$(pwd)
pasta="$atual/efemerides"
listasat=$(cat $pasta/satelites)
listapasta=$(cat lista_pastas)
totallinhas=$(wc -l < lista_pastas)
numtotal=$(echo "$totallinhas * 10"|bc)
count=0
(
for noite in $listapasta;do
  cd $noite
  echo "noite $noite"
  for satelite in $listasat;do
    if [ ! -d $satelite ];then
      mkdir $satelite
    fi
    if [ ! -e $satelite/output ];then
      while read -a linha; do
        n1=${linha[13]}
        grep "$n1" $pasta/$satelite/targets_$satelite > /dev/null;
        if [ $? -eq 0 ]; then
          grep "$n1" output >> $satelite/output
        fi
      done < output
    fi
    if [ -e output ];then 
      cd $satelite
      cp $pasta/$satelite/targets_$satelite .
      cp $pasta/$satelite/PRAIA_header_edit_20_01.dat .
      cp $atual/PRAIA_header_edit_20_01 .
      ./PRAIA_header_edit_20_01 < PRAIA_header_edit_20_01.dat
      if [ -e output2 ];then
        rm output
        mv output2 output
      fi
    fi
    echo "Satelite $noite/$satelite convertido"
    rm targets_$satelite
    count=$[$count + 1]
    porcentagemtot=$(echo "$count * 100 / $numtotal"|bc -l)
    porcentagem=$(printf "%.2f" $porcentagemtot)
    echo "# Noite $noite; Satelite $satelite; Total $porcentagem %"
    echo "$porcentagemtot"
    cd $atual/$noite
  done
  echo "noite $noite terminada"
  cd $atual
done
) |
zenity --progress \
  --title="Gerando arquivos e convertendo" \
  --text="Gerando..." \
  --percentage=0
