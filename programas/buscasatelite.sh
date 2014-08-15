atual=$(pwd)
pasta="$atual/efemerides"
listasat=$(cat $pasta/satelites)
listapasta=$(cat lista_pastas)
for noite in $listapasta;do
  cd $noite
  echo "noite $noite"
  for satelite in $listasat;do
    cd $satelite
    while read -a linha
      do
      n1=${linha[13]}
      grep "$n1" $pasta/$satelite/targets_$satelite >> targets_$satelite;
    done < output
    ls *2mass.red.xy > lista_xy_2mass
    ls *ucac2.red.xy > lista_xy_ucac2
    ls *ucac4.red.xy > lista_xy_ucac4
    cp $atual/PRAIA_targets_search_20_02 .
    cp $pasta/$satelite/PRAIA_targets_search_20_02.dat .
    ./PRAIA_targets_search_20_02 < PRAIA_targets_search_20_02.dat
    cd $atual/$noite
  done
  echo "noite $noite terminada"
  cd $atual
done
