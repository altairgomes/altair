atual=$(pwd)
pasta="$atual/efemerides/IAG"
listasat=$(cat $pasta/satelites)
listapasta=$(cat lista_pastas)
for noite in $listapasta;do
  cd $noite
  for satelite in $listasat;do
    cd $satelite
    if [ -s output ]; then
      linhasoutput=$(wc -l < output)
      linhasastrometry=$(wc -l < astrometry_photometry_praia_star1)
      if [ $(echo "$linhasastrometry < $linhasoutput" | bc) -ne 0 ]; then
        echo -e "$noite/$satelite\t$linhasastrometry\t$linhasoutput" >> $atual/verificar_noites
      fi
    fi
    cd $atual/$noite
  done
  cd $atual
done
