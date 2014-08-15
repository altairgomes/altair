atual=$(pwd)
pasta="$atual/efemerides"
listasat=$(cat $pasta/satelites)
listapasta=$(cat lista_pastas)
progresso="$atual/progresso"
datain=$(date +%D-%T)
echo -e "Teste de processamento" >> $progresso
totimg=0
for noite in $listapasta;do
  cd $noite
  imgnoite=0
  echo -e "\n\nNoite $noite" >> $progresso
  datanoitein=$(date +%T)
  for satelite in $listasat;do
    cp $satelite/output .
    totallinhas=$(wc -l < output)
    totimg=$[$totimg + $totallinhas]
    imgnoite=$[$imgnoite + $totallinhas]
    datasatin=$(date +%T)
    ./PRAIA_astrometry_20_08 < PRAIA_astrometry_20_08.dat
    mv *.xy $satelite
    mv *.reg $satelite
    mv astrometry_* $satelite
    datasatfim=$(date +%T)
    echo "Tempo de processamento: inicio:$datasatin fim:$datasatfim; $totallinhas imagens do satelite $satelite" >> $progresso
  done
  datanoitefim=$(date +%T)
  echo -e "\nTotal na noite $noite inicio:$datanoitein fim:$datanoitefim total de imagens: $imgnoite" >> $progresso
  cd $atual
done
datafim=$(date +%D-%T)
echo -e "\n\nTotal processamento inicio:$datain fim:$datafim total de imagens: $totimg" >> $progresso
