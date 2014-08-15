atual=$(pwd)
pastas=$(cat lista_pastas)
for noite in $pastas;do
  egrep "$noite" astrometry_findscale_01 > /dev/null;
  if [ $? -eq 1 ]; then
    echo $noite >> lista_pastas3;
  fi
done
