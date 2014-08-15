for arquivo in search*;do
  while read -a linha;do
    echo "${linha[43]} JD" >> julian_date_${arquivo:15}
  done < $arquivo
done
