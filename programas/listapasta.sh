atual=$(pwd)
for i in *;do
 if [ -d $i ];then
  l=0
  cd $i
  for k in *;do
   if [ -d $k ]; then
    l=1
    echo "$i/$k" >> $atual/lista_pastas
   fi
  done
  if [ "$l" -eq 0 ];then
   echo "$i" >> $atual/lista_pastas
  fi
  cd ..
 fi
done
