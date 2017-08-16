while read line;do
  ds9 ${line:114:35} & 
  pid1=$!
  yad --text="${line:114:35}" --button="Boa":1 --button="Ruim":0
  res=$?
  kill -9 $pid1
  if [ "$res" = "1" ];then
    echo " $line" >> output2
  fi
done < output
