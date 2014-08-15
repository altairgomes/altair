 n=0
 listacourse=(101 102 201 202)
 lim=7
 lim2=4
 unit=1
 level=1
 for course in ${listacourse[*]}; do
  if [ "$course" -gt "200" ]; then
  lim2=$[$lim2 - 1]
  level=2
  fi
  while [ $unit != $lim2 ]; do
   if [ $unit != 1 ]; then
   lim=6
   fi
   lesson=1
   while [ $lesson != $lim ]; do
    suite=1
    while [ $suite != 11 ]; do
     item=1
     while [ $item != 5 ]; do
      count=$[(suite*4)-4+$item]
      n=$[$n + 1]
      item=$[$item + 1]
     done
     suite=$[$suite + 1]
    done
    lesson=$[$lesson + 1]
    cd ..
   done
   unit=$[$unit + 1]
   cd ..
  done
  lim2=$[$lim2 + 3]
  cd ..
 done
 cd ..
echo "$n"
