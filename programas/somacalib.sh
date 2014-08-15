entrada="photometry_Chariklo.dat"
saida="photometry_Chariklo_med.dat"
count=0
bin=$1
num=$(echo "($bin - 1)/2" | bc)
while read a b c d e f g h i j k;do
  ba[$count]="$a"
  bb[$count]="$b"
  bc[$count]="$c"
  bd[$count]="$d"
  be[$count]="$e"
  bf[$count]="$f"
  bg[$count]="$g"
  bh[$count]="$h"
  bi[$count]="$i"
  bj[$count]="$j"
  bk[$count]="$k"
  count=$(echo "$count + 1" | bc)
done < $entrada
quant=${#ba[@]}
init=$(echo "$num" | bc)
fini=$(echo "$quant - $num - 1" | bc)
n=$init
while [ $n -le $fini ];do
  echo $n
  come=$(echo "$n - $num" | bc)
  term=$(echo "$n + $num" | bc)
  x=$come
  somacalib=0.0
  while [ $x -le $term ];do
    somacalib=$(calc "$somacalib + ${bh[$x]}" < /dev/null | sed 's/\~//')
    x=$(echo "$x + 1" | bc)
  done
  calibmed=$(calc "$somacalib / $bin" < /dev/null | sed 's/\~//')
  printf " %3d %3d %16.8f %5.3f %16.5f %10.4f %6.3f %16.5f %10.4f %6.3f %s\n" ${ba[$n]} ${bb[$n]} ${bc[$n]} ${bd[$n]} ${be[$n]} ${bf[$n]} ${bg[$n]} $calibmed ${bi[$n]} ${bj[$n]} ${bk[$n]} >> $saida
  n=$(echo "$n + 1" | bc)
done
