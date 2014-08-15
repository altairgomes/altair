entrada="curve_chariklo_29abr2014.dat"
saida="curve_5pts"
binagem=$1
count=0
somac="0.0"
somad="0.0"
somaj="0.0"
somal="0.0"
somam="0.0"
while read a b c d e f g h i j k l m;do
  count=$(echo "$count + 1" | bc)
  somac=$(calc "$somac + $c" < /dev/null | sed 's/\~//')
  somad=$(calc "$somad + $d" < /dev/null | sed 's/\~//')
  somaj=$(calc "$somaj + $j" < /dev/null | sed 's/\~//')
  somal=$(calc "$somal + $l" < /dev/null | sed 's/\~//')
  somam=$(calc "$somam + $m" < /dev/null | sed 's/\~//')
  if [ "$count" = "$binagem" ];then
    saidac=$(calc "$somac / $count" < /dev/null | sed 's/\~//')
    saidad=$(calc "$somad / $count" < /dev/null | sed 's/\~//')
    saidaj=$(calc "$somaj / $count" < /dev/null | sed 's/\~//')
    saidal=$(calc "$somal / $count" < /dev/null | sed 's/\~//')
    saidam=$(calc "$somam / $count" < /dev/null | sed 's/\~//')
    printf " %3d %5.3f %16.8f %16.8f %16.14f %6.3f %6.3f %9.5f %9.5f %16.11f %16.11f %16.5f %16.5f\n" $a $b $saidac $saidad $e $f $g $h $i $saidaj $k $saidal $saidam >> $saida
    count=0
    somac=0.0
    somad=0.0
    somaj=0.0
    somal=0.0
    somam=0.0
  fi 
done < $entrada 
