inicial=2448622.5
final=2456293.5
step=1.0
soma=2448622.5
while [ "$soma" != "$final" ] ;do
 echo "$soma JD" >> julian_date_orbit
 soma=$(echo "$soma + $step" | bc -l)
done
