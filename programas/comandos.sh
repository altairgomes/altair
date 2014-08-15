i=$(cat lista_pastas2);k=`pwd`; for j in $i;do cd $k/$j; ./PRAIA_header_extraction_20_05 < PRAIA_header_extraction_20_05.dat;done

i=$(cat lista_pastas); for j in $i;do if [ -e $j/output ];then s=0; else echo "$j" >> lista_pastas2;fi;done

wget --progress=bar:force "http://www.campinas.sp.gov.br/uploads/pdf/371642022.pdf" -O"teste" 2>&1 | zenity --title="teste" --progress

while read COL_1 COL_2 COL_3
do
    <some code using your read variables COL_1, COL_2, COL_3>
done < table.dat

let C="$A + $B"

C=$(echo "34 * 45"|bc)

(
echo "10" ; sleep 1
echo "# Updating mail logs" ; sleep 1
echo "20" ; sleep 1
echo "# Resetting cron jobs" ; sleep 1
echo "50" ; sleep 1
echo "This line will just be ignored" ; sleep 1
echo "75" ; sleep 1
echo "# Rebooting system" ; sleep 1
echo "100" ; sleep 1
) |
zenity --progress \
  --title="Update System Logs" \
  --text="Scanning mail logs..." \
  --percentage=0


if [ $(echo "$dist > $limitedist" | bc) -ne 0 ] && [ $(echo "$elev > $limiteelev" | bc) -ne 0 ];then
   grep "${linha[2]}" $entrada >> $saida1
 else
   grep "${linha[2]}" $entrada >> $saida2
fi
