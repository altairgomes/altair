#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"
#
# Inserir Objetos
#
insert_object(){
objects=$(yad --form \
  --image "accessories-text-editor" \
  --title "Inserir objetos" \
  --width=400 \
  --height=200 \
  --field="Nome do Objeto" "" \
  --field="Ascencao reta" "" \
  --field="Declinacao" "" \
  --field="Informações" "" \
  )
}
#
# Ler Objetos
#
read_object(){
arqalvos2="alvos2"
(while read nome afa afb afc dca dcb dcc inf;do
  objeto+=("$nome")
  ascencao+=("$afa $afb $afc")
  declinacao+=(" $dca $dcb $dcc")
  objinf+=("$inf")
  nastes+=("0")
  poetes+=("0")
done < $arqalvos
i=0
while [ $i -lt "${#objeto[@]}" ];do
  grupoo+="${objeto[$i]}+"
  read rah rab rac <<< ${ascencao[$i]}
  grupoar+=(" ${ascencao[$i]}")
  grupoarr+=(" $(calc "($rah + ($rab/60) + ($rac/3600))*15*($pi/180)" < /dev/null | sed 's/\~//')")
  read dea deb dec <<< ${declinacao[$i]}
  grupode+=(" ${declinacao[$i]}")
  grupo+=("\n\t\t${objeto[$i]} (${objinf[$i]})\n\t\t  RA: $rah $rab $rac; DEC: $dea $deb $dec;")
  if [ $(echo " $dea < 0" | bc) -ne 0 -o " $dea" = " -0" -o " $dea" = " -00" ];then
    decx=$(calc " $dea - ($deb/60) - ($dec/3600)" < /dev/null | sed 's/\~//')
    grupoder+=("$(calc " $decx*($pi/180)" < /dev/null | sed 's/\~//')")
  else
    decx=$(calc " $dea + ($deb/60) + ($dec/3600)" < /dev/null | sed 's/\~//')
    grupoder+=("$(calc " $decx*($pi/180)" < /dev/null | sed 's/\~//')")
  fi
  j=$[$i+1]
  distmax="0"
  alfamed="0"
  decmed="0"
  grupovar="0"
  while [ $j -lt "${#objeto[@]}" ];do
    read rabh rabb rabc <<< ${ascencao[$j]}
    aretab=$(calc "($rabh + ($rabb/60) + ($rabc/3600))*15*($pi/180)" < /dev/null | sed 's/\~//')
    read deba debb debc <<< ${declinacao[$j]}
    if [ $(echo " $deba < 0" | bc) -ne 0 -o " $deba" = " -0" -o " $deba" = " -00" ];then
      decy=$(calc " $deba - ($debb/60) - ($debc/3600)" < /dev/null | sed 's/\~//')
      decb=("$(calc " $decy*($pi/180)" < /dev/null | sed 's/\~//')")
    else
      decy=$(calc " $deba + ($debb/60) + ($debc/3600)" < /dev/null | sed 's/\~//')
      decb=("$(calc " $decy*($pi/180)" < /dev/null | sed 's/\~//')")
    fi
    k=0
    fazer="1"
    while [ $k -lt "${#grupo[@]}" ];do
      alfadif=$(calc "$aretab-${grupoarr[$k]}" < /dev/null | sed 's/\~//')
      distang=$(calc "acos(sin($decb)*sin(${grupoder[$k]})+cos($decb)*cos(${grupoder[$k]})*cos($alfadif))" < /dev/null | sed 's/\~//')
      distang=$(calc "$distang*(180/$pi)*60" < /dev/null | sed 's/\~//')
      if [ $(echo "$distang > $limdist" | bc) -ne 0 -a "$fazer" = "1" ];then
        fazer="0"
      fi
      if [ $(echo "$distang < $limdist" | bc) -ne 0 -a $(echo "$distang > $distmax" | bc) -ne 0 -a "$fazer" = "1" ];then
        distmax="$distang"
        alfamed=$(calc "($aretab+${grupoarr[$k]})/2" < /dev/null | sed 's/\~//')
        decmed=$(calc "($decb+${grupoder[$k]})/2" < /dev/null | sed 's/\~//')
      fi
      k=$[$k+1]
    done
    if [ $fazer = "1" ];then
      echo "${grupo[@]}" | egrep "${objeto[$j]}" > /dev/null
      if [ $? -ne 0 ];then
        grupo+=("\n\t\t${objeto[$j]} (${objinf[$j]})\n\t\t  RA: $rabh $rabb $rabc; DEC: $deba $debb $debc;")
        grupoo+="${objeto[$j]}+"
        grupoar+=("${ascencao[$j]}")
        grupoarr+=("$aretab")
        grupode+=(" ${declinacao[$j]}")
        grupoder+=(" $decb")
        grupovar="1"
      fi
    fi
    j=$[$j+1]
  done
  alfahmed=$(calc "($alfamed*(180/$pi))/15" < /dev/null | sed 's/\~//')
  read alfamedh <<< $(calc "int($alfahmed)" | sed 's/\~//')
  read alfamedm <<< $(calc "int(($alfahmed-$alfamedh)*60)" | sed 's/\~//')
  read alfameds <<< $(printf "%5.2f" $(calc "(($alfahmed-$alfamedh)*60-$alfamedm)*60" | sed 's/\~//'))
  dechmed=$(calc " $decmed*(180/$pi)" < /dev/null | sed 's/\~//')
  read decmedh <<< $(calc "int($dechmed)" | sed 's/\~//')
  read decmedm <<< $(calc "abs(int(($dechmed - $decmedh)*60))" | sed 's/\~//')
  read decmeds <<< $(printf "%5.2f" $(calc "((abs($dechmed - $decmedh)*60)-$decmedm)*60" | sed 's/\~//'))
  teste=$(cat .aux1 | grep "$grupoo ")
  if [  "$teste" != ""  ];then
    grupovar="0"
  fi
  if [ "$grupovar" = "1" ];then
    echo "$grupoo $alfamedh $alfamedm $alfameds $decmedh $decmedm $decmeds ${grupo[@]}" >> $arqalvos2
  else
    teste2=$(cat .aux1 | grep "${objeto[$i]}+")
    if [  "$teste2" = ""  ];then
      echo "${objeto[$i]} ${ascencao[$i]} ${declinacao[$i]} ${objinf[$i]}" >> $arqalvos2
    fi
  fi
  echo "$grupoo" >> .aux1
  i=$[$i+1]
  read per <<< $(calc "($i*100)/${#objeto[@]}" | sed 's/\~//')
  echo "$per"
  echo "#$(printf "%5.2f" $per)%"
  unset grupo; unset grupoo; unset grupoar; unset grupoarr; unset grupode; unset grupoder
done
unset objeto; unset ascencao; unset declinacao; unset objinf; unset nascer; unset poente; unset nastes; unset poetes;
) | yad --progress \
  --width=300\
  --title="Lendo Objetos" \
  --text="Lendo Objetos" \
  --progress-text=0.00% \
  --percentage=0 \
  --auto-close
while read -r nome afa afb afc dca dcb dcc inf;do
  objeto+=("$nome")
  ascencao+=("$afa $afb $afc")
  declinacao+=(" $dca $dcb $dcc")
  objinf+=("$inf")
  nastes+=("0")
  poetes+=("0")
done < $arqalvos2
rm $arqalvos2
rm .aux1
}
#
# Gerar arquivo de planejamento
#
gerar_plan(){
read_object
(arquivo="plano_$datein"
read lta ltb ltc <<< $latitude
read lga lgb lgc <<< $longitude
if [ $(echo " $lta < 0" | bc) -ne 0 -o " $lta" = " -0" -o " $lta" = " -00" ];then
  latdeg=$(calc " $lta - ($ltb/60) - ($ltc/3600)" | sed 's/\~//')
else
  latdeg=$(calc " $lta + ($ltb/60) + ($ltc/3600)" | sed 's/\~//')
fi
if [ $(echo " $lga < 0" | bc) -ne 0 -o " $lga" = " -0" -o " $lga" = " -00" ];then
  londeg=$(calc " $lga - ($lgb/60) - ($lgc/3600)" | sed 's/\~//')
else
  londeg=$(calc " $lga + ($lgb/60) + ($lgc/3600)" | sed 's/\~//')
fi
if [ $(echo " $londeg < 0" | bc) ];then londeg=$(calc " $londeg + 360" | sed 's/\~//');fi
latit=$(calc " $latdeg*($pi/180)" | sed 's/\~//')
lonhor=$(calc " $londeg/15" | sed 's/\~//')
diacomeco=$(date --date=$datein +%_j)
diafinal=$(date --date=$datefi +%_j)
tempoi=$(calc "$horain*60" | sed 's/\~//')
tempof=$(calc "24*60*($diafinal - $diacomeco) + $horafi*60" | sed 's/\~//')
dif=$(calc "$tempof-$tempoi" | sed 's/\~//')
int=$intervalo
min=0
hor=$(calc "$horain - $fuso" | sed 's/\~//')
dia=$diacomeco
tempo=$tempoi
ano=$(date --date=$datein +%Y)
mes=$(date --date=$datein +%_j)
tau=$(calc "$ano+$mes/365.25-2000" | sed 's/\~//')
tsec=$(calc "($ano+$mes/365.25-2000)/100" | sed 's/\~//')
presm=$(calc "3.07496+0.00186*$tsec" | sed 's/\~//')
presn=$(calc "1.33621-0.00057*$tsec" | sed 's/\~//')
presndeg=$(calc "20.0431-0.0085*$tsec" | sed 's/\~//')
echo -e "Plano de observação da noite $datein\n" > $arquivo
echo -e "Latitude: $latitude  Longitude: $longitude\nLimite de altura: $altlimº\nTamanho do campo: $limdist arcmin\n\n" >> $arquivo
while [ $tempo -le $tempof ];do
  if [ $(echo "$min > 59" | bc) -ne 0 ]; then
    min=$(calc "$min-60" | sed 's/\~//')
    hor=$(calc "$hor+1" | sed 's/\~//')
  fi
  if [ $(echo "$hor > 23" | bc) -ne 0 ]; then
    hor=$(calc "$hor-24" | sed 's/\~//')
    dia=$(calc "$dia+1" | sed 's/\~//')
  fi
  hora=$(calc "$hor + $min/60" | sed 's/\~//')
  tsl=$(calc "6.654 + 0.06570982441908*$dia + 1.00273791*$hora + $lonhor" | sed 's/\~//')
  tsl=$(calc " $tsl%24" | sed 's/\~//')
  hora2=$(calc "$hor + $fuso" < /dev/null | sed 's/\~//')
  if [ $(echo "$hora2 < 0" | bc) -ne 0 ];then hora2=$(calc "$hora2+24" < /dev/null | sed 's/\~//');fi
  if [ $(echo "$hora2 > 23" | bc) -ne 0 ];then hora2=$(calc "$hora2-24" < /dev/null | sed 's/\~//');fi
  read hour <<< $hor
  read hora2 <<< $hora2
  read minuto <<< $min
  i=0
  for item in "${objeto[@]}"; do
    if [ "${objeto[$i]}" != "" ];then
    read afa afb afc <<< ${ascencao[$i]}
    read dca dcb dcc <<< ${declinacao[$i]}
    afhor=$(calc "$afa + ($afb/60) + ($afc/3600)" < /dev/null | sed 's/\~//')
    afdeg=$(calc "$afhor*15" < /dev/null | sed 's/\~//')
    afrad=$(calc "$afdeg*($pi/180)" < /dev/null | sed 's/\~//')
    if [ $(echo " $dca < 0" | bc) -ne 0 -o " $dca" = " -0" -o " $dca" = " -00" ];then
      decdeg=$(calc " $dca - ($dcb/60) - ($dcc/3600)" < /dev/null | sed 's/\~//')
    else
      decdeg=$(calc " $dca + ($dcb/60) + ($dcc/3600)" < /dev/null | sed 's/\~//')
    fi
    dec=$(calc " $decdeg*($pi/180)" < /dev/null | sed 's/\~//')
    dalfa=$(calc "($presm+$presn*sin($afrad)*tan($dec))*$tau" < /dev/null | sed 's/\~//')
    ddelta=$(calc "$tau*($presndeg)*cos($afrad)" < /dev/null | sed 's/\~//')
    alfapre=$(calc "$afhor+($dalfa/3600)" < /dev/null | sed 's/\~//')
    decpre=$(calc " $dec+($ddelta/3600)*($pi/180)" < /dev/null | sed 's/\~//')
    hourangle=$(calc "$tsl - $alfapre" < /dev/null | sed 's/\~//')
    houranglerad=$(calc "$hourangle*15*($pi/180)" < /dev/null | sed 's/\~//')
    distzenrad=$(calc "acos(sin($decpre)*sin($latit)+cos($decpre)*cos($latit)*cos($houranglerad))" < /dev/null | sed 's/\~//')
    distzen=$(calc "$distzenrad*(180/$pi)" < /dev/null | sed 's/\~//')
    altura=$(calc "90-$distzen" < /dev/null | sed 's/\~//')
    read height <<< $altura
    height[$i]=$(echo ${height::5})
    if [ $(echo "$altura > $altlim" | bc) -ne 0 ];then
      coorda[$i]=$(printf "%3s%3s %4s" $afa $afb $afc)
      coordd[$i]=$(printf "%3s%3s %4s" $dca $dcb $dcc)
      maxalt=$(calc "($alfapre - 6.654 - 0.06570982441908*$dia - $lonhor)/1.00273791 + $fuso" | sed 's/\~//')
      maxalt=$(calc " $maxalt%24" | sed 's/\~//')
      if [ $(echo " $maxalt > 23" | bc) -ne 0 ]; then
        maxalt=$(calc " $maxalt - 24" | sed 's/\~//')
      elif [ $(echo " $maxalt < 0" | bc) -ne 0 ]; then
        maxalt=$(calc " $maxalt + 24" | sed 's/\~//')
      fi
      read maxalth <<< $(calc "int($maxalt)" | sed 's/\~//')
      read maxaltm <<< $(calc "int(($maxalt-$maxalth)*60)" | sed 's/\~//')
      read maxalts <<< $(calc "int((($maxalt-$maxalth)*60-$maxaltm)*60)" | sed 's/\~//')
      maxaltob[$i]="$maxalth:$maxaltm:$maxalts"
      maxdist=$(calc "(90-$altlim)*($pi/180)" < /dev/null | sed 's/\~//')
      anghorrad=$(calc "acos(cos($maxdist)*sec($latit)*sec($decpre) - tan($latit)*tan($decpre))" < /dev/null | sed 's/\~//')
      anghor=$(calc "($anghorrad/15)*(180/$pi)" < /dev/null | sed 's/\~//')
      limhor=$(calc "$anghor - $hourangle" < /dev/null | sed 's/\~//')
      limhor=$(calc " $limhor%24" | sed 's/\~//')
      if [ $(echo " $limhor > 23" | bc) -ne 0 ]; then
        limhor=$(calc " $limhor - 24" | sed 's/\~//')
      elif [ $(echo " $limhor < 0" | bc) -ne 0 ]; then
        limhor=$(calc " $limhor + 24" | sed 's/\~//')
      fi
      read limhorh <<< $(calc "int($limhor)" | sed 's/\~//')
      read limhorm <<< $(calc "int(($limhor-$limhorh)*60)" | sed 's/\~//')
      read limhors <<< $(calc "int((($limhor-$limhorh)*60-$limhorm)*60)" | sed 's/\~//')
      echo "${objeto[$i]} $limhorh:$limhorm:$limhors" >> .aux0
    fi
    if [ $(echo "${nastes[$i]} < $altlim" | bc) -ne 0 -a $(echo "$altura > $altlim" | bc) -ne 0 ];then
      nastes[$i]="${height[$i]}"
      nascer[$i]="$hora2:$minuto(${height[$i]})"
    fi
    if [ $(echo "${nastes[$i]} > $altlim" | bc) -ne 0 -a $(echo "$altura > $altlim" | bc) -ne 0 ];then
      poente[$i]="$hora2:$minuto(${height[$i]})"
    fi
    let i++
    fi
  done
  lista=$(sort -n -k 2.1,2.2 -k 2.4,2.5 -k 2.7,2.8 .aux0)
  echo "$lista" > .aux0
  while read object texto;do
    i=0
    for item in "${objeto[@]}"; do
      if [ "$object" = "$item" ];then
        if [ "$(echo "${objinf[$i]}" | egrep "DEC")" != "" ];then
          target+=("${objeto[$i]%+*}\n\tRAmedio:${coorda[$i]}\tDECmedio: ${coordd[$i]}\n\tAltura: ${height[$i]}\n\tPassagem Meridiana: ${maxaltob[$i]} TL\n\tTempo restante para atingir limite: $texto ${objinf[$i]}")
        else
          target+=("${objeto[$i]}  (${objinf[$i]})\n\tRA:${coorda[$i]}\tDEC: ${coordd[$i]}\n\tAltura: ${height[$i]}\n\tPassagem Meridiana: ${maxaltob[$i]} TL\n\tTempo restante para atingir limite: $texto")
        fi
      fi
      let i++
    done
  done < .aux0
  rm .aux0
  echo -e "\n--- TL=$hora2:$minuto (UT=$hour:$minuto) ------------------------------------------------------------------------------------------------------------------------------------" >> $arquivo
  for x in "${target[@]}";do
    echo -e "$x" >> $arquivo
  done
  unset target
  tempo=$(calc "$tempo + $int" < /dev/null | sed 's/\~//')
  read porce <<< $(calc "($tempo-$tempoi)*100/$dif" < /dev/null | sed 's/\~//')
  echo "$porce"
  porce=$(printf "%5.2f" $porce)
  echo "#Gerando $porce%"
  min=$(calc "$min + $int" < /dev/null | sed 's/\~//')
done
echo -e "\n\n\n--------------------------------------------------------------------------------------------------------------------------------------------------------------\n" >> $arquivo
echo -e "Visibilidade dos alvos\n" >> $arquivo
i=0
for item in "${objeto[@]}"; do
  rightas=$(printf "%3s%3s%8s" ${ascencao[$i]})
  decs=$(printf "%3s%3s%7s" ${declinacao[$i]})
  printf "%30s (A.R.:$rightas; DEC:$decs): inicio=${nascer[$i]} --> fim=${poente[$i]}\n" ${objeto[$i]}  >> $arquivo
  let i++
done) | yad --progress \
--width=400 \
--title="Plano de Observacao..." \
--text="Gerando Planejamento..." \
--progress-text="Gerando 0.00%"   \
--percentage=0 #--auto-close
unset objeto; unset ascencao; unset declinacao; unset objinf; unset nascer; unset poente; unset nastes; unset poetes;
}
## funcao para calcular pra hora atual
dateatual(){
read_object
read lta ltb ltc <<< $latitude
read lga lgb lgc <<< $longitude
if [ $(echo " $lta < 0" | bc) -ne 0 -o " $lta" = " -0" -o " $lta" = " -00" ];then
  latdeg=$(calc " $lta - ($ltb/60) - ($ltc/3600)" | sed 's/\~//')
else
  latdeg=$(calc " $lta + ($ltb/60) + ($ltc/3600)" | sed 's/\~//')
fi
if [ $(echo " $lga < 0" | bc) -ne 0 -o " $lga" = " -0" -o " $lga" = " -00" ];then
  londeg=$(calc " $lga - ($lgb/60) - ($lgc/3600)" | sed 's/\~//')
else
  londeg=$(calc " $lga + ($lgb/60) + ($lgc/3600)" | sed 's/\~//')
fi
if [ $(echo " $londeg < 0" | bc) ];then londeg=$(calc " $londeg + 360" | sed 's/\~//');fi
latit=$(calc " $latdeg*($pi/180)" | sed 's/\~//')
lonhor=$(calc " $londeg/15" | sed 's/\~//')
horain=$(date +%_H)
tempof="25"
ano=$(date +%Y)
mes=$(date +%_m)
tau=$(calc "$ano+$mes/365.25-2000" | sed 's/\~//')
tsec=$(calc "($ano+$mes/365.25-2000)/100" | sed 's/\~//')
presm=$(calc "3.07496+0.00186*$tsec" | sed 's/\~//')
presn=$(calc "1.33621-0.00057*$tsec" | sed 's/\~//')
presndeg=$(calc "20.0431-0.0085*$tsec" | sed 's/\~//')
arquivo2=".listaatual"
cancelar="0"
while [ "$cancelar" != "1" ];do
  horain=$(date +%_H)
  min=$(date +%_M)
  seg=$(date +%_S)
  hor=$(calc "$horain - $fuso" | sed 's/\~//')
  dia=$(date +%_j)
  data=$(date +%F)
  if [ $(echo "$hor > 23" | bc) -ne 0 ]; then
    hor=$(calc "$hor-24" | sed 's/\~//')
    dia=$(calc "$dia+1" | sed 's/\~//')
  fi
  hora=$(calc " $hor + $min/60 + $seg/3600" | sed 's/\~//')
  tsl=$(calc "6.654 + 0.06570982441908*$dia + 1.00273791*$hora + $lonhor" | sed 's/\~//')
  tsl=$(calc " $tsl%24" | sed 's/\~//')
  i=0
  for item in "${objeto[@]}"; do
    if [ "${objeto[$i]}" != "" ];then
    read afa afb afc <<< ${ascencao[$i]}
    read dca dcb dcc <<< ${declinacao[$i]}
    afhor=$(calc "$afa + ($afb/60) + ($afc/3600)" < /dev/null | sed 's/\~//')
    afdeg=$(calc "$afhor*15" < /dev/null | sed 's/\~//')
    afrad=$(calc "$afdeg*($pi/180)" < /dev/null | sed 's/\~//')
    if [ $(echo " $dca < 0" | bc) -ne 0 -o " $dca" = " -0" -o " $dca" = " -00" ];then
      decdeg=$(calc " $dca - ($dcb/60) - ($dcc/3600)" < /dev/null | sed 's/\~//')
    else
      decdeg=$(calc " $dca + ($dcb/60) + ($dcc/3600)" < /dev/null | sed 's/\~//')
    fi
    dec=$(calc " $decdeg*($pi/180)" < /dev/null | sed 's/\~//')
    dalfa=$(calc "($presm+$presn*sin($afrad)*tan($dec))*$tau" < /dev/null | sed 's/\~//')
    ddelta=$(calc "$tau*($presndeg)*cos($afrad)" < /dev/null | sed 's/\~//')
    alfapre=$(calc "$afhor+($dalfa/3600)" < /dev/null | sed 's/\~//')
    decpre=$(calc " $dec+($ddelta/3600)*($pi/180)" < /dev/null | sed 's/\~//')
    hourangle=$(calc "$tsl - $alfapre" < /dev/null | sed 's/\~//')
    houranglerad=$(calc "$hourangle*15*($pi/180)" < /dev/null | sed 's/\~//')
    distzenrad=$(calc "acos(sin($decpre)*sin($latit)+cos($decpre)*cos($latit)*cos($houranglerad))" < /dev/null | sed 's/\~//')
    distzen=$(calc "$distzenrad*(180/$pi)" < /dev/null | sed 's/\~//')
    altura=$(calc "90-$distzen" < /dev/null | sed 's/\~//')
    read height <<< $altura
    height[$i]=$(echo ${height::5})
    if [ $(echo "$altura > $altlim" | bc) -ne 0 ];then
      coorda[$i]=$(printf "%3s%3s %4s" $afa $afb $afc)
      coordd[$i]=$(printf "%3s%3s %4s" $dca $dcb $dcc)
      maxalt=$(calc "($alfapre - 6.654 - 0.06570982441908*$dia - $lonhor)/1.00273791 + $fuso" | sed 's/\~//')
      maxalt=$(calc " $maxalt%24" | sed 's/\~//')
      if [ $(echo " $maxalt > 23" | bc) -ne 0 ]; then
        maxalt=$(calc " $maxalt - 24" | sed 's/\~//')
      elif [ $(echo " $maxalt < 0" | bc) -ne 0 ]; then
        maxalt=$(calc " $maxalt + 24" | sed 's/\~//')
      fi
      read maxalth <<< $(calc "int($maxalt)" | sed 's/\~//')
      read maxaltm <<< $(calc "int(($maxalt-$maxalth)*60)" | sed 's/\~//')
      read maxalts <<< $(calc "int((($maxalt-$maxalth)*60-$maxaltm)*60)" | sed 's/\~//')
      maxaltob[$i]="$maxalth:$maxaltm:$maxalts"
      maxdist=$(calc "(90-$altlim)*($pi/180)" < /dev/null | sed 's/\~//')
      anghorrad=$(calc "acos(cos($maxdist)*sec($latit)*sec($decpre) - tan($latit)*tan($decpre))" < /dev/null | sed 's/\~//')
      anghor=$(calc "($anghorrad/15)*(180/$pi)" < /dev/null | sed 's/\~//')
      limhor=$(calc "$anghor - $hourangle" < /dev/null | sed 's/\~//')
      limhor=$(calc " $limhor%24" | sed 's/\~//')
      if [ $(echo " $limhor > 23" | bc) -ne 0 ]; then
        limhor=$(calc " $limhor - 24" | sed 's/\~//')
      elif [ $(echo " $limhor < 0" | bc) -ne 0 ]; then
        limhor=$(calc " $limhor + 24" | sed 's/\~//')
      fi
      read limhorh <<< $(calc "int($limhor)" | sed 's/\~//')
      read limhorm <<< $(calc "int(($limhor-$limhorh)*60)" | sed 's/\~//')
      read limhors <<< $(calc "int((($limhor-$limhorh)*60-$limhorm)*60)" | sed 's/\~//')
      echo "${objeto[$i]} $limhorh:$limhorm:$limhors" >> .aux0
    fi
    let i++
    fi
  done
  read tslh <<< $(calc "int($tsl)" | sed 's/\~//')
  read tslm <<< $(calc "int(($tsl-$tslh)*60)" | sed 's/\~//')
  read tsls <<< $(calc "int((($tsl-$tslh)*60-$tslm)*60)" | sed 's/\~//')
  read fus <<< $(calc "int($fuso)" | sed 's/\~//')
  lista=$(sort -n -k 2.1,2.2 -k 2.4,2.5 -k 2.7,2.8 .aux0)
  echo "$lista" > .aux0
  while read object texto;do
    i=0
    for item in "${objeto[@]}"; do
      if [ "$object" = "$item" ];then
        if [ "$(echo "${objinf[$i]}" | egrep "DEC")" != "" ];then
          target+=("${objeto[$i]%+*}\n\tRAmedio:${coorda[$i]}\tDECmedio: ${coordd[$i]}\n\tAltura atual: ${height[$i]}\n\tPassagem Meridiana: ${maxaltob[$i]} TL\n\tTempo restante para atingir limite: $texto ${objinf[$i]}")
        else
          target+=("${objeto[$i]}  (${objinf[$i]})\n\tRA:${coorda[$i]}\tDEC: ${coordd[$i]}\n\tAltura atual: ${height[$i]}\n\tPassagem Meridiana: ${maxaltob[$i]} TL\n\tTempo restante para atingir limite: $texto")
        fi
      fi
      let i++
    done
  done < .aux0
  rm .aux0
  echo -e "\nDATA: $data  HORA: $horain:$min:$seg   Fuso: $fus\nTempo Sideral Local: $tslh:$tslm:$tsls" > $arquivo2
  echo -e "\nLatitude: $latitude  Longitude: $longitude\nQuantidade de objetos observaveis: ${#target[@]}\nLimite de altura: $altlimº\n\n" >> $arquivo2
  for x in "${target[@]}";do
    echo -e "$x\n" >> $arquivo2
  done
  yad --height=600 --width=600 --text-info --filename="$arquivo2" --button="Sair":1 --button="Atualizar":0
  ret=$?
  if [ "$ret" = "1" ];then cancelar="1";fi
  unset target
done
i=0
unset objeto; unset ascencao; unset declinacao; unset objinf; unset nascer; unset poente; unset nastes; unset poetes;
}
## funcao inicial
cancel=0
n=0
arqalvos="alvos"
pi=$(calc "4*atan(1)" | sed 's/\~//')
until [ "$cancel" = "1" ];do
  count=$(wc -l < $arqalvos)
  alvos=$(sort -n --key=2 --key=3 $arqalvos)
  datas=$(yad --form \
  --image "accessories-text-editor" \
  --title "Planejar observações" \
  --width=500 \
  --height=600 \
  --button=Sair:1 --button="Gerar planejamento":2 \
  --text="Já existem $count objeto(s) inserido(s)" \
  --field="Data inicial de observação":DT "$datein" \
  --field="Hora inicial da observação (TL)":NUM 18!12..24!1 \
  --field="Data do fim de observação":DT "$datefi" \
  --field="Hora final da observação (TL)":NUM 6!0..12!1 \
  --field="Ver informações a cada x minutos":NUM 30!10..120!10 \
  --field="Calcular pra hora atual":CHK FALSE \
  --field="Longitude" "314 25 2.5" \
  --field="Latitude" " -22 32 7.8" \
  --field="Fuso":NUM [0[!-12..12[!1]]] \
  --field="Limite de altura":NUM 30!0..90!5 \
  --field="Limite de distancia dos objetos (arcmin)":NUM 10!1..15!1 \
  --field="Objetos":TXT "$alvos" \
  --date-format="%F" \
  )
  res=$?
  if [ ! -z "$datas" ];then
    datein=$(echo "$datas" | cut -d"|" -f 1)
    horain=$(echo "$datas" | cut -d"|" -f 2)
    horain=$(echo $horain | tr ',' '.')
    datefi=$(echo "$datas" | cut -d"|" -f 3)
    horafi=$(echo "$datas" | cut -d"|" -f 4)
    horafi=$(echo $horafi | tr ',' '.')
    intervalo=$(echo "$datas" | cut -d"|" -f 5)
    intervalo=$(echo $intervalo | tr ',' '.')
    atual=$(echo "$datas" | cut -d"|" -f 6)
    longitude=$(echo "$datas" | cut -d"|" -f 7)
    latitude=$(echo "$datas" | cut -d"|" -f 8)
    fuso=$(echo "$datas" | cut -d"|" -f 9)
    altlim=$(echo "$datas" | cut -d"|" -f 10)
    read altlim <<< $(calc "int($altlim)" | sed 's/\~//')
    limdist=$(echo "$datas" | cut -d"|" -f 11)
    read limdist <<< $(calc "int($limdist)" | sed 's/\~//')
    alvos=$(echo "$datas" | cut -d"|" -f 12)
    echo -e "$alvos" > $arqalvos
    alvos=$(sort -n --key=2 --key=3 --key=4 $arqalvos)
    echo -e "$alvos" > $arqalvos
    if [ "$atual" = "TRUE" ];then res="3";fi
  fi
  case $res in
    1)cancel=1;;
    2)cancel=0; gerar_plan;;
    3)cancel=0; dateatual;;
  esac
done
