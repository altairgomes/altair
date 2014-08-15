atual=$(pwd)
objects="$atual/.objects"                    #arquivo com os dados dos objeto
calibracao="$atual/.calibra"                 #arquivo com os dados das calibrações
categoriax="$atual/.categoriax"              #arquivo com os dados das categorias
telescopes="$atual/.telescopes"              #arquivo com os dados dos telescopio
cameras="$atual/.cameras"                    #arquivo com os dados das cameras
binagens="$atual/.binagens"                  #arquivo com os dados das binagens
filtros="$atual/.filtros"                    #arquivo com os dados dos filtros
calendario="$atual/.datas"                   #arquivo com os dados das datas
pessoas="$atual/.observadores"               #arquivo com os dados das pessoas
addreadme="$atual/readme_adiciona"           #arquivo readme da funcao adiciona
copiareadme="$atual/readme_copia"            #arquivo readme da funcao copia
anareadme="$atual/readme_analise"            #arquivo readme da funcao analise
redreadme="$atual/readme_reducao"            #arquivo readme da funcao reducao

export LC_NUMERIC="en_US.UTF-8"
export LC_TIME="en_US.UTF-8"


############### Leitura dos arquivos de dados #######################################################################################
read_categoria(){
if [ -f $categoriax ];then
  cate=$(sort $categoriax)
  for j in $cate;do
    echo -n "$j!" >> .aux11
    echo -n "FALSE $j " >> .aux12
  done
  VAR1=$(cat .aux11)
  opcaocategoria=$(echo ${VAR1%!*})
  rm .aux11
  opcaocategoria2=$(cat .aux12)
  rm .aux12
fi
}
read_object(){
if [ -f $objects ];then
  sort $objects > .aux2
  while read COL_1 COL_2 COL_3;do
    echo -n "!$COL_1" >> .aux21
    echo -n "FALSE $COL_1 " >> .aux22
    echo "$COL_1" >> .aux23
    if [ -e .aux24 ];then cat .aux24 | egrep "$COL_1" > /dev/null;fi
    if [ "${COL_1:0:4}" != "Est_" ];then
      echo -n "!$COL_1" >> .aux24
    fi
  done < .aux2
  rm .aux2
  opcaoobjeto2=$(cat .aux22)
  rm .aux22
  listobjeto=$(cat .aux23)
  rm .aux23
  opcaoobjeto3=$(cat .aux24)
  rm .aux24
fi
echo -n "!BIAS!FLAT" >> .aux21
opcaoobjeto=$(cat .aux21)
rm .aux21
}
read_telescope(){
if [ -f $telescopes ];then
  teles=$(sort $telescopes)
  for j in $teles;do
    echo -n "!$j" >> .aux31
    echo -n "FALSE $j " >> .aux32
  done
  opcaotelescopio=$(cat .aux31)
  rm .aux31
  opcaotelescopio2=$(cat .aux32)
  rm .aux32
fi
}
read_camera(){
if [ -f $cameras ];then
  came=$(sort $cameras)
  for j in $came;do
    echo -n "!$j" >> .aux41
    echo -n "FALSE $j " >> .aux42
  done
  opcaocamera=$(cat .aux41)
  rm .aux41
  opcaocamera2=$(cat .aux42)
  rm .aux42
fi
}
read_binagem(){
if [ -f $binagens ];then
  bina=$(sort $binagens)
  for j in $bina;do
    echo -n "!$j" >> .aux51
    echo -n "FALSE $j " >> .aux52
  done
  opcaobinagem=$(cat .aux51)
  rm .aux51
  opcaobinagem2=$(cat .aux52)
  rm .aux52
fi
}
read_filtro(){
if [ -f $filtros ];then
  filt=$(sort $filtros)
  for j in $filt;do
    echo -n "!$j" >> .aux61
    echo -n "FALSE $j " >> .aux62
  done
  opcaofiltro=$(cat .aux61)
  rm .aux61
  opcaofiltro2=$(cat .aux62)
  rm .aux62
fi
}
read_data(){
if [ -f $calendario ];then
  dat=$(sort -k 1.1,1.4 -k 1.6,1.8M -k 1.10,1.11 $calendario)
  for j in $dat;do
    echo -n "!$j" >> .aux71
    echo -n "FALSE $j " >> .aux72
  done
  opcaodata=$(cat .aux71)
  rm .aux71
  opcaodata2=$(cat .aux72)
  rm .aux72
fi
}
read_pessoas(){
if [ -f $pessoas ];then
  pes=$(sort $pessoas)
  for j in $pes;do
    echo -n "!$j" >> .aux81
    echo -n "FALSE $j " >> .aux82
  done
  opcaopessoa=$(cat .aux81)
  rm .aux81
  opcaopessoa2=$(cat .aux82)
  rm .aux82
fi
}

#########  inserir novo objeto  #####################################################################################################
insert_new() {
  if [ "$cancel" != "1" ];then
    if [ "$category" = "" ];then
      categoria=$( \
        yad --form \
          --center \
          --title="CADASTRO DE IMAGENS" \
          --width=400 \
          --height=100 \
          --image="accessories-text-editor" \
          --text="Objeto: $adobjeto" \
          --field="Corpo Principal":CBE "$opcaocategoria" \
        )
      category=$(echo "$categoria" | cut -d"|" -f 1)
      if [ -z $category ];then cancel=1;fi
    fi
  fi
######### inserir nova categoria#####################################
  if [ "$cancel" != "1" ];then
    if [ -f $categoriax ];then
      cat $categoriax | egrep "$category" > /dev/null;
      if [ $? -eq 1 ]; then
        echo $category >> $categoriax
      fi
    else
      echo $category >> $categoriax
    fi
    subcategoriax="$atual/$category/.subcategoriax"              #arquivo com os dados das subcategorias
    if [ -f $subcategoriax ];then
      subcate=$(sort $subcategoriax)
      for j in $subcate;do
        echo -n "$j!" >> .aux0
      done
      VAR2=$(cat .aux0)
      opcaosubcategoria=$(echo ${VAR2%!*})
      rm .aux0
    fi
    subcategoria=$( \
      yad --form \
        --center \
        --title="CADASTRO DE IMAGENS" \
        --width=400 \
        --height=100 \
        --image="accessories-text-editor" \
        --text="Objeto: $adobjeto" \
        --field="Corpo Secundario":CBE "$opcaosubcategoria" \
      )
    subcategory=$(echo "$subcategoria" | cut -d"|" -f 1)
######## inserir nova subcategoria ####################################
    if [ -z $subcategory ];then cancel=1;fi
  fi
  if [ "$cancel" != "1" ];then
    if [ -f $subcategoriax ];then
      cat $subcategoriax | egrep "$subcategory" > /dev/null;
      if [ $? -eq 1 ]; then
        echo $subcategory >> $subcategoriax
      fi
    else
      echo $subcategory >> $subcategoriax
    fi
    echo -e "$adobjeto\t$category\t$subcategory" >> $objects
    if [ "$mainbody" != "" ];then
      echo -e "$mainbody\t$category\t$subcategory" >> $objects
    fi
    if [ ! -d "$atual/$category/$subcategory" ]; then
      mkdir $atual/$category/$subcategory
    fi
    if [ ! -d "$atual/$category/$subcategory/$adobjeto" ]; then
      mkdir $atual/$category/$subcategory/$adobjeto
    fi
  fi
}

########## adicionar imagens ao banco de dados ######################################################################################
adiciona() {
  informacoes="Observadores:\nCondicoes da Noite:"
  cancel=0
  until [ $cancel = "1" ]; do
    res=1
    mainbody=""
    novo1="sim"
    novo2="sim"
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    read_data;
    VAR_FORM=$( \
      yad --form \
        --center \
        --title="INSERCAO DE IMAGENS NO BANCO DE DADOS" \
        --width=400 \
        --height=700 \
        --image="accessories-text-editor" \
        --text="Insira o nome do objeto por * ou ** ou ***" \
        --field="*Objeto":CBE "$adobjeto$opcaoobjeto" \
        --field="**Estrela de ocultacao Est_Objeto_data":CHK FALSE \
        --field="**Objeto Ocultante":CBE "$opcaoobjeto3" \
        --field="**Data da Ocultacao":DT "" \
        --field="***Varios alvos no mesmo campo":CHK FALSE \
        --field="***Escolha o corpo principal dos objetos":CBE "$addmainbody!$opcaocategoria" \
        --field="Telescopio":CBE "$adtelescopio$opcaotelescopio" \
        --field="CCD":CBE "$adcamera$opcaocamera" \
        --field="Filtro":CBE "$adfiltro$opcaofiltro" \
        --field="Binagem":CBE "$adbinagem$opcaobinagem" \
        --field="Data":DT "$addata" --date-format="%Y_%b_%d" \
        --field="Imagens":MFL "" \
        --field="Renomear as imagens no formato \"Objeto_numero\"":CHK FALSE \
        --field="Relatorio":TXT "$informacoes"  \
        --field="README":BTN "yad --text-info --width=750 --height=400 --filename=$addreadme" \
      )
    adobjeto=$(echo "$VAR_FORM" | cut -d"|" -f 1)      ######
    estocc=$(echo "$VAR_FORM" | cut -d"|" -f 2)        ######
    mainbody=$(echo "$VAR_FORM" | cut -d"|" -f 3)      ######
    datocc=$(echo "$VAR_FORM" | cut -d"|" -f 4)        ######
    addmanyobj=$(echo "$VAR_FORM" | cut -d"|" -f 5)    ######
    addmainbody=$(echo "$VAR_FORM" | cut -d"|" -f 6)   ######
    adtelescopio=$(echo "$VAR_FORM" | cut -d"|" -f 7)  ######
    adcamera=$(echo "$VAR_FORM" | cut -d"|" -f 8)      ######
    adfiltro=$(echo "$VAR_FORM" | cut -d"|" -f 9)      ######
    adbinagem=$(echo "$VAR_FORM" | cut -d"|" -f 10)    ######
    addata=$(echo "$VAR_FORM" | cut -d"|" -f 11)       ######
    imagens=$(echo "$VAR_FORM" | cut -d"|" -f 12)      ######
    renomear=$(echo "$VAR_FORM" | cut -d"|" -f 13)     ######
    informacoes=$(echo "$VAR_FORM" | cut -d"|" -f 14)  ######
    if [ "$adobjeto" = "" -a "$estocc" = "FALSE" -a "$addmanyobj" = "FALSE" ];then cancel=1;fi
    if [ -z "$adtelescopio" ];then cancel=1;fi
    if [ -z "$adcamera" ];then cancel=1;fi
    if [ "$adobjeto" != "BIAS" -a "$adfiltro" = "" ];then cancel=1;fi
    if [ -z "$adbinagem" ];then cancel=1;fi
    if [ -z "$addata" ];then cancel=1;fi
    if [ -z "$imagens" ];then cancel=1;fi
    if [ "$cancel" != "1" ];then
      if [ "$addmanyobj" = "TRUE" ];then
        category="$addmainbody"
        subcategory="MIXED"
        if [ ! -d "$atual/$addmainbody" ]; then
          mkdir $atual/$addmainbody
        fi
        if [ ! -d "$atual/$addmainbody/MIXED" ]; then
          mkdir $atual/$addmainbody/MIXED
        fi
        if [ -f $objects ];then
          sort $objects > .auxadd
          while read COL_1 COL_2 COL_3;do
            if [ "$COL_2" = "$addmainbody" ];then
              echo -n "!$COL_1" >> .auxadd2
            fi
          done < .auxadd
          rm .auxadd
          opcaolista=$(cat .auxadd2)
          rm .auxadd2
        fi
        manyobjects=$(yad --form \
          --center \
          --width=400 \
          --height=200 \
          --text="Escolha os objetos que estao no campo\nA ordem nao importa" \
          --field="Objeto 1":CBE "$opcaolista" \
          --field="Objeto 2":CBE "$opcaolista" \
          --field="Objeto 3":CBE "$opcaolista" \
          --field="Objeto 4":CBE "$opcaolista" \
          --field="Objeto 5":CBE "$opcaolista" \
          --field="Objeto 6":CBE "$opcaolista" \
          --field="Objeto 7":CBE "$opcaolista" \
          --field="Objeto 8":CBE "$opcaolista" \
          --field="Objeto 9":CBE "$opcaolista" \
          --field="Objeto 10":CBE "$opcaolista" \
        )
        i=0
        j=0
        adobjetoa=""
        while [ $i -lt "10" ];do
          saida=$(echo "$manyobjects" | cut -d"|" -f $[$i+1])
          k=0
          varadd="sim"
          while [ $k -lt "${#mix[@]}" ];do
            if [ "${mix[$k]}" = "$saida" ];then
              varadd="nao"
              k="${#mix[@]}"
            fi
            k=$[$k+1]
          done
          if [ "$varadd" = "sim" -a "$saida" != "" ];then
            mix[$j]=$saida
            adobjetoa+="${mix[$j]}+"
            j=$[$j+1]
            novoob="sim"
            while read COL_1 COL_2 COL_3; do
              if [ "$COL_1" = "$saida" -a "$novoob" = "sim" ]; then
                novoob="nao"
              fi
            done < $objects
            if [ "$novoob" = "sim" ]; then
              adobjeto="$saida"; insert_new;
            fi
          fi
          i=$[$i+1]
        done
        adobjeto=$(echo ${adobjetoa%+*})
        existe="0"
        for m in $atual/$addmainbody/MIXED/*;do
          if [ "$(wc -l < $m/list_objects)" = "${#mix[*]}" -a "$existe" = "0" ];then
            listexiste=$(cat "$m/list_objects")
            n=0
            countlst="0"
            while [ $n -lt "${#mix[@]}" ];do
              for z in $listexiste;do
                if [ "$z" = "${mix[$n]}" ];then
                  countlst=$[$countlst+1]
                fi
              done
              n=$[$n+1]
            done
          fi
          if [ "$countlst" = "${#mix[*]}" -a "$existe" = "0" ];then
            adobjeto=${m##*/}
            existe="1"
          fi
        done
        if [ ! -d "$atual/$addmainbody/MIXED/$adobjeto" ]; then
          mkdir $atual/$addmainbody/MIXED/$adobjeto
        fi
        category="$addmainbody"
        subcategory="MIXED"
        if [ "$existe" != "1" ];then
          k=0
          while [ $k -lt "${#mix[@]}" ];do
            echo "${mix[$k]}" >> $atual/$category/$subcategory/$adobjeto/list_objects
            k=$[$k+1]
          done
        fi
        unset mix
      elif [ "$estocc" = "TRUE" ];then
        if [ "$mainbody" = "" -o "$datocc" = "" ];then
          perg=$(yad --form \
            --center \
            --width=400 \
            --height=200 \
            --text="Voce selecionou estrela de ocultacao.
Porem uma das informacoes abaixo esta vazia
Deseja continuar assim mesmo?
Objeto = Est_(Objeto Ocultante)_(Data da ocultacao)" \
            --field="**Objeto Ocultante":CBE "$mainbody$opcaoobjeto3" \
            --field="**Data da Ocultacao":DT "$datocc" --date-format="%Y_%b_%d" \
          )
          mainbody=$(echo "$perg" | cut -d"|" -f 1)      ######
          datocc=$(echo "$perg" | cut -d"|" -f 2)        ######
          if [ -z $perg ];then cancel=1;fi
        fi
        adobjeto=$(echo "Est_$mainbody""_$datocc")
      fi
      if  [ "$cancel" != "1" -a "$addmanyobj" = "FALSE" ];then
        if [ "$adobjeto" = "BIAS" ];then
          category="CALIBRACAO"
          subcategory="BIAS"
        elif [ "$adobjeto" = "FLAT" ];then
          category="CALIBRACAO"
          subcategory="FLAT"
        else
          while read COL_1 COL_2 COL_3
          do
            if [ "$COL_1" = "$adobjeto" -o "$COL_1" = "$mainbody" ]; then
              category="$COL_2"
              subcategory="$COL_3"
            fi
            if [ "$COL_1" = "$mainbody" ]; then
              novo2="nao"
            fi
            if [ "$COL_1" = "$adobjeto" ];then
              novo1="nao"
            fi
          done < $objects
          if [ "$novo1" = "sim" -a "$novo2" = "sim" ]; then
            category=""; insert_new;
          elif [ "$novo1" = "nao" -a "$novo2" = "sim" -a "$mainbody" != "" ]; then
            echo -e "$mainbody\t$category\t$subcategory" >> $objects
          elif [ "$novo1" = "sim" -a "$novo2" = "nao" ];then
            echo -e "$adobjeto\t$category\t$subcategory" >> $objects
          fi
        fi
      fi
      if [ "$cancel" != "1" ];then
        if [ ! -d "$atual/$category" ]; then
          mkdir $atual/$category
        fi
        if [ ! -d "$atual/$category/$subcategory" ]; then
          mkdir $atual/$category/$subcategory
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto
        fi
        cat $telescopes | egrep "$adtelescopio" > /dev/null;
        if [ $? -eq 1 ]; then
          echo -e "$adtelescopio" >> $telescopes
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio
        fi
        cat $cameras | egrep "$adcamera" > /dev/null;
        if [ $? -eq 1 ]; then
          echo -e "$adcamera" >> $cameras
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera
        fi
        cat $filtros | egrep "$adfiltro" > /dev/null;
        if [ $? -eq 1 ]; then
          echo -e "$adfiltro" >> $filtros
        fi
        if [ "$adobjeto" != "BIAS" ];then
          if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfiltro" ]; then
            mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfiltro
          fi
        fi
        x=$(echo ${adbinagem:1:1})
        if [ "$x" != "x" ]; then
          n=$(echo ${#adbinagem})
          if [ "$n" = "1" ]; then
            adbinagem="$adbinagem""x$adbinagem"
          elif [ "$n" = "2" ]; then
            if [ "${adbinagem:0:1}" = "${adbinagem:1:1}" ];then
              adbinagem="${adbinagem:0:1}x${adbinagem:0:1}"
            else
              adbinagem="$adbinagem""x$adbinagem"
            fi
          fi
        fi
        cat $binagens | egrep "$adbinagem" > /dev/null;
        if [ $? -eq 1 ]; then
          echo -e "$adbinagem" >> $binagens
        fi
        adfilter="$adfiltro"
        if [ "$adobjeto" = "BIAS" ];then
          adfilter=""
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem
        fi
        cat $calendario | egrep "$addata" > /dev/null;
        if [ $? -eq 1 ]; then
          echo -e "$addata" >> $calendario
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata
        fi
        if [ ! -d "$atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata/RAW" ]; then
          mkdir $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata/RAW
        fi
        arquivos=$(echo "$imagens" | tr '!' ' ')
        for i in $arquivos;do
          echo $i >> .aux
        done
        count=$(wc -l < .aux)
        rm .aux
        x=0
        quero="0"
        if [ "$count" != "0" ]; then
          if [ "$(ls -A $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata/RAW)" ];then
            yad --image "dialog-question" --width=300 --title "Verificacao" --button=Nao:1 --button=Sim:0 \
--text "Ja existem imagens para Objeto=$adobjeto; Telescopio:$adtelescopio; CCD:$adcamera; Filtro:$adfiltro; Binagem:$adbinagem; \
Data:$addata no banco de dados. Continuar pode sobescrever imagens. Deseja continuar?"
            subs=$?
            if [ "$subs" = "0" ];then
              quero="1"
            fi
          else
            quero="1"
          fi
          if [ "$quero" = "1" ];then
            (for file in $arquivos;do
              x=$[$x + 1]
              nome=""
              if [ "$renomear" = "TRUE" ];then
                q=${#x}
                zeros="0000"
                nx=${zeros:$q}
                nome="$adobjeto""_$nx""$x.fits"
              fi
              cp $file $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata/RAW/$nome
              porcentagemtot=$(echo "$x * 100 / $count"|bc -l)
              porcentagem=$(printf "%5.2f" $porcentagemtot)
              echo "#Objeto: $adobjeto; Telescopio: $adtelescopio; CCD: $adcamera; Filtro: $adfilter; Binagem: $adbinagem; Data: $addata; $porcentagem%"
              echo "$porcentagemtot"
            done) | yad --progress \
            --title="Copiando..." \
            --text="Copiando Imagens..." \
            --progress-text="Objeto: $adobjeto; Telescopio: $adtelescopio; CCD: $adcamera; Filtro: $adfilter; Binagem: $adbinagem; Data: $addata; 0.00%"   \
            --percentage=0 --auto-close
          fi
        fi
        echo -e "$informacoes" > $atual/$category/$subcategory/$adobjeto/$adtelescopio/$adcamera/$adfilter/$adbinagem/$addata/Relatorio
        yad --image "dialog-question" --title "COPIA COMPLETA" --button=Sair:1 --button="Analise":3 --button="Reducao":2 --button="Imagens":0 \
--text "Deseja adicionar mais imagens, arquivos de analise ou arquivos de reducao ao banco de dados?"
        res=$?
        case $res in
          0)cancel=0;;
          1)cancel=1;;
          2)cancel=1; redtel="$adtelescopio"; redcam="$adcamera"; redfil="$adfilter"; redbin="$adbinagem"; reddata="$addata"; insert_reduction;;
          3)cancel=1; anatel="$adtelescopio"; anacam="$adcamera"; anafil="$adfilter"; anabina="$adbinagem"; anadata="$addata"; insert_analise;;
        esac
      fi
    fi
  done
}

######### copiar imagens para algum lugar ###########################################################################################
copia() {
  res=0
  diretorio=$HOME
  until [ $res = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    read_data;
    cancelcopia=0
    askcopia=$( \
      yad --form \
        --center \
        --title="INSERCAO DE IMAGENS NO BANCO DE DADOS" \
        --width=500 \
        --height=400 \
        --image="accessories-text-editor" \
        --field="Diretório onde salvar as imagens":DIR "$diretorio" \
        --field="Separar por Objeto":CHK TRUE \
        --field="Todos os Telescopios":CHK TRUE \
        --field="Separar por Telecopio":CHK TRUE \
        --field="Todos os CCDs":CHK TRUE \
        --field="Separar por CCD":CHK TRUE \
        --field="Todos os Filtros":CHK TRUE \
        --field="Separar por Filtro":CHK TRUE \
        --field="Todas as Binagens":CHK TRUE \
        --field="Separar por Binagem":CHK TRUE \
        --field="Todas as datas disponiveis":CHK TRUE \
        --field="Separar por Data":CHK TRUE \
        --field="Copiar Bias e Flats para cada objeto":CHK FALSE \
        --field="Copiar Arquivos de Reducao":CHK FALSE \
        --field="Copiar Arquivos de Analise":CHK FALSE \
        --field="README":BTN "yad --text-info --width=650 --height=400 --filename=$copiareadme" \
      )
    diretorio=$(echo "$askcopia" | cut -d"|" -f 1)   ########
    ressepobj=$(echo "$askcopia" | cut -d"|" -f 2)   ########
    restel=$(echo "$askcopia" | cut -d"|" -f 3)      ########
    resseptel=$(echo "$askcopia" | cut -d"|" -f 4)   ########
    rescam=$(echo "$askcopia" | cut -d"|" -f 5)      ########
    ressepcam=$(echo "$askcopia" | cut -d"|" -f 6)   ########
    resfil=$(echo "$askcopia" | cut -d"|" -f 7)      ########
    ressepfil=$(echo "$askcopia" | cut -d"|" -f 8)   ########
    resbin=$(echo "$askcopia" | cut -d"|" -f 9)      ########
    ressepbin=$(echo "$askcopia" | cut -d"|" -f 10)  ########
    resdata=$(echo "$askcopia" | cut -d"|" -f 11)    ########
    ressepdata=$(echo "$askcopia" | cut -d"|" -f 12) ########
    resbifl=$(echo "$askcopia" | cut -d"|" -f 13)    ########
    resred=$(echo "$askcopia" | cut -d"|" -f 14)     ########
    resana=$(echo "$askcopia" | cut -d"|" -f 15)     ########
    if [ -z "$askcopia" ];then cancelcopia=1; res="1";fi
    if [ "$cancelcopia" = "0" ];then
      askobjeto=$(yad --text="Selecione o(s) objeto(s)" --list --checklist --column selecionar \
$opcaoobjeto2 --column opcao --hide-header --separator=" " --width=400 --height=400 --center --print-column=2)
      objeto=$(echo "$askobjeto" | cut -d"|" -f 1)
    fi
    if [ -z "$objeto" ];then cancelcopia=1; res="1";fi
    if [ "$restel" = "TRUE" -a "$cancelcopia" = "0" ];then
      telescopio=$(cat $telescopes)
    elif [ "$restel" = "FALSE" -a "$cancelcopia" = "0" ];then
      asktelescopio=$(yad --text="Selecione o(s) telescopio(s)" --list --checklist --column selecionar \
$opcaotelescopio2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      telescopio=$(echo "$asktelescopio" | cut -d"|" -f 1)
    fi
    if [ -z "$telescopio" ];then cancelcopia=1; res="1";fi
    if [ "$rescam" = "TRUE" -a "$cancelcopia" = "0" ];then
      camera=$(cat $cameras)
    elif [ "$rescam" = "FALSE" -a "$cancelcopia" = "0" ];then
      askcamera=$(yad --text="Selecione o(s) CCD(s)" --list --checklist --column selecionar \
$opcaocamera2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      camera=$(echo "$askcamera" | cut -d"|" -f 1)
    fi
    if [ -z "$camera" ];then cancelcopia=1; res="1";fi
    if [ "$resfil" = "TRUE" -a "$cancelcopia" = "0" ];then
      filtro=$(cat $filtros)
    elif [ "$resfil" = "FALSE" -a "$cancelcopia" = "0" ];then
      askfiltro=$(yad --text="Selecione o(s) filtro(s)" --list --checklist --column selecionar \
$opcaofiltro2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      filtro=$(echo "$askfiltro" | cut -d"|" -f 1)
    fi
    if [ -z "$filtro" ];then cancelcopia=1; res="1";fi
    if [ "$resbin" = "TRUE" -a "$cancelcopia" = "0" ];then
      binagem=$(cat $binagens)
    elif [ "$resbin" = "FALSE" -a "$cancelcopia" = "0" ];then
      askbinagem=$(yad --text="Selecione a(s) binagem(ns)" --list --checklist --column selecionar \
$opcaobinagem2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      binagem=$(echo "$askbinagem" | cut -d"|" -f 1)
    fi
    if [ -z "$binagem" ];then cancelcopia=1; res="1";fi
    if [ "$resdata" = "TRUE" -a "$cancelcopia" = "0" ];then
      data=$(cat $calendario)
    elif [ "$resdata" = "FALSE" -a "$cancelcopia" = "0" ];then
      askdata=$(yad --text="Selecione a(s) data(s)" --list --checklist --column selecionar \
$opcaodata2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      data=$(echo "$askdata" | cut -d"|" -f 1)
    fi
    if [ -z "$data" ];then cancelcopia=1; res="1";fi
    o2=""
    t2=""
    c2=""
    f2=""
    b2=""
    d2=""
    if [ "$cancelcopia" = "0" ];then
      for o in $objeto;do
        while read COL_1 COL_2 COL_3
        do
          if [ "$COL_1" = "$o" ]; then
            category="$COL_2"
            subcategory="$COL_3"
          fi
        done < $objects
        if [ "$ressepobj" = "TRUE" ];then
          o2="$o/"
        fi
        for t in $telescopio;do
          if [ "$resseptel" = "TRUE" ];then
            t2="$t/"
          fi
          for c in $camera;do
            if [ "$ressepcam" = "TRUE" ];then
              c2="$c/"
            fi
            for f in $filtro;do
              if [ "$ressepfil" = "TRUE" ];then
                f2="$f/"
              fi
              for b in $binagem;do
                if [ "$ressepbin" = "TRUE" ];then
                  b2="$b/"
                fi
                for d in $data;do
                  if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d" ];then
                    if [ "$ressepdata" = "TRUE" ];then
                      d2="$d/"
                    fi
                    if [ ! -d $diretorio/$o2 ];then mkdir $diretorio/$o2; fi
                    if [ ! -d $diretorio/$o2$t2 ];then mkdir $diretorio/$o2$t2; fi
                    if [ ! -d $diretorio/$o2$t2$c2 ];then mkdir $diretorio/$o2$t2$c2; fi
                    if [ ! -d $diretorio/$o2$t2$c2$f2 ];then mkdir $diretorio/$o2$t2$c2$f2; fi
                    if [ ! -d $diretorio/$o2$t2$c2$f2$b2 ];then mkdir $diretorio/$o2$t2$c2$f2$b2; fi
                    if [ ! -d $diretorio/$o2$t2$c2$f2$b2$d2 ];then mkdir $diretorio/$o2$t2$c2$f2$b2$d; fi
                    if [ ! -d $diretorio/$o2$t2$c2$f2$b2$d2"RAW" ];then mkdir $diretorio/$o2$t2$c2$f2$b2$d2"RAW"; fi
                    ls $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/RAW > .auxb
                    imagens=$(cat .auxb)
                    rm .auxb
                    if [ -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/"REDUCTION" ];then
                      if [ "$resred" = "TRUE" ];then
                        if [ ! -d $diretorio/$o2$t2$c2$f2$b2$d2"REDUCTION" ];then mkdir $diretorio/$o2$t2$c2$f2$b2$d2"REDUCTION";fi
                        echo -e "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/* $diretorio/$o2$t2$c2$f2$b2$d2""REDUCTION" >> .auxt
                      fi
                    fi
                    if [ -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/"ANALISE" ];then
                      if [ "$resana" = "TRUE" ];then
                        if [ ! -d $diretorio/$o2$t2$c2$f2$b2$d/ANALISE ];then mkdir $diretorio/$o2$t2$c2$f2$b2$d/ANALISE; fi
                        echo -e "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/* $diretorio/$o2$t2$c2$f2$b2$d2""ANALISE" >> .auxt
                      fi
                    fi
                    for j in $imagens;do
                      echo -e "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/RAW/$j $diretorio/$o2$t2$c2$f2$b2$d2""RAW" >> .auxt
                    done
                    if [ "$resbifl" = "TRUE" ];then
                      if [ -d "$atual/CALIBRACAO/BIAS/BIAS/$t/$c/$b/$d" ];then
                        if [ ! -d "$diretorio/$o2$t2$c2$f2$b2$d2""BIAS" ];then mkdir "$diretorio/$o2$t2$c2$f2$b2$d2""BIAS"; fi
                        ls $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/RAW > .auxbi
                        list_bias=$(cat .auxbi)
                        rm .auxbi
                        for j in $list_bias;do
                          cat .auxt | egrep "'$diretorio/$o2$t2$c2$f2$b2$d2''BIAS/$j'" > /dev/null;
                          if [ $? -eq 1 ]; then
                            echo -e "$atual/CALIBRACAO/BIAS/BIAS/$t/$c/$b/$d/RAW/$j $diretorio/$o2$t2$c2$f2$b2$d2""BIAS" >> .auxt
                          fi
                        done
                      fi
                      if [ -d "$atual/CALIBRACAO/FLAT/FLAT/$t/$c/$f/$b/$d" ];then
                        if [ ! -d "$diretorio/$o2$t2$c2$f2$b2$d2""FLAT" ];then mkdir "$diretorio/$o2$t2$c2$f2$b2$d2""FLAT"; fi
                        ls $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/RAW > .auxfl
                        list_flat=$(cat .auxfl)
                        rm .auxfl
                        for j in $list_flat;do
                          cat .auxt | egrep "'$diretorio/$o2$t2$c2$f2$b2$d2''FLAT/$j'" > /dev/null;
                          if [ $? -eq 1 ]; then
                            echo -e "$atual/CALIBRACAO/FLAT/FLAT/$t/$c/$f/$b/$d/RAW/$j $diretorio/$o2$t2$c2$f2$b2$d2""FLAT" >> .auxt
                          fi
                        done
                      fi
                    fi
                  fi
                done
              done
            done
          done
        done
      done
      count=$(wc -l < .auxt)
      x=0
      if [ "$count" != "0" ]; then
        (while read inicio fim
        do
          cp $inicio $fim
          x=$[$x + 1]
          porcentagemtot=$(echo "$x * 100 / $count"|bc -l)
          porcentagem=$(printf "%5.2f" $porcentagemtot)
          echo "# Copiando imagens: $porcentagem %"
          echo "$porcentagemtot"
        done < .auxt) | yad --progress \
        --title="Copia" \
        --text="Copiando imagens para diretorio $diretorio..." \
        --progress-text="Copiando imagens: 0.00 %" \
        --percentage=0 --auto-close
        rm .auxt
      fi
      yad --image "dialog-question" --title "COPIA COMPLETA" --button=Nao:1 --button=Sim:0 \
--text "Deseja copiar mais imagens para o seu diretorio?"
      res=$?
    fi
  done
}

######### adiciona arquivos de analise ##############################################################################################
insert_analise(){
  cancel=0
  anacomment="Programa e versao utilizada na analise:\nComentarios gerais sobre o que foi feito"
  until [ $cancel = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    read_data;
    read_pessoas;
    ANA_FORM=$( \
      yad --form \
        --center \
        --title="INSERCAO DE ARQUIVOS DE ANALISE NO BANCO DE DADOS" \
        --width=500 \
        --height=600 \
        --image="accessories-text-editor" \
        --field="Todos os Objetos":CHK FALSE \
        --field="Telescopio":CB "$anatel$opcaotelescopio" \
        --field="CCD":CB "$anacam$opcaocamera" \
        --field="Todos os Filtros":CHK FALSE \
        --field="Filtro":CB "$anafilt$opcaofiltro" \
        --field="Todas as Binagens":CHK FALSE \
        --field="Binagem":CB "$anabina$opcaobinagem" \
        --field="Data da observacao":CB "$anadata$opcaodata" \
        --field="Autor da Analise":CBE "$anapes$opcaopessoa" \
        --field="Ano e mes da Analise" "Ano_Mes" \
        --field="Diretório onde estao os arquivos de analise":DIR "$directoryana" \
        --field="Comentarios sobre a analise":TXT "$anacomment" \
        --field="README":BTN "yad --text-info --width=600 --height=400 --filename=$anareadme" \
      )
    if [ -z "$ANA_FORM" ];then cancel="1"; fi
    anaobj=$(echo "$ANA_FORM" | cut -d"|" -f 1)
    anatel=$(echo "$ANA_FORM" | cut -d"|" -f 2)
    anacam=$(echo "$ANA_FORM" | cut -d"|" -f 3)
    anafil=$(echo "$ANA_FORM" | cut -d"|" -f 4)
    anafilt=$(echo "$ANA_FORM" | cut -d"|" -f 5)
    anabin=$(echo "$ANA_FORM" | cut -d"|" -f 6)
    anabina=$(echo "$ANA_FORM" | cut -d"|" -f 7)
    anadata=$(echo "$ANA_FORM" | cut -d"|" -f 8)
    anapes=$(echo "$ANA_FORM" | cut -d"|" -f 9)
    anaano=$(echo "$ANA_FORM" | cut -d"|" -f 10)
    directoryana=$(echo "$ANA_FORM" | cut -d"|" -f 11)
    anacomment=$(echo "$ANA_FORM" | cut -d"|" -f 12)
    objeto=$listobjeto
    if [ -z "$anatel" ];then cancel=1;fi
    if [ -z "$anacam" ];then cancel=1;fi
    if [ "$anafil" = "FALSE" -a -z "$anafilt" ];then cancel=1;fi
    if [ "$anabin" = "FALSE" -a -z "$anabina" ];then cancel=1;fi
    if [ -z "$anadata" ];then cancel=1;fi
    if [ -z "$anapes" ];then cancel=1;fi
    if [ -z "$directoryana" ];then cancel=1;fi
    if [ "$cancel" != "1" ];then
      if [ "$anaobj" = "FALSE" ];then
        askobjeto=$(yad --text="Selecione os objetos" --list --checklist --column selecionar \
$opcaoobjeto2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
        objeto=$(echo "$askobjeto" | cut -d"|" -f 1)
      fi
      anapasta="$anapes""_$anaano"
      cat $pessoas | egrep "$anapes" > /dev/null;
      if [ $? -eq 1 ]; then
        echo -e "$anapes" >> $pessoas;
      fi
      for o in $objeto;do
        while read COL_1 COL_2 COL_3
        do
          if [ "$COL_1" = "$o" ]; then
            category="$COL_2"
            subcategory="$COL_3"
          fi
        done < $objects
        t=$anatel
        if [ -d "$atual/$category/$subcategory/$o/$t" ];then
          c=$anacam
          if [ -d "$atual/$category/$subcategory/$o/$t/$c" ];then
            filt=$(cat $filtros)
            if [ "$anafil" = "FALSE" ];then filt="$anafilt";fi
            for f in $filt;do
              if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f" ];then
                bins=$(cat $binagens)
                if [ "$anabin" = "FALSE" ];then bins=$anabina;fi
                for b in $bins;do
                  if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b" ];then
                    d=$anadata
                    if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d" ];then
                      if [ ! -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE ];then
                        mkdir $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE
                      fi
                      if [ ! -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta ];then
                        mkdir $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta
                        cp -R $directoryana/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta/
                      else
                        if [ "$(ls -A $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta)" ];then
                          yad --image "dialog-question" --width=300 --title "Verificacao" --button=Nao:1 --button=Sim:0 \
--text "Ja existem arquivos de analise para Objeto=$o; Telescopio:$t; CCD:$c; Filtro:$f; Binagem:$b; \
Data:$d feita por $anapasta. Deseja continuar e substituir os dados?"
                          subs=$?
                          if [ "$subs" = "0" ];then
                            rm $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta/*
                            cp -R $directoryana/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta
                          fi
                        else
                          cp -R $directoryana/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta/
                        fi
                      fi
                      echo -e "$anacomment" > $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/ANALISE/$anapasta/Comentarios
                    fi
                  fi
                done
              fi
            done
          fi
        fi
      done
      yad --image "dialog-question" --title "COPIA COMPLETA" --button=Nao:1 --button=Sim:0 \
--text "Deseja copiar mais arquivos de analise para o seu banco de dados?"
      cancel=$?
    fi
  done
}

######### adiciona arquivos de reducao ##############################################################################################
insert_reduction(){
  cancel=0
  redcomment="Programa e versao utilizada na reducao:\nComentarios gerais sobre o que foi feito"
  until [ $cancel = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    read_data;
    read_pessoas;
    RED_FORM=$( \
      yad --form \
        --center \
        --title="INSERCAO DE ARQUIVOS DE REDUCAO NO BANCO DE DADOS" \
        --width=500 \
        --height=600 \
        --image="accessories-text-editor" \
        --field="Todos os Objeto":CHK FALSE \
        --field="Telescopio":CB "$redtel$opcaotelescopio" \
        --field="CCD":CB "$redcam$opcaocamera" \
        --field="Todos os Filtros":CHK FALSE \
        --field="Filtro":CB "$redfilt$opcaofiltro" \
        --field="Todas as Binagens":CHK FALSE \
        --field="Binagem":CB "$redbin$opcaobinagem" \
        --field="Data da observacao":CB "$reddata$opcaodata" \
        --field="Autor da Reducao":CBE "$redpes$opcaopessoa" \
        --field="Ano e mes da Reducao" "Ano_Mes" \
        --field="Diretório onde estao os arquivos de reducao":DIR "$directoryred" \
        --field="Comentarios sobre a analise":TXT "$redcomment" \
        --field="README":BTN "yad --text-info --width=620 --height=400 --filename=$redreadme" \
      )
    if [ -z "$RED_FORM" ];then cancel="1"; fi
    redobj=$(echo "$RED_FORM" | cut -d"|" -f 1)
    redtel=$(echo "$RED_FORM" | cut -d"|" -f 2)
    redcam=$(echo "$RED_FORM" | cut -d"|" -f 3)
    redfil=$(echo "$RED_FORM" | cut -d"|" -f 4)
    redfilt=$(echo "$RED_FORM" | cut -d"|" -f 5)
    redbin=$(echo "$RED_FORM" | cut -d"|" -f 6)
    redbina=$(echo "$RED_FORM" | cut -d"|" -f 7)
    reddata=$(echo "$RED_FORM" | cut -d"|" -f 8)
    redpes=$(echo "$RED_FORM" | cut -d"|" -f 9)
    redano=$(echo "$RED_FORM" | cut -d"|" -f 10)
    directoryred=$(echo "$RED_FORM" | cut -d"|" -f 11)
    redcomment=$(echo "$ANA_FORM" | cut -d"|" -f 12)
    objeto=$listobjeto
    if [ -z "$redtel" ];then cancel=1;fi
    if [ -z "$redcam" ];then cancel=1;fi
    if [ "$redfil" = "FALSE" -a -z "$redfilt" ];then cancel=1;fi
    if [ "$redbin" = "FALSE" -a -z "$redbina" ];then cancel=1;fi
    if [ -z "$reddata" ];then cancel=1;fi
    if [ -z "$redpes" ];then cancel=1;fi
    if [ -z "$directoryred" ];then cancel=1;fi
    if [ "$cancel" != "1" ];then
      if [ "$redobj" = "FALSE" ];then
        askobjeto=$(yad --text="Selecione os objetos" --list --checklist --column selecionar \
$opcaoobjeto2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
        objeto=$(echo "$askobjeto" | cut -d"|" -f 1)
      fi
      redpasta="$redpes""_$redano"
      cat $pessoas | egrep "$redpes" > /dev/null;
      if [ $? -eq 1 ]; then
        echo -e "$redpes" >> $pessoas;
      fi
      for o in $objeto;do
        while read COL_1 COL_2 COL_3
        do
          if [ "$COL_1" = "$o" ]; then
            category="$COL_2"
            subcategory="$COL_3"
          fi
        done < $objects
        t=$redtel
        if [ -d "$atual/$category/$subcategory/$o/$t" ];then
          c=$redcam
          if [ -d "$atual/$category/$subcategory/$o/$t/$c" ];then
            filt=$(cat $filtros)
            if [ "$redfil" = "FALSE" ];then filt=$redfilt;fi
            for f in $filt;do
              if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f" ];then
                bins=$(cat $binagens)
                if [ "$redbin" = "FALSE" ];then bins=$redbina;fi
                for b in $bins;do
                  if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b" ];then
                    d=$reddata
                    if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$d" ];then
                      if [ ! -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION ];then
                        mkdir $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION
                      fi
                      if [ ! -d $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta ];then
                        mkdir $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta
                        cp -R $directoryana/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta
                        else
                        if [ "$(ls -A $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta)" ];then
                          yad --image "dialog-question" --width=300 --title "Verificacao" --button=Nao:1 --button=Sim:0 \
--text "Ja existem arquivos de reducao para Objeto=$o; Telescopio:$t; CCD:$c; Filtro:$f; Binagem:$b; \
Data:$d feita por /$redpasta. Deseja continuar e substituir os dados?"
                          subs=$?
                          if [ "$subs" = "0" ];then
                            rm $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta/*
                            cp -R $directoryred/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta
                          fi
                        else
                          cp -R $directoryred/* $atual/$category/$subcategory/$o/$t/$c/$f/$b/$d/REDUCTION/$redpasta
                        fi
                      fi
                    fi
                  fi
                done
              fi
            done
          fi
        fi
      done
      yad --image "dialog-question" --title "COPIA COMPLETA" --button=Nao:1 --button=Sim:0 \
--text "Deseja copiar mais arquivos de reducao para o seu banco de dados?"
      cancel=$?
    fi
  done
}

######### gera uma estatistica do banco de dados ####################################################################################
estatistica_database(){
arq_estat="$atual/estatistica"
#(
  countobj=$(wc -l < $objects)
  counttel=$(wc -l < $telescopes)
  countcam=$(wc -l < $cameras)
  countbin=$(wc -l < $binagens)
  countfil=$(wc -l < $filtros)
  countdat=$(wc -l < $calendario)
  sort -k 2 -k 3 -k 1 $objects > .aux00
  t=$(sort $telescopes)
  c=$(sort $cameras)
  f=$(sort $filtros)
  b=$(sort $binagens)
  d=$(sort -k 1.1,1.4 -k 1.6,1.8M -k 1.10,1.11 $calendario)
  anos=""
  for obser in $d;do
    echo "$anos" | egrep "${obser%%_*}" >/dev/null
    if [ $? -ne 0 ];then anos="$anos ${obser%%_*}";fi
  done
  echo "Estatistica do Banco de Dados" > $arq_estat
  while read COL_1 COL_2 COL_3
  do
    obje="$COL_1"
    categ="$COL_2"
    subcateg="$COL_3"
    find "$atual/$categ/$subcateg/$obje" -name *.fits >> .aux01
  done < .aux00
  imagens=$(wc -l < .aux01)
  rm .aux01
  find "$atual/CALIBRACAO/BIAS/BIAS/" -name *.fits > .aux02
  bias=$(wc -l < .aux02)
  rm .aux02
  find "$atual/CALIBRACAO/FLAT/FLAT/" -name *.fits > .aux03
  flat=$(wc -l < .aux03)
  rm .aux03
  echo -e "\n\nNumero total de imagens no banco de dados: $imagens" >> $arq_estat
  echo -e "Numero de imagens Bias no banco de dados: $bias" >> $arq_estat
  echo -e "Numero de imagens Flat no banco de dados: $flat" >> $arq_estat
  declare -A tldtob tldtcm
  i=0; for k in $t;do tl[$i]="$k"; j=0; for l in $d;do tldtob["$i:$j"]=""; tldtcm["$i:$j"]=""; j=$[j+1]; done; i=$[i+1];done
  i=0; for k in $c;do cm[$i]="$k"; i=$[i+1];done
  i=0; for k in $f;do flt[$i]="$k"; i=$[i+1];done
  i=0; for k in $b;do bng[$i]="$k"; i=$[i+1];done
  i=0; for k in $d;do dt[$i]="$k"; i=$[i+1];done
  while read COL_1 COL_2 COL_3
  do
    obje="$COL_1"
    categ="$COL_2"
    subcateg="$COL_3"
    imobj=0
    i=0; for k in $t;do imobtel[$i]=0; i=$[i+1];done
    i=0; for k in $c;do imobcam[$i]=0; i=$[i+1];done
    i=0; for k in $f;do imobfil[$i]=0; i=$[i+1];done
    i=0; for k in $b;do imobbin[$i]=0; i=$[i+1];done
    i=0; for k in $d;do imobdat[$i]=0; i=$[i+1];done
    obdat=""
    j=0
    declare -A tlandt
    for te in $t;do
      k=0
      for year in $anos; do
        tlandt["$j:$year"]=""
      done
      for ca in $c;do
        l=0
        for fi in $f;do
          m=0
          for bi in $b;do
            n=0
            for da in $d;do
              if [ -d "$atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da" ];then
                ls $atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da/RAW > .auxb
                numim=$(wc -l < .auxb)
                rm .auxb
                imobj=$(echo "$imobj + $numim" | bc -l)
                imobtel[$j]=$(echo "${imobtel[$j]} + $numim" | bc -l)
                imobcam[$k]=$(echo "${imobcam[$k]} + $numim" | bc -l)
                imobfil[$l]=$(echo "${imobfil[$l]} + $numim" | bc -l)
                imobbin[$m]=$(echo "${imobbin[$m]} + $numim" | bc -l)
                imobdat[$n]=$(echo "${imobdat[$n]} + $numim" | bc -l)
                echo "${tldtob["$j:$n"]}" | egrep "$obje;" >/dev/null
                if [ $? -ne 0 ];then tldtob["$j:$n"]="${tldtob["$j:$n"]} $obje;";fi
                echo "${tldtcm["$j:$n"]}" | egrep "$ca" >/dev/null
                if [ $? -ne 0 ];then tldtcm["$j:$n"]="${tldtcm["$j:$n"]}$ca";fi
                cat ".aux$te${da%%_*}" | egrep "${da#*_}" >/dev/null
                if [ $? -ne 0 ];then echo "${da#*_}" >> ".aux$te${da%%_*}";fi
                if [ -d $atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da/"REDUCTION" ];then
                  reduclist=$(ls $atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da/"REDUCTION")
                fi
                if [ -d $atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da/"ANALISE" ];then
                  analist=$(ls $atual/$categ/$subcateg/$obje/$te/$ca/$fi/$bi/$da/"ANALISE")
                fi
                if [ -d $atual/CALIBRACAO/BIAS/BIAS/$te/$ca/$bi/$da/RAW ];then
                  biasexi=1
                fi
                if [ -d $atual/CALIBRACAO/FLAT/FLAT/$te/$ca/$fi/$bi/$da/RAW ];then
                  flatexi=1
                fi
              fi
              n=$[$n+1]
            done
            m=$[$m+1]
          done
          l=$[$l+1]
        done
        k=$[$k+1]
      done
      j=$[$j+1]
    done
    if [ "$imobj" != "0" ];then
      echo -e "\n\nObjeto:$obje\n\tCorpo Principal: $categ; Corpo Secundario: $subcateg" >> $arq_estat
      echo -e "\tQuantidade total de imagens: $imobj" >> $arq_estat
      i=0
      for k in $t;do
        if [ "${imobtel[$i]}" != "0" ];then
          echo -e "\tImagens no telescopio ${tl[$i]} = ${imobtel[$i]}" >> $arq_estat
          for year in $anos;do
            if [ -e ".aux$k$year" ];then
              tlandt["$i:$year"]=$(sort -Mk 1 -k2 --field-separator="_" ".aux$k$year")
              rm ".aux$k$year"
              tlandt["$i:$year"]=${tlandt["$i:$year"]//
/; }
              echo -e "\t\t$year(${tlandt["$i:$year"]})"  >> $arq_estat
            fi
          done
        fi;
        i=$[i+1]
      done
    fi
  done < .aux00
  rm .aux00
  echo -e "\n\nObjetos observados por noite por Telescopio:" >> $arq_estat
  i=0;
  for k in $t;do
    echo -e "\n${tl[$i]}: " >> $arq_estat
    j=0;
    for l in $d;do
      if [ "${tldtob["$i:$j"]}" != "" ];then echo -e "\t${dt[$j]} = (${tldtcm["$i:$j"]}) ${tldtob["$i:$j"]}" >> $arq_estat;fi;
      j=$[j+1];
    done
    i=$[i+1];
  done
#) | yad --progress \
#--title="Estatistica..." \
#--text="Gerando estatistica do banco de dados..." \
#--progress-text="Contando imagens: 0.00%"   \
#--pulsate --auto-close
yad --text-info --width=600 --height=400 --filename="$arq_estat"
}

#
#Pergunta ao usuario sobre a tarefa a ser feita #####################################################################################
#
fazer=1
until [ "$fazer" = "0" ]; do
  inicial=$(yad --text="Escolha a opçao que deseja executar" --list --radiolist --column selecionar \
  --width=440 --height=250 --center --title="DATABASE LNA" \
FALSE 1 "Adicionar imagens no banco de dados" FALSE 2 "Adicionar arquivos de analise no banco de dados" \
FALSE 3 "Adicionar arquivos de reducao no banco de dados" FALSE 4 "Copiar imagens para algum diretorio" \
FALSE 5 "Estatistica do banco de dados" \
--column opcao --column tarefa --hide-column=2 --hide-header --print-column=2)
  tarefa=$(echo "$inicial" | cut -d"|" -f 1)
  case $tarefa in
    1) adiciona;;
    2) insert_analise;;
    3) insert_reduction;;
    4) copia;;
    5) estatistica_database;;
  esac
 if [ -z "$inicial" ];then
  fazer="0"
 fi
done
