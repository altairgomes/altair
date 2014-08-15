atual=$(pwd)
objects="$atual/.objects"                    #arquivo com os dados dos objeto
calibracao="$atual/.calibra"                 #arquivo com os dados das calibrações
categoriax="$atual/.categoriax"              #arquivo com os dados das categorias
telescopes="$atual/.telescopes"              #arquivo com os dados dos telescopio
cameras="$atual/.cameras"                    #arquivo com os dados das cameras
binagens="$atual/.binagens"                  #arquivo com os dados das binagens
filtros="$atual/.filtros"                    #arquivo com os dados dos filtros

read_categoria(){
if [ -e $categoriax ];then
  cate=$(cat $categoriax)
  for j in $cate;do
    echo -n "$j!" >> .aux1
    echo -n "FALSE $j " >> .aux12
  done
  opcaocategoria=$(cat .aux1)
  rm .aux1
  opcaocategoria2=$(cat .aux12)
  rm .aux12
fi
}
read_object(){
if [ -e $objects ];then
  while read COL_1 COL_2 COL_3;do
    echo -n "!$COL_1" >> .aux2
    echo -n "FALSE $COL_1 " >> .aux22
  done < $objects
  echo -n "!BIAS!FLAT" >> .aux2
  opcaoobjeto=$(cat .aux2)
  rm .aux2
  opcaoobjeto2=$(cat .aux22)
  rm .aux22
fi
}
read_telescope(){
if [ -e $telescopes ];then
  teles=$(cat $telescopes)
  for j in $teles;do
    echo -n "!$j" >> .aux3
    echo -n "FALSE $j " >> .aux32
  done
  opcaotelescopio=$(cat .aux3)
  rm .aux3
  opcaotelescopio2=$(cat .aux32)
  rm .aux32
fi
}
read_camera(){
if [ -e $cameras ];then
  came=$(cat $cameras)
  for j in $came;do
    echo -n "!$j" >> .aux4
    echo -n "FALSE $j " >> .aux42
  done
  opcaocamera=$(cat .aux4)
  rm .aux4
  opcaocamera2=$(cat .aux42)
  rm .aux42
fi
}
read_binagem(){
if [ -e $binagens ];then
  bina=$(cat $binagens)
  for j in $bina;do
    echo -n "!$j" >> .aux5
    echo -n "FALSE $j " >> .aux52
  done
  opcaobinagem=$(cat .aux5)
  rm .aux5
  opcaobinagem2=$(cat .aux52)
  rm .aux52
fi
}
read_filtro(){
if [ -e $filtros ];then
  filt=$(cat $filtros)
  for j in $filt;do
    echo -n "!$j" >> .aux6
    echo -n "FALSE $j " >> .aux62
  done
  opcaofiltro=$(cat .aux6)
  rm .aux6
  opcaofiltro2=$(cat .aux62)
  rm .aux62
fi
}

#########  inserir novo objeto  ########################################################################################################
insert_new() {
  new_object=$objeto
  categoria=$( \
    yad --form \
      --center \
      --title="CADASTRO DE IMAGENS" \
      --width=400 \
      --height=100 \
      --image="accessories-text-editor" \
      --field="Categoria":CBE "$opcaocategoria" \
    )
  category=$(echo "$categoria" | cut -d"|" -f 1)
######### inserir nova categoria#####################################
  if [ -e $categoriax ];then
    cat $categoriax | egrep "$category" > /dev/null;
    if [ $? -eq 1 ]; then
      new_category=$category
      mkdir $atual/$new_category
      echo $new_category >> $categoriax
    fi
  else
    new_category=$category
    mkdir $atual/$new_category
    echo $new_category >> $categoriax
  fi
  subcategoriax="$atual/$category/.categoriax"              #arquivo com os dados das subcategorias
  if [ -e $subcategoriax ];then
    subcate=$(cat $subcategoriax)
    for j in $subcate;do
      echo -n "$j!" >> .aux0
    done
    opcaosubcategoria=$(cat .aux0)
    rm .aux0
  fi
  subcategoria=$( \
    yad --form \
      --center \
      --title="CADASTRO DE IMAGENS" \
      --width=400 \
      --height=100 \
      --image="accessories-text-editor" \
      --field="Subcategoria":CBE "$opcaocategoria" \
    )
  subcategory=$(echo "$subcategoria" | cut -d"|" -f 1)
######## inserir nova subcategoria ####################################
  if [ -e $subcategoriax ];then
    cat $subcategoriax | egrep "$subcategory" > /dev/null;
    if [ $? -eq 1 ]; then
      new_subcategory=$subcategory
      mkdir $atual/$category/$new_subcategory
      echo $new_subcategory >> $subcategoriax
    fi
  else
    new_subcategory=$subcategory
    mkdir $atual/$category/$new_subcategory
    echo $new_subcategory >> $subcategoriax
  fi
  mkdir $atual/$category/$subcategory/$new_object
  echo -e "$new_object\t$category\t$subcategory" >> $objects
}

#########  inserir novo telescopio  ############################################################################################
insert_telescope() {
  new_telescope=$telescopio
  echo -e "$new_telescope" >> $telescopes
}

#########  inserir nova camera  ########################################################################################################
insert_cam() {
  new_camera=$camera
  echo -e "$new_camera" >> $cameras
}

#########  inserir nova binagem  ################################################################################################
insert_bin() {
  new_bin=$binagem
  echo -e "$new_bin" >> $binagens
}

#########  inserir novo filtro  ###################################################################################################
insert_fil() {
  new_filter=$filtro
  echo -e "$new_filter" >> $filtros
}

########## adicionar imagens ao banco de dados ######################################################################################
adiciona() {
  res=0
  informacoes="Observadores:\nCondicoes da Noite:"
  cancel=0
  until [ $res = "1" ]; do
    res=1
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    VAR_FORM=$( \
      yad --form \
        --center \
        --title="INSERCAO DE IMAGENS NO BANCO DE DADOS" \
        --width=400 \
        --height=400 \
        --image="accessories-text-editor" \
        --field="Objeto":CBE "$objeto$opcaoobjeto" \
        --field="Telescopio":CBE "$telescopio$opcaotelescopio" \
        --field="Camera":CBE "$camera$opcaocamera" \
        --field="Filtro":CBE "$filtro$opcaofiltro" \
        --field="Binagem":CBE "$binagem$opcaobinagem" \
        --field="Data":DT "$data" --date-format=%F \
        --field="Imagens":MFL "$imagens" \
        --field="Observadores":TXT "$informacoes"  \
      )
    objeto=$(echo "$VAR_FORM" | cut -d"|" -f 1)
    telescopio=$(echo "$VAR_FORM" | cut -d"|" -f 2)
    camera=$(echo "$VAR_FORM" | cut -d"|" -f 3)
    filtro=$(echo "$VAR_FORM" | cut -d"|" -f 4)
    binagem=$(echo "$VAR_FORM" | cut -d"|" -f 5)
    data=$(echo "$VAR_FORM" | cut -d"|" -f 6)
    imagens=$(echo "$VAR_FORM" | cut -d"|" -f 7)
    informacoes=$(echo "$VAR_FORM" | cut -d"|" -f 8)
    if [ -z $objeto ];then cancel=1;fi
    if [ -z $telescopio ];then cancel=1;fi
    if [ -z $camera ];then cancel=1;fi
    if [ -z $binagem ];then cancel=1;fi
    if [ -z $data ];then cancel=1;fi
    if [ "$cancel" != "1" ];then
      if [ "$objeto" = "BIAS" ];then
        category="CALIBRACAO"
        subcategory="BIAS"
      elif [ "$objeto" = "FLAT" ];then
        category="CALIBRACAO"
        subcategory="FLAT"
      else
        echo "$opcaoobjeto" | egrep "$objeto" > /dev/null;
        if [ $? -eq 1 ]; then
          insert_new;
        fi
        while read COL_1 COL_2 COL_3
        do
          if [ "$COL_1" = "$objeto" ]; then
            category="$COL_2"
            subcategory="$COL_3"
          fi
        done < $objects
      fi
      if [ ! -d "$atual/$category" ]; then
        mkdir $atual/$category
      fi
      if [ ! -d "$atual/$category/$subcategory" ]; then
        mkdir $atual/$category/$subcategory
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto" ]; then
        mkdir $atual/$category/$subcategory/$objeto
      fi
      cat $telescopes | egrep "$telescopio" > /dev/null;
      if [ $? -eq 1 ]; then
        insert_telescope;
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio
      fi
      cat $cameras | egrep "$camera" > /dev/null;
      if [ $? -eq 1 ]; then
        insert_cam;
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera
      fi
      cat $filtros | egrep "$filtro" > /dev/null;
      if [ $? -eq 1 ]; then
        insert_fil;
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro
      fi
      cat $binagens | egrep "$binagem" > /dev/null;
      if [ $? -eq 1 ]; then
        insert_bin;
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data
      fi
      if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data/RAW" ]; then
        mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data/RAW
      fi
      arquivos=$(echo "$imagens" | tr '!' ' ')
      for i in $arquivos;do
        echo $i >> .aux
      done
      count=$(wc -l < .aux)
      rm .aux
      x=0
      if [ "$count" != "0" ]; then
        (for file in $arquivos;do
          cp $file $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data/RAW
          x=$[$x + 1]
          porcentagemtot=$(echo "$x * 100 / $count"|bc -l)
          porcentagem=$(printf "%.2f" $porcentagemtot)
          echo "#Objeto: $objeto; Telescopio: $telescopio; Camera: $camera; Filtro: $filtro; Binagem: $binagem; Data: $data; $porcentagem%"
          echo "$porcentagemtot"
        done) | yad --progress \
        --title="Copiando..." \
        --text="Copiando Imagens..." \
        --progress-text="Objeto: $objeto; Telescopio: $telescopio; Camera: $camera; Filtro: $filtro; Binagem: $binagem; Data: $data; 0.00%"   \
        --percentage=0
      fi
      echo -e "$informacoes" > $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data/Informacoes
      yad --image "dialog-question" --title "COPIA COMPLETA" --button=Nao:1 --button=Sim:0 \
--text "Deseja adicionar mais imagens ao banco de dados?"
      res=$?
    fi
  done
}

######### copiar imagens para algum lugar ##############################################################################################
copia() {
  res=0
  until [ $res = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    askcopia=$( \
      yad --form \
        --center \
        --title="INSERCAO DE IMAGENS NO BANCO DE DADOS" \
        --width=400 \
        --height=400 \
        --image="accessories-text-editor" \
        --field="Diretório onde salvar as imagens":DIR "$diretorio" \
        --field="Copiar Bias e Flats para cada objeto":CHK FALSE \
        --field="Separar por objeto":CHK FALSE \
        --field="Todos os Telescopios":CHK FALSE \
        --field="Todas as Cameras":CHK FALSE \
        --field="Separar por Camera":CHK FALSE \
        --field="Todos os Filtros":CHK FALSE \
        --field="Separar por Filtro":CHK FALSE \
        --field="Todas as Binagens":CHK FALSE \
        --field="Separar por Binagem":CHK FALSE \
        --field="Todas as datas disponiveis":CHK FALSE \
        --field="Copiar Arquivos de Reducao":CHK FALSE \
        --field="Copiar Arquivos de Analise":CHK FALSE \
      )
    diretorio=$(echo "$askcopia" | cut -d"|" -f 1) ########
    resbifl=$(echo "$askcopia" | cut -d"|" -f 2)
    ressepobj=$(echo "$askcopia" | cut -d"|" -f 3)
    restel=$(echo "$askcopia" | cut -d"|" -f 4)    ########
    rescam=$(echo "$askcopia" | cut -d"|" -f 5)    ########
    ressepcam=$(echo "$askcopia" | cut -d"|" -f 6)
    resfil=$(echo "$askcopia" | cut -d"|" -f 7)    ########
    ressepfil=$(echo "$askcopia" | cut -d"|" -f 8)
    resbin=$(echo "$askcopia" | cut -d"|" -f 9)    ########
    ressepbin=$(echo "$askcopia" | cut -d"|" -f 10)
    resdata=$(echo "$askcopia" | cut -d"|" -f 11)  ########
    resred=$(echo "$askcopia" | cut -d"|" -f 12)
    resana=$(echo "$askcopia" | cut -d"|" -f 13)
    askobjeto=$(yad --text="Selecione os objetos" --list --checklist --column selecionar \
$opcaoobjeto2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
    objeto=$(echo "$askobjeto" | cut -d"|" -f 1)
    if [ "$restel" = "TRUE" ];then
      telescopio=$(cat $telescopes)
    else
      asktelescopio=$(yad --text="Selecione o telescopio" --list --checklist --column selecionar \
$opcaotelescopio2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      telescopio=$(echo "$asktelescopio" | cut -d"|" -f 1)
    fi
    if [ "$rescam" = "TRUE" ];then
      camera=$(cat $cameras)
    else
      askcamera=$(yad --text="Selecione a camera" --list --checklist --column selecionar \
$opcaocamera2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      camera=$(echo "$askcamera" | cut -d"|" -f 1)
    fi
    if [ "$resfil" = "TRUE" ];then
      filtro=$(cat $filtros)
    else
      askfiltro=$(yad --text="Selecione o filtro" --list --checklist --column selecionar \
$opcaofiltro2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      filtro=$(echo "$askfiltro" | cut -d"|" -f 1)
    fi
    if [ "$resbin" = "TRUE" ];then
      binagem=$(cat $binagens)
    else
      askbinagem=$(yad --text="Selecione a binagem" --list --checklist --column selecionar \
$opcaobinagem2 --column opcao --hide-header --separator=" " --width=400 --height=200 --center --print-column=2)
      binagem=$(echo "$askbinagem" | cut -d"|" -f 1)
    fi
########################################################################################################
    for o in $objeto;do
      while read COL_1 COL_2 COL_3
      do
        if [ "$COL_1" = "$o" ]; then
          category="$COL_2"
          subcategory="$COL_3"
        fi
      done < $objects
      for t in $telescopio;do
        if [ -d "$atual/$category/$subcategory/$o/$t" ];then
          for c in $camera;do
            if [ -d "$atual/$category/$subcategory/$o/$t/$c" ];then
              for f in $filtro;do
                if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f" ];then
                  for b in $binagem;do
                    if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b" ];then
                      ls $atual/$category/$subcategory/$o/$t/$c/$f/$b > .auxa
                      data=$(cat .auxa)
                      rm .auxa
                      if [ $resdata = "FALSE" ];then
                        for z in $data; do
                          echo "FALSE $z" >> .auxd
                        done
                        opcaodata=$(cat .auxd)
                        rm .auxd
                        data=$(yad --text="Selecione a data: Objeto:$o; Telescopio:$t; Camera:$c; Filtro:$f Binagem:$b" \
 --list --checklist --column selecionar $opcaodata --column opcao --hide-header --separator=" " \
 --width=400 --height=200 --center --print-column=2)
                      fi
                      if [ ! -z $data ];then
                        if [ ! -d $diretorio/$o ];then mkdir $diretorio/$o; fi
                        if [ ! -d $diretorio/$o/$t ];then mkdir $diretorio/$o/$t; fi
                        if [ ! -d $diretorio/$o/$t/$c ];then mkdir $diretorio/$o/$t/$c; fi
                        if [ ! -d $diretorio/$o/$t/$c/$f ];then mkdir $diretorio/$o/$t/$c/$f; fi
                        if [ ! -d $diretorio/$o/$t/$c/$f/$b ];then mkdir $diretorio/$o/$t/$c/$f/$b; fi
                      fi
                      for i in $data;do
                        if [ ! -d $diretorio/$o/$t/$c/$f/$b/$i ];then mkdir $diretorio/$o/$t/$c/$f/$b/$i; fi
                        if [ ! -d $diretorio/$o/$t/$c/$f/$b/$i/RAW ];then mkdir $diretorio/$o/$t/$c/$f/$b/$i/RAW; fi
                        ls $atual/$category/$subcategory/$o/$t/$c/$f/$b/$i/RAW > .auxb
                        imagens=$(cat .auxb)
                        rm .auxb
                        if [ "$resred" = "TRUE" ];then
                          if [ ! -d $diretorio/$o/$t/$c/$f/$b/$i/REDUCTION ];then mkdir $diretorio/$o/$t/$c/$f/$b/$i/REDUCTION; fi
                          cp $diretorio/$o/$t/$c/$f/$b/$i/REDUCTION/* $diretorio/$o/$t/$c/$f/$b/$i/REDUCTION
                        fi
                        if [ "$resana" = "TRUE" ];then
                          
                        fi
                        for j in $imagens;do
                          echo "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$i/RAW/$j $diretorio/$o/$t/$c/$f/$b/$i/RAW" >> .auxt
                        done  
                      done
                    fi
                  done
                fi
              done
            fi
          done
        fi
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
        porcentagem=$(printf "%.2f" $porcentagemtot)
        echo "# Copiando imagens: $porcentagem %"
        echo "$porcentagemtot"
      done < .auxt) | yad --progress \
      --title="Copia" \
      --text="Copiando imagens para diretorio $diretorio..." \
      --progress-text="Copiando imagens: 0.00 %"
      --percentage=0
      rm .auxt
    fi
    yad --image "dialog-question" --title "COPIA COMPLETA" --button=Nao:1 --button=Sim:0 \
--text "Deseja copiar mais imagens para o seu diretorio?"
    res=$?
  done
}

#
#Pergunta ao usuario sobre a tarefa a ser feita ###########################################################################3###
#
fazer=1
if [ "$fazer" = "1" ]; then
  inicial=$(yad --text="O que deseja fazer?" --list --radiolist --column selecionar \
  --width=400 --height=200 --center \
  FALSE 1 "Copiar imagens para algum diretorio" FALSE 2 "Adicionar imagens no banco de dados" \
  --column opcao --column tarefa --hide-column=2 --hide-header --print-column=2)
  tarefa=$(echo "$inicial" | cut -d"|" -f 1)
  case $tarefa in
    1) copia; ;;
    2) adiciona; ;;
  esac
fi


