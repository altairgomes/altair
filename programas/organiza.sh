#!/bin/bash

atual=$(pwd)
objects="$atual/.objects"                    #arquivo com os dados dos objeto
categoriax="$atual/.categoriax"              #arquivo com os dados das categorias
telescopes="$atual/.telescopes"              #arquivo com os dados dos telescopio
cameras="$atual/.cameras"                    #arquivo com os dados das cameras
binagens="$atual/.binagens"                  #arquivo com os dados das binagens
filtros="$atual/.filtros"                    #arquivo com os dados dos filtros

read_categoria(){
  if [ -e $categoriax ];then
    cate=$(cat $categoriax)
    for j in $cate;do
      echo -e "FALSE $j" >> .aux1
    done
    opcaocategoria=$(cat .aux1)
    rm .aux1
  fi
}
read_object(){
if [ -e $objects ];then
  while read COL_1 COL_2 COL_3;do
    echo -e "FALSE $COL_1" >> .aux2
  done < $objects
  opcaoobjeto=$(cat .aux2)
  rm .aux2
fi
}
read_telescope(){
if [ -e $telescopes ];then
  teles=$(cat $telescopes)
  for j in $teles;do
    echo -e "FALSE $j" >> .aux3
  done
  opcaotelescopio=$(cat .aux3)
  rm .aux3
fi
}
read_camera(){
if [ -e $cameras ];then
  came=$(cat $cameras)
  for j in $came;do
    echo -e "FALSE $j" >> .aux4
  done
  opcaocamera=$(cat .aux4)
  rm .aux4
fi
}
read_binagem(){
if [ -e $binagens ];then
  bina=$(cat $binagens)
  for j in $bina;do
    echo -e "FALSE $j" >> .aux5
  done
  opcaobinagem=$(cat .aux5)
  rm .aux5
fi
}
read_filtro(){
if [ -e $filtros ];then
  filt=$(cat $filtros)
  for j in $filt;do
    echo -e "FALSE $j" >> .aux6
  done
  opcaofiltro=$(cat .aux6)
  rm .aux6
fi
}

#########  inserir novo objeto  ############
insert_new() {
  new_object=$(zenity --title="Novo Objeto" --entry --text="Digite o nome do novo objeto sem espacos")
  objeto=$new_object
  category=$(zenity --text="Selecione a categoria" --list --radiolist --column selecionar \
$opcaocategoria FALSE "Nova Categoria" --column opcao --hide-header)
######### inserir nova categoria
  if [ "$category" = "Nova Categoria" ]; then
    new_category=$(zenity --title="Nova Categoria" --entry --text="Digite o nome da nova categoria sem espacos")
    category=$new_category
    mkdir $atual/$new_category
    echo $new_category >> $categoriax
  fi
  subcategoriax="$atual/$category/.categoriax"              #arquivo com os dados das subcategorias
  if [ -e $subcategoriax ];then
    subcate=$(cat $subcategoriax)
    for j in $subcate;do
      echo -e "FALSE $j" >> .aux0
    done
    opcaosubcategoria=$(cat .aux0)
    rm .aux0
  fi
  subcategoria=$(zenity --text="Selecione a subcategoria" --list --radiolist --column selecionar \
$opcaosubcategoria FALSE "Nova Subcategoria" --column opcao --hide-header)
######## inserir nova subcategoria ######################
  if [ "$subcategoria" = "Nova Subcategoria" ]; then
    new_subcategory=$(zenity --title="Nova Subcategoria" --entry --text="Digite o nome da nova subcategoria sem espacos")
    subcategoria=$new_subcategory
    mkdir $atual/$category/$new_subcategory
    echo $new_subcategory >> $subcategoriax
  fi
  mkdir $atual/$category/$subcategoria/$new_object
  echo -e "$new_object\t$category\t$subcategoria" >> $objects
}

#########  inserir novo telescopio  ############
insert_telescope() {
  new_telescope=$(zenity --title="Novo Telescopio" --entry --text="Digite o nome do novo Telescopio sem espacos")
  telescopio=$new_telescope
  echo -e "$new_telescope" >> $telescopes
}

#########  inserir nova camera  ############
insert_cam() {
  new_camera=$(zenity --title="Nova Camera" --entry --text="Digite o nome do nova camera sem espacos")
  camera=$new_camera
  echo -e "$new_camera" >> $cameras
}

#########  inserir nova binagem  ############
insert_bin() {
  new_bin=$(zenity --title="Nova Binagem" --entry --text="Digite o nome da nova binagem sem espacos")
  binagem=$new_bin
  echo -e "$new_bin" >> $binagens
}

#########  inserir novo filtro  ############
insert_fil() {
  new_filter=$(zenity --title="Novo Filtro" --entry --text="Digite o nome do novo novo filtro sem espacos")
  filtro=$new_filter
  echo -e "$new_filter" >> $filtros
}

######### copiar imagens para algum lugar ##############
copia() {
  res=0
  ret=0
  until [ $res = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    anw=1
    if [ $ret = "1" ]; then
      repetir=$(zenity --title="Pergunta" --question --text="Copiar para o mesmo diretorio: $diretorio?" \
--ok-label="SIM" --cancel-label="NAO")
      anw=$?
    fi
    if [ $anw = "1" ]; then
      diretorio=$(zenity --title="Selecione a pasta onde irá gravar as imagens" --file-selection --directory)
    fi
    objeto=$(zenity --text="Selecione o objeto" --list --checklist --column selecionar \
$opcaoobjeto --column opcao --hide-header --separator=" ")
    telescopio=$(zenity --text="Selecione o telescopio" --list --checklist --column selecionar \
$opcaotelescopio --column opcao --hide-header --separator=" ")
    camera=$(zenity --text="Selecione a camera" --list --checklist --column selecionar \
$opcaocamera --column opcao --hide-header --separator=" ")
    binagem=$(zenity --text="Selecione a binagem" --list --checklist --column selecionar \
$opcaobinagem --column opcao --hide-header --separator=" ")
    filtro=$(zenity --text="Selecione o filtro" --list --checklist --column selecionar \
$opcaofiltro --column opcao --hide-header --separator=" ")
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
                      datas=$(cat .auxa)
                      rm .auxa
                      for z in $datas; do
                        if [ -d "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$z" ];then
                          echo "FALSE $z" >> .auxd
                        fi
                      done
                      opcaodata=$(cat .auxd)
                      rm .auxd
                      data=$(zenity --text="Selecione a data: Objeto:$o; Telescopio:$t; Camera:$c; Filtro:$f Binagem:$b" \
 --list --checklist --column selecionar $opcaodata --column opcao --hide-header --separator=" ")
                      if [ ! -z $data ];then
                        if [ ! -d $diretorio/$o ];then mkdir $diretorio/$o; fi
                        if [ ! -d $diretorio/$o/$t ];then mkdir $diretorio/$o/$t; fi
                        if [ ! -d $diretorio/$o/$t/$c ];then mkdir $diretorio/$o/$t/$c; fi
                        if [ ! -d $diretorio/$o/$t/$c/$f ];then mkdir $diretorio/$o/$t/$c/$f; fi
                        if [ ! -d $diretorio/$o/$t/$c/$f/$b ];then mkdir $diretorio/$o/$t/$c/$f/$b; fi
                      fi
                      for i in $data;do
                        if [ ! -d $diretorio/$o/$t/$c/$f/$b/$i ];then mkdir $diretorio/$o/$t/$c/$f/$b/$i; fi
                        ls $atual/$category/$subcategory/$o/$t/$c/$f/$b/$i > .auxb
                        imagens=$(cat .auxb)
                        rm .auxb
                        for j in $imagens;do
                          echo "$atual/$category/$subcategory/$o/$t/$c/$f/$b/$i/$j $diretorio/$o/$t/$c/$f/$b/$i" >> .auxt
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
      done < .auxt) | zenity --progress \
      --title="Copiando..." \
      --text="Copiando imagens para diretorio $diretorio..." \
      --percentage=0
      rm .auxt
    fi
    zenity --title="Copia Completa" --question --text="Deseja copiar mais imagens para o seu diretorio?" \
--ok-label="SIM" --cancel-label="NAO"
    res=$?
    ret=1
  done
}

########## adicionar imagens ao banco de dados #################
adiciona() {
  res=0
  ret=0
  until [ $res = "1" ]; do
    read_categoria;
    read_object;
    read_telescope;
    read_camera;
    read_filtro;
    read_binagem;
    if [ $ret = "1" ]; then
      repetir=$(zenity --text="Quais informacoes deseja manter?" --list --checklist --column selecionar \
FALSE Objeto "$objeto" FALSE Telescopio "$telescopio" FALSE Camera "$camera" FALSE Filtro "$filtro" FALSE Binagem "$binagem" FALSE Data "$data" \
--column opcao --column valor --hide-header --separator=" " --multiple)
    fi
    echo "$repetir" | egrep "Objeto" > /dev/null;
    if [ $? -eq 1 ]; then
      objeto=$(zenity --text="Selecione o objeto" --list --radiolist --column selecionar \
$opcaoobjeto FALSE "Novo Objeto" --column objeto --hide-header)
      if [ "$objeto" = "Novo Objeto" ]; then
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
    echo "$repetir" | egrep "Telescopio" > /dev/null;
    if [ $? -eq 1 ]; then
      telescopio=$(zenity --text="Selecione o telescopio" --list --radiolist --column selecionar \
$opcaotelescopio FALSE "Novo Telescopio" --column opcao --hide-header)
      if [ "$telescopio" = "Novo Telescopio" ]; then
        insert_telescope;
      fi
    fi
    if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio" ]; then
      mkdir $atual/$category/$subcategory/$objeto/$telescopio
    fi
    echo "$repetir" | egrep "Camera" > /dev/null;
    if [ $? -eq 1 ]; then
      camera=$(zenity --text="Selecione a camera" --list --radiolist --column selecionar \
$opcaocamera FALSE "Nova Camera" --column opcao --hide-header)
      if [ "$camera" = "Nova Camera" ]; then
        insert_cam;
      fi
    fi
    if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera" ]; then
      mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera
    fi
    echo "$repetir" | egrep "Filtro" > /dev/null;
    if [ $? -eq 1 ]; then
      filtro=$(zenity --text="Selecione o filtro" --list --radiolist --column selecionar \
$opcaofiltro FALSE "Novo Filtro" --column opcao --hide-header)
      if [ "$filtro" = "Novo Filtro" ]; then
        insert_fil;
      fi
    fi
    if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro" ]; then
      mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro
    fi
    echo "$repetir" | egrep "Binagem" > /dev/null;
    if [ $? -eq 1 ]; then
      binagem=$(zenity --text="Selecione a binagem" --list --radiolist --column selecionar \
$opcaobinagem FALSE "Nova Binagem" --column opcao --hide-header)
      if [ "$binagem" = "Nova Binagem" ]; then
        insert_bin;
      fi
    fi
    if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem" ]; then
      mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem
    fi
    echo "$repetir" | egrep "Data" > /dev/null;
    if [ $? -eq 1 ]; then
      data=$(zenity --calendar --date-format=%F)
    fi
    if [ ! -d "$atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data" ]; then
      mkdir $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data
    fi
    arquivos=$(zenity --title="Selecao de imagens: Objeto: $objeto; Telescopio: $telescopio; Camera: $camera; Filtro: $filtro; Binagem: $binagem; Data: $data" \
--file-selection --multiple --separator=" ")
    for i in $arquivos;do
      echo $i >> .aux
    done
    count=$(wc -l < .aux)
    rm .aux
    x=0
    if [ "$count" != "0" ]; then
      (for file in $arquivos;do
        cp $file $atual/$category/$subcategory/$objeto/$telescopio/$camera/$filtro/$binagem/$data/
        x=$[$x + 1]
        porcentagemtot=$(echo "$x * 100 / $count"|bc -l)
        porcentagem=$(printf "%.2f" $porcentagemtot)
        echo "# Copiando imagens: Objeto: $objeto; Telescopio: $telescopio; Camera: $camera; Filtro: $filtro; Binagem: $binagem; Data: $data; $porcentagem %"
        echo "$porcentagemtot"
      done) | zenity --progress \
      --title="Copiando..." \
      --text="Copiando imagens: Objeto: $objeto; Telescopio: $telescopio; Camera: $camera; Filtro: $filtro; Binagem: $binagem; Data: $data; 0.00 %" \
      --percentage=0
    fi
    zenity --title="Copia Completa" --question --text="Deseja adicionar mais imagens ao banco de dados?" \
--ok-label="Sim" --cancel-label="Nao"
    res=$?
    ret=1
  done
}

#
#Pergunta ao usuario sobre a tarefa a ser feita
#
fazer=1
if [ "$fazer" = "1" ]; then
  tarefa=$(zenity --text="O que deseja fazer?" --list --radiolist --column selecionar \
  FALSE 1 "Copiar imagens para algum diretorio" FALSE 2 "Adicionar imagens no banco de dados" \
  --column opcao --column tarefa --hide-column=2 --hide-header)

  case $tarefa in
    1) copia; ;;
    2) adiciona; ;;
  esac
  if [ ! -z $data ];then
    fazer=0
  fi
fi
