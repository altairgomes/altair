#!/bin/bash

download() {
#Pergunta ao usuario sobre o idioma
#

opcao1=$(zenity --entry --text "Escolha os idiomas que voce deseja fazer o download. De um espaco entre os numeros\n\
Voce pode pausar esse programa a qualquer momento, quando voce rodar novamente ele continuara do local onde parou\n\
Se voce pausar no meio do download de um audio, ele devera ser apagado antes de rodar novamente\n\
Provavelmente o nome desse arquivo sera item.mp3\n\
Eu aconselho a rodar o programa duas vezes no idioma selecionado\n\
Pois por um motivo desconhecido ele ja pulou alguns audios\n\n\
1) Japones\n\
2) Frances\n\
3) Mandarim\n\
4) Ingles\n\
5) Espanhol\n\
6) Alemao\n\
7) Coreano\n\
8) Russo\n\
9) Italiano\n\
10) Hindi\n\
11) Portugues Brasil" )

for language in ${opcao1[*]}; do
 case $language in
   1) idioma="ja"; pastaidioma="Japones"; aux="Ja" ;;
   2) idioma="fr-fr"; pastaidioma="Frances"; aux="Fr" ;;
   3) idioma="zh-cn"; pastaidioma="Mandarim"; aux="Ma" ;;
   4) idioma="en-us"; pastaidioma="Ingles"; aux="In" ;;
   5) idioma="es-mx"; pastaidioma="Espanhol"; aux="Es" ;;
   6) idioma="de-de"; pastaidioma="Alemao"; aux="Al" ;;
   7) idioma="ko"; pastaidioma="Coreano"; aux="Co" ;;
   8) idioma="ru"; pastaidioma="Russo"; aux="Ru" ;;
   9) idioma="it"; pastaidioma="Italiano"; aux="It" ;;
   10) idioma="hi-in"; pastaidioma="Hindi"; aux="Hi" ;;
   11) idioma="pt-br"; pastaidioma="Portugues-BR"; aux="Pt" ;;
 esac
#
# Definicao de Variaveis
#
 (
 n=0
 listacourse=(101 102 201 202)
 lim=7
 lim2=4
 unit=1
 level=1
 if [ -e "$pastaidioma" ]; then
  cd "$pastaidioma"
 else
  mkdir "$pastaidioma"
  cd "$pastaidioma"
 fi
 for course in ${listacourse[*]}; do
  if [ "$course" -gt "200" ]; then
  lim2=$[$lim2 - 1]
  level=2
  fi
  if [ -e "Course$course" ]; then
   cd "Course$course"
  else
   mkdir "Course$course"
   cd "Course$course"
  fi
  while [ $unit != $lim2 ]; do
   if [ -e "Unit$unit" ]; then
    cd "Unit$unit"
   else
    mkdir "Unit$unit"
    cd "Unit$unit"
   fi
   if [ $unit != 1 ]; then
   lim=6
   fi
   lesson=1
   while [ $lesson != $lim ]; do
    if [ -e "Lesson$lesson" ]; then
     cd "Lesson$lesson"
    else
     mkdir "Lesson$lesson"
     cd "Lesson$lesson"
    fi
    suite=1
    while [ $suite != 11 ]; do
     item=1
     while [ $item != 5 ]; do
      count=$[(suite*4)-4+$item]
      if [ ! -e "$count.mp3" ] && [ ! -e "0$count.mp3" ]; then
       if [ -e "item.mp3" ]; then
        rm item.mp3
       fi
       wget http://prod.static.livemocha.com.s3.amazonaws.com/s3/1.39/LessonPlan/$idioma/Level_$level/Course_$course/Unit_$unit/Lesson_$lesson/Suite_$suite/Item_$item/item.mp3
       if [ "$count" -lt "10" ]; then
        mv item.mp3 0$count.mp3
       else
        mv item.mp3 $count.mp3
       fi
      fi
      n=$[$n + 1]
      porcentagemtot=$(echo "$n * 100 / 2040"|bc -l)
      porcentagem=$(printf "%.2f" $porcentagemtot)
      echo "# Download $porcentagem % de $pastaidioma, Nivel:$level, Curso:$course Unidade:$unit, Licao:$lesson, audio: $count"
      echo "$porcentagemtot"
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
 ) |zenity --progress \
    --title="Download de audios do Livemocha" \
    --text="Idioma = $pastaidioma" \
    --percentage=0
done
clear
echo "Obrigado por usar este software     Divulguem!!!!!"
}

concatena(){
#Pergunta ao usuario sobre o idioma
#

opcao1=$(zenity --entry --text "Escolha os idiomas que voce deseja concatenar num so arquivo\n\
Nele tera varios arquivos da forma Tudo_01_02\n\
Onde esta 01 estara o numero da Unidade do curso e onde esta 02 estara o numero da Licao da Unidade\n\
Todos os 40 audios da Licao dos idiomas selecionados estarao em ordem no arquivo da seguinte forma:\n\
Primeiro audio da licao do primeiro idioma escolhido, primeiro audio da licao do segundo idioma escolhido\n\
Segundo audio da licao do primeiro idioma escolhido, segundo audio da licao do segundo idioma escolhido e assim por diante\n\
Voce pode colocar quantos idiomas desejar\n\n\
1) Japones\n\
2) Frances\n\
3) Mandarim\n\
4) Ingles\n\
5) Espanhol\n\
6) Alemao\n\
7) Coreano\n\
8) Russo\n\
9) Italiano\n\
10) Hindi\n\
11) Portugues Brasil" )

#
# Definicao de Variaveis
#
 (
 mkdir Tudo
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
      for language in ${opcao1[*]}; do
       case $language in
        1) idioma="Japones" ;;
        2) idioma="Frances" ;;
        3) idioma="Mandarim" ;;
        4) idioma="Ingles" ;;
        5) idioma="Espanhol" ;;
        6) idioma="Alemao" ;;
        7) idioma="Coreano" ;;
        8) idioma="Russo" ;;
        9) idioma="Italiano" ;;
        10) idioma="Hindi" ;;
        11) idioma="Portugues-BR" ;;
       esac
       if [ $unit != 10 ]; then
        if [ "$count" -lt "10" ]; then
         cat $idioma/Course$course/Unit$unit/Lesson$lesson/0$count.mp3 >> Tudo/Tudo_0"$unit"_"$lesson".mp3
        else
         cat $idioma/Course$course/Unit$unit/Lesson$lesson/$count.mp3 >> Tudo/Tudo_0"$unit"_"$lesson".mp3
        fi
       else
        if [ "$count" -lt "10" ]; then
         cat $idioma/Course$course/Unit$unit/Lesson$lesson/0$count.mp3 >> Tudo/Tudo_"$unit"_"$lesson".mp3
        else
         cat $idioma/Course$course/Unit$unit/Lesson$lesson/$count.mp3 >> Tudo/Tudo_"$unit"_"$lesson".mp3
        fi
       fi
      done
      n=$[$n + 1]
      porcentagemtot=$(echo "$n * 100 / 2040"|bc -l)
      porcentagem=$(printf "%.2f" $porcentagemtot)
      echo "# Download $porcentagem % de $pastaidioma, Nivel:$level, Curso:$course Unidade:$unit, Licao:$lesson, audio: $count"
      echo "$porcentagem"
      item=$[$item + 1]
     done
     suite=$[$suite + 1]
    done
    lesson=$[$lesson + 1]
   done
   unit=$[$unit + 1]
  done
  lim2=$[$lim2 + 3]
 done
 ) |zenity --progress \
    --title="Download de audios do Livemocha" \
    --text="Idioma = $pastaidioma" \
    --percentage=0
clear
echo "Obrigado por usar este software     Divulguem!!!!!"
}


#
#Pergunta ao usuario sobre a tarefa
#
tarefa=$(zenity --entry --text "Escolha a opção que deseja fazer\n\
1) Fazer download de audios do livemocha\n\
2) Concatenar todos os audios de uma mesma licao de um ou mais idiomas diferentes em um mesmo arquivo\n\
A opcao 2 é indicada para quem deseja passar os arquivos para um mp3 ou celular\n\
e pode ser usada pra concatenar varios idimoas alternadamente\n\
Para essa opcao é necessario que os audios ja tenham sido baixados\n" )

case $tarefa in
  1) download; ;;
  2) concatena; ;;
esac
