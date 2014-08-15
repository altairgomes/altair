atual=$(pwd)
pastas=$(cat lista_pastas)
erro=(bi BI Bi FL fla flt Fl FR Fr fr FV Fv fv fb FB Fb esc dk DK dark Dark ZR Zr zr zer Zer ZER zz ZZ);
for noite in $pastas;do
  if [ -e $noite/output ];then
    cd $noite;
    ls *.fits > temp;
    fits=(`cat temp`);
    d=1;
    until [ $d = "0" ];do
      escolha=(`echo "$(($RANDOM % ${#fits[@]}))"`);
      arquivo=${fits[$escolha]};
      d=0;
      x=0;
      while [ $x != ${#erro[@]} ];do
        echo "$arquivo" | egrep "${erro[$x]}" > /dev/null;
        if [ $? -eq 0 ]; then
          d=1;
        fi
        let "x = x + 1";
      done
    done
    text=("`grep "$arquivo" output`")
    ds9 $arquivo & OPCAO=$(zenity                                  \
      --title="LISTA DE PARAMETROS"                                \
      --text="Lista de parametros da noite: $noite \n imagem $arquivo"                \
      --list                                                       \
      --checklist --column "problema"                              \
      --column "parametros" FALSE alfa "(${text:0:14})" FALSE delta "(${text:15:13})" FALSE hora "(${text:30:11})" FALSE exptime "(${text:87:4})" FALSE data "(${text:42:10})" FALSE filter "(${text:92:21})" FALSE outro  \
      --column "valores"    \
      --width=800           \
      --height=300);
    if [ ! -z $OPCAO ];then
      comment=$(zenity --title="Comentarios" --entry --text="escreva algum comentario sobre o problema")
      echo "$noite = $OPCAO --- $comment" >> $atual/problema.txt
    fi
    rm temp;
    cd $atual;
  else
    echo "$noite = No Header" >> problema.txt;
  fi
done
