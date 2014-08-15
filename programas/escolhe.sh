atual=$(pwd)
pastas=$(cat lista_pastas2)
erro=(bi BI Bi FL fla flt Fl FR fC Fr fr fR FV Fv fv fb FB Fb foco esc dk DK dark Dark ZR Zr zr zer Zer ZER zz ZZ);
for noite in $pastas;do
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
    echo "$noite/$arquivo" >> $atual/lista_imagens_teste
    rm temp;
    cd $atual;
done
