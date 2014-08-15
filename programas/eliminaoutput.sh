atual=$(pwd)
pastas=$(cat lista_pastas)
for noite in $pastas;do
    cd $noite;
    gedit output;
    cd $atual;
done
