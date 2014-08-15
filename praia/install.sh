#
##### Instalar YAD #######################
#
install_yad(){
  echo "Checking YAD"
  yad=$(which yad)
  yadnewversion="0.26.1"
  yadoldversion=$(yad --version)
  if [ -z "$yad" ];then
    echo "YAD not found. Installing..."
    echo "sudo apt-get install intltool"
    sudo apt-get install intltool
    echo "sudo apt-get install gtk+-3.0"
    sudo apt-get install gtk+-3.0
    cd yad
    echo "Installing YAD"
    chmod u+x configure
    ./configure
    echo "sudo make"
    sudo make
    echo "sudo make install"
    sudo make install
    cd ..
  elif [ "$yadnewversion" != "$yadoldversion" ];then
    echo "YAD found, but it\'s old version..."
    cd yad
    echo "Upgrading YAD"
    chmod u+x configure
    echo "Configuring YAD files"
    ./configure
    echo "sudo make"
    sudo make
    echo "sudo make install"
    sudo make install
    cd ..
  else
    echo "YAD found in the newest version"
  fi
}
#
###### Instalar gfortran ###########
#
install_gfortran(){
  echo "Checking Gfortran"
  gfortran=$(which gfortran)
  if [ -z "$yad" ];then
    echo "Gfortran not found. Installing"
    echo "sudo apt-get install gfortran"
    sudo apt-get install gfortran
  else
    echo "Gfortran found"
  fi
}
#
######## Instalar PRAIA
#
install_praia(){
  diretorio="/usr/lib/praia"
  sudo echo "Copying PRAIA"
  sudo cp praia.desktop /usr/share/applications
  sudo cp praia-0.1.1.sh /usr/bin/praia
  echo "Copying icon"
  sudo cp icon_praia.png /usr/share/icons/gnome/48x48/apps
  sudo chmod a+rx /usr/bin/praia
  if [ ! -d $diretorio ];then
    sudo mkdir $diretorio
  fi
  echo "Copying Language Files"
  sudo cp -R language $diretorio
  echo "Copying Manuals"
  sudo cp praia.1 /usr/local/man/man1
  sudo cp -R manual $diretorio
  if [ ! -d "$diretorio/headers" ];then
    sudo mkdir "$diretorio/headers"
  fi
  if [ ! -d "$diretorio/script" ];then
    sudo mkdir "$diretorio/script"
  fi
  sudo cp -Rn configuration $diretorio
  cd fontes
  echo "Compiling PRAIA"
  j=$(ls *.f)
  for i in $j;do
    echo "$i"
    gfortran -O3 $i -o ${i%.f}
  done
  cd ..
  sudo cp -R fontes $diretorio
  sudo chmod a+rwx -R $diretorio
}
#
#### Inicio
#
install_yad
install_gfortran
install_praia
