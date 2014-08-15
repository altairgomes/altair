#!/bin/bash
#
###################### Definindo variaveis ################################################################################
#
export LC_NUMERIC="en_US.UTF-8"    ## define sistema numerico americano pro sistema
export LC_TIME="en_US.UTF-8"       ## define sistema de data americano pro sistema
#
export praiaversion="0.1.1"        ## define a versao do programa
export dirpadrao="/usr/lib/praia"   ## define o diretorio padrao das arquivos PRAIA
dirfontes="$dirpadrao/fontes"       ## define o diretorio onde estao as fontes do PRAIA
imagem="/usr/share/icons/gnome/48x48/apps/icon_praia.png"      ## define o caminho da imagem icone do PRAIA
diraux="$dirpadrao"             ## define um diretorio auxiliar
dirscript="$dirpadrao/script"   ## define o diretorio onde estao os scripts do usuario
dirmanual="$dirpadrao/manual"   ## define o diretorio onde estao os os manuais das tasks
dirlanguage="$dirpadrao/language"   ## define o diretorio onde estao os arquivos de idioma
dirheaders="$dirpadrao/headers"    ## define o diretorio onde serao salvos os headers
dirconf="$dirpadrao/configuration"   ## define o diretorio onde estarao os arquivos de configuracoes
arqconf="$dirconf/configuration.sh"  ## arquivos de configuracoes gerais
intervalo="                                                  |"
intervalo2="        |-----------------------------------------|"
intervalo3="  |-----------------------------------------------|"
intervalo4="                                                                                                                                  |"
fimdedat="***************************************************************************************************************************************"
#
########################## Funcao Inicio ################################################################################
#
inicio(){
  ret=1
  until [ "$ret" = "0" ];do
#
#### lendo informacoes do usuario ####
#
  arqconfini="$dirconf/inicio.sh"
  . $arqconfini
    inicioescolha="$inicioastitle!$iniciobttitle!$iniciohetitle!$iniciohenewtitle!\
$iniciohedittitle!$inicioebtitle!\
$iniciotjtitle!$iniciotstitle!$iniciounititle!$inicioscript!$iniciocftitle"
#    inicioescolha="$inicioastitle!$iniciobttitle!$iniciocotitle!$iniciogrtitle!$iniciohetitle!$iniciohenewtitle!\
#$iniciohedittitle!$iniciofstitle!$iniciojdtitle!$inicioebtitle!$iniciolctitle!$inicioostitle!$iniciopmtitle!$iniciordtitle!\
#$iniciosttitle!$iniciotjtitle!$iniciotstitle!$iniciotrtitle!$iniciounititle!$inicioscript!$iniciocftitle"
    FORM_PRAIA=$( \
      yad --form \
      --window-icon=$imagem \
      --center \
      --title="PRAIA" \
      --width=200 \
      --height=200 \
      --image="$imagem" \
      --field="PRAIA":CB "${inicioescolha//$task/^$task}" \
      --field="$iniciobtnreadme":BTN "yad --text-info --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$inicioreadme" \
      )
    task=$(echo "$FORM_PRAIA" | cut -d"|" -f 1)
#
#### Chamando as funcoes correspondentes ####
#
    case $task in
      "$inicioastitle") PRAIA_astrometry;;
      "$iniciobttitle") PRAIA_big_table;;
      "$iniciocotitle") PRAIA_coronografy;;
      "$iniciogrtitle") PRAIA_global_reduction;;
      "$iniciohetitle") PRAIA_header_extraction;;
      "$iniciohenewtitle") PRAIA_header_extraction_30;;
      "$iniciohedittitle") PRAIA_header_edit;;
      "$iniciofstitle") PRAIA_findscale;;
      "$iniciojdtitle") PRAIA_jd_gregorian;;
      "$inicioebtitle") PRAIA_JPL_ephem_batch;;
      "$iniciolctitle") PRAIA_light_curve_numerical_fit;;
      "$inicioostitle") PRAIA_occ_star_search;;
      "$iniciopmtitle") PRAIA_proper_motions;;
      "$iniciordtitle") PRAIA_radec_dec_hex;;
      "$iniciosttitle") PRAIA_statistics;;
      "$iniciotjtitle") PRAIA_targets_JPL;;
      "$iniciotstitle") PRAIA_targets_search;;
      "$iniciotrtitle") PRAIA_trajectory;;
      "$iniciounititle") PRAIA_uniform;;
      "$inicioscript") Rodar_script;;
      "$iniciocftitle") configuration;;
    esac
    if [ -z "$FORM_PRAIA" ]; then ret="0";fi
    if [ ! -z "$FORM_PRAIA" ]; then echo "export task=\"$task\"" > $arqconfini;fi
  done
}

inicio2(){
  FORM_PRAIA=$( \
    yad --form \
    --window-icon=$imagem \
    --center \
    --title="PRAIA" \
    --width=400 \
    --height=400 \
    --image="$imagem" \
    --columns=2 \
    --field="$inicioastitle":BTN "" \
    --field="$iniciobttitle":BTN "" \
    --field="$iniciocotitle":BTN "" \
    --field="$iniciogrtitle":BTN "" \
    --field="$iniciohetitle":BTN "" \
    --field="$iniciohedittitle":BTN "" \
    --field="$iniciojdtitle":BTN "" \
    --field="$inicioebtitle":BTN "" \
    --field="$iniciolctitle":BTN "" \
    --field="$inicioostitle":BTN "" \
    --field="$iniciopmtitle":BTN "" \
    --field="$iniciordtitle":BTN "" \
    --field="$iniciosttitle":BTN "" \
    --field="$iniciotjtitle":BTN "" \
    --field="$iniciotstitle":BTN "" \
    --field="$iniciotrtitle":BTN "" \
    --field="$iniciocftitle":BTN "bash -c configuration" \
    )
}
#
########################## Funcao PRAIA Astrometry ################################################################################
#
PRAIA_astrometry(){
  arqconfas="$dirconf/astrometry.sh"
 . $arqconfas
  praiaastrometry="PRAIA_astrometry_20_12"
  praiaastrometrydat="PRAIA_astrometry_20_12.dat"
  if [ "$1" = "-e" ];then
    if [ -s "$praiaastrometrydat" ];then
      roda_PRAIA $praiaastrometry $praiaastrometrydat $(pwd);
    else
      echo "dat file not found"
    fi
  else
    cancel="0"
    until [ "$cancel" = "1" ];do
      yad --plug=300 \
        --tabnum=1 \
        --form \
        --center \
        --field="$astrom2mass":RO "$twomass" \
        --field="$astromucac2":RO "$ucac2" \
        --field="$astromucac4":RO "$ucac4" \
        --field="$astromusecatask":CHK "$astckusercat" \
        --field="$astromusecat":FL "$astusercat" \
        --field="$astrom2massindex":CHK "$ast2massindex" \
        --field="$astromhead:" "$asthead" \
        --field="$astromtargets:" "$asttargets" \
        --field="$astromfdp:" "$astfdp" \
        --field="$astromfdpnum:":NUM $astfdpnum!0..10!1 \
        --field="$astrombpx:" "$astbpx" \
        --field="$astromckdat":CHK FALSE \
        --field="$astrombtnreadme":BTN "yad --text-info --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$astromreadme" &> .auxast1 &
      pid1=$!
      yad --plug=300 \
        --tabnum=2 \
        --form \
        --field="$astrom2mcommom:":NUM $ast2mcommom!0..20!1 \
        --field="$astromtargetid:":NUM $asttargetid!20..400!20 \
        --field="$astrompixsca:" "$astpixsca" \
        --field="$astrompixscaer:" "$astpixscaer" \
        --field="$astromsamepixe":CHK "$astsamepixe" \
        --field="$astrommaxcount:" "$astmaxcount" \
        --field="$astrommincount:" "$astmincount" \
        --field="$astrombtnreadme":BTN "yad --text-info --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$astromreadme" &> .auxast2 &
      pid2=$!
      yad --plug=300 \
        --tabnum=3 \
        --form \
        --field="$astrombackflat:":NUM $astbackflat!0..15!1 \
        --field="$astrombackfact:" "$astbackfact" \
        --field="$astromcroide:" "$astcroide" \
        --field="$astrompolin:":NUM $astpolin!0..3!1 \
        --field="$astromrad3":CHK "$astrad3" \
        --field="$astromrad5":CHK "$astrad5" \
        --field="$astrombtnreadme":BTN "yad --text-info --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$astromreadme" &> .auxast3 &
      pid3=$!
      yad --notebook \
        --title="$astrometrytitle" \
        --width=700 \
        --height=500 \
        --center \
        --window-icon="$imagem" \
        --key=300 --tab="$astromtab1" --tab="$astromtab2" --tab="$astromtab3"
#      while kill -0 $pid1 || kill -0 $pid2 || kill -0 $pid3; do
#        /bin/true
#      done
      astrometry1=$(cat .auxast1)
      astrometry2=$(cat .auxast2)
      astrometry3=$(cat .auxast3)
      rm .auxast1
      rm .auxast2
      rm .auxast3
      if [ -z "$astrometry1" ];then cancel="1";fi
      if [ ! -z "$astrometry1" ]; then
        if [ "$1" != "-e" ];then
#### tab1
          astckusercat=$(echo "$astrometry1" | cut -d"|" -f 4)
          astusercat=$(echo "$astrometry1" | cut -d"|" -f 5)
          ast2massindex=$(echo "$astrometry1" | cut -d"|" -f 6)
          asthead=$(echo "$astrometry1" | cut -d"|" -f 7)
          asttargets=$(echo "$astrometry1" | cut -d"|" -f 8)
          astfdp=$(echo "$astrometry1" | cut -d"|" -f 9)
          astfdpnum=$(echo "$astrometry1" | cut -d"|" -f 10)
          astbpx=$(echo "$astrometry1" | cut -d"|" -f 11)
          astckdat=$(echo "$astrometry1" | cut -d"|" -f 12)
          if [ "$astckusercat" = "TRUE" ];then astuc="1";else astuc="2";fi
          if [ "$ast2massindex" = "TRUE" ];then ast2mi="1";else ast2mi="2";fi
#### tab2
          ast2mcommom=$(echo "$astrometry2" | cut -d"|" -f 1)
          asttargetid=$(echo "$astrometry2" | cut -d"|" -f 2)
          astpixsca=$(echo "$astrometry2" | cut -d"|" -f 3)
          astpixscaer=$(echo "$astrometry2" | cut -d"|" -f 4)
          astsamepixe=$(echo "$astrometry2" | cut -d"|" -f 5)
          astmaxcount=$(echo "$astrometry2" | cut -d"|" -f 6)
          astmincount=$(echo "$astrometry2" | cut -d"|" -f 7)
          if [ "$astsamepixe" = "TRUE" ];then astsp="1";else astsp="2";fi
#### tab3
          astbackflat=$(echo "$astrometry3" | cut -d"|" -f 1)
          astbackfact=$(echo "$astrometry3" | cut -d"|" -f 2)
          astcroide=$(echo "$astrometry3" | cut -d"|" -f 3)
          astpolin=$(echo "$astrometry3" | cut -d"|" -f 4)
          astrad3=$(echo "$astrometry3" | cut -d"|" -f 5)
          astrad5=$(echo "$astrometry3" | cut -d"|" -f 6)
          if [ "$astrad3" = "TRUE" ];then astr3="3";else astr3="0";fi
          if [ "$astrad5" = "TRUE" ];then astr5="5";else astr5="0";fi
          if [ "$astckdat" = "FALSE" ];then  ## verifica se ja existe dat diretorio
            geradat_astrometry;
          fi
        #roda_PRAIA $praiaastrometry $praiaastrometrydat $d;
        #salva_astrometry;
        fi
      fi
    done
  fi
}
#
##### gerando arquivo dat ####
#
geradat_astrometry(){
  echo "$twomass${intervalo:${#twomass}} root directory under which the 2MASS catalogue sub-directories lay" > $praiaastrometrydat           ##
  echo "$ucac2${intervalo:${#ucac2}} root directory under which the UCAC2 catalogue sub-directories lay" >> $praiaastrometrydat              ##
  echo "$ucac4${intervalo:${#ucac4}} root directory under which the UCAC4 catalogue sub-directories lay" >> $praiaastrometrydat              ##
  echo "$astusercat${intervalo:${#astusercat}} User reference catalogue (PRAIA format)" >> $praiaastrometrydat    ##
  echo "$astuc${intervalo:${#astuc}} Make reduction with user's reference catalogue ? 1 - yes;  2 - no" >> $praiaastrometrydat              ##
  echo "$ast2mi${intervalo:${#ast2mi}} Use 2MASS catalog extraction acceleration index: 1 - yes; 2 - no" >> $praiaastrometrydat            ##
  echo "$asthead${intervalo:${#asthead}} extracted header data from fits images" >> $praiaastrometrydat          ##
  echo "$asttargets${intervalo:${#asttargets}} targets input file: (RA,Dec), JD, target name" >> $praiaastrometrydat    ##
  echo "$astfdp${intervalo:${#astfdp}} Field Distortion Pattern data file" >> $praiaastrometrydat            ##
  echo "$astfdpnum${intervalo:${#astfdpnum}} n nearest points for FDP computations" >> $praiaastrometrydat      ##
  echo "$astbpx${intervalo:${#astbpx}} Bad pixel mask  (xmin xmax ymin ymax)  (common to all treated images)" >> $praiaastrometrydat            ##
  echo "astrometry_photometry_praia_star11                | photometric statistics of each field" >> $praiaastrometrydat
  echo "astrometry_reduction_2mu2_praia_star11            | reduction statistics of each field for UCAC2 and 2MASS (original)" >> $praiaastrometrydat
  echo "astrometry_reduction_2mu4_praia_star11            | reduction statistics of each field for UCAC4 and 2MASS (original)" >> $praiaastrometrydat
  echo "astrometry_reduction_mp_praia_star11              | reduction statistics of each field for UCAC2 and 2MASS (t.p.+common p.m.)" >> $praiaastrometrydat
  echo "astrometry_reduction_mp_med_praia_star11          | reduction statistics of each field for UCAC2 and 2MASS (t.p.+common and non-common p.m.)" >> $praiaastrometrydat
  echo "astrometry_reduction_2mus_praia_star11            | reduction statistics of each field for User catalogue and 2MASS (original)" >> $praiaastrometrydat
  echo "astrometry_2MASS_target_praia_star11              | target statistics for 2MASS (original)" >> $praiaastrometrydat
  echo "astrometry_UCAC2_target_praia_star11              | target statistics for UCAC2" >> $praiaastrometrydat
  echo "astrometry_UCAC4_target_praia_star11              | target statistics for UCAC4" >> $praiaastrometrydat
  echo "astrometry_2MASS_target_mpu2_praia_star11         | target statistics for 2MASS (t.p.+common p.m.)" >> $praiaastrometrydat
  echo "astrometry_2MASS_target_mpu2_med_praia_star11     | target statistics for 2MASS (t.p.+common and non-common p.m.)" >> $praiaastrometrydat
  echo "astrometry_USER_target_praia_star11               | target statistics for User reference catalogue" >> $praiaastrometrydat
  echo "ucac2.red.xy                                      | extension of \"xy\" output files, UCAC2 reduction" >> $praiaastrometrydat
  echo "ucac4.red.xy                                      | extension of \"xy\" output files, UCAC4 reduction" >> $praiaastrometrydat
  echo "2mass.red.xy                                      | extension of \"xy\" output files  2MASS reduction" >> $praiaastrometrydat
  echo "2mass.rmp.xy                                      | extension of \"xy\" output files  2MASS (t.p.+common p.m.) reduction" >> $praiaastrometrydat
  echo "2mass.rme.xy                                      | extension of \"xy\" output files  2MASS (t.p.+common and non-common p.m.) reduction" >> $praiaastrometrydat
  echo "wfi.red.xy                                        | extension of \"xy\" output files, User reference catalog reduction" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} field area for tp+pm reduction (key,dg,min,arcsec): key=1-2 CCD size;  key=0 areax=dg,m,s areay=dg,m,s" >> $praiaastrometrydat
  echo "$ast2mcommom${intervalo:${#ast2mcommom}} maximum number of non-common 2MASS stars for (t.p.+common and non-common p.m.) reduction" >> $praiaastrometrydat  ##
  echo "$asttargetid${intervalo:${#asttargetid}} error radius for target identification (arcsec)" >> $praiaastrometrydat  ##
  echo "$astpixsca${intervalo:${#astpixsca}} pixel scale (arcsec/pixel)" >> $praiaastrometrydat      ##
  echo "$astpixscaer${intervalo:${#astpixscaer}} +/- error of pixel scale (arcsec/pixel)" >> $praiaastrometrydat  ##
  echo "$astsp${intervalo:${#astsp}} Automatic catalogue <-> (x,y) identification: 1 - same N-S/E-W orientation & pixel scale;  2 - mixed images" >> $praiaastrometrydat              ##
  echo "$astmaxcount${intervalo:${#astmaxcount}} ADU maximum counting cutoff for non-linear pixels (~saturation) (ex.: 32000, 65000, ...)" >> $praiaastrometrydat  ##
  echo "$astmincount${intervalo:${#astmincount}} ADU minimum counting cutoff for non-linear pixels (~sky background) (ex.: -10, 0 , ...)" >> $praiaastrometrydat  ##
  echo "$twomass${intervalo:${#twomass}} Pixel physical counts: 0 = from image header or 1 = from user (here)" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} Pixel physical counts: bscale;  Pixel = bscale * matrix + bzero" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} Pixel physical counts: bzero ;  Pixel = bscale * matrix + bzero" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} bitpix: -99 reads from image header; otherwise use 16, 32, 64, -32, -64 following FITS conventions" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} litteendian x bigendian: byte-swap (0 = automatic; 1 = don't swap ; 2 = swap bytes)" >> $praiaastrometrydat
  echo "$astbackflat${intervalo:${#astbackflat}} sky background flattening: degree of complete bi-variate 2-D polynomial model (1 - 15) (0 = no flattening)" >> $praiaastrometrydat  ##
  echo "$twomass${intervalo:${#twomass}} N for smoothing filter of (2N+1) channels in sky background computations (0 = no filter) (suggestion: N=5)" >> $praiaastrometrydat
  echo "$astbackfact${intervalo:${#astbackfact}} sky background theshold factor: theshold = sky + FACTOR * sigma (objects ID)" >> $praiaastrometrydat  ##
  echo "$twomass${intervalo:${#twomass}} minimum maximum FWHM (range of FWHMs) (objects ID)" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} No. of iterations in trace image fitting" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} No. of brightest 2MASS stars for cross-identification with brightest measured (x,y) objects" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} No. of brightest measured (x,y) objects for cross-identification with brightest 2MASS stars" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} (RA, Dec) area for brightest 2MASS stars search (arcmin)" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} error radius (arcsec) for cross-identification between brightest catalogue/measured objs" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} (O-C) cutoff for outliers in (RA,DEC) reductions with 2MASS (original)" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} (O-C) cutoff for outliers in (RA,DEC) reductions with UCAC2, UCAC4 & 2MASS corrected versions" >> $praiaastrometrydat
  echo "$twomass${intervalo:${#twomass}} (O-C) cutoff for outliers in (RA,DEC) reductions with user's reference catalogue" >> $praiaastrometrydat
  echo "$astpolin${intervalo:${#astpolin}} polynomial (x,y) <-> (X,Y) in (RA,DEC) reductions: 0 = 4 Ctes; 1 to 3 = complete order" >> $praiaastrometrydat        ##
  echo "$astr3${intervalo:${#astr3}} radial distortion of 3rd order (x,y) <-> (X,Y) in (RA,DEC) reduction: 0 = no; 3 = yes" >> $praiaastrometrydat              ##
  echo "$astr5${intervalo:${#astr5}} radial distortion of 5th order (x,y) <-> (X,Y) in (RA,DEC) reduction: 0 = no; 5 = yes" >> $praiaastrometrydat              ##
  echo "00 00                                             | range of fits images to reduce from extracted header data file" >> $praiaastrometrydat
  echo "$fimdedat" >> $praiaastrometrydat
}
#
#### Salvando variaveis ####
#
salva_astrometry(){
  echo "export astckusercat=\"$astckusercat\"" > $arqconfas
  echo "export astusercat=\"$astusercat\"" >> $arqconfas
  echo "export ast2massindex=\"$ast2massindex\"" >> $arqconfas
  echo "export asthead=\"$asthead\"" >> $arqconfas
  echo "export asttargets=\"$asttargets\"" >> $arqconfas
  echo "export astfdp=\"$astfdp\"" >> $arqconfas
  echo "export astfdpnum=\"$astfdpnum\"" >> $arqconfas
  echo "export astbpx=\"$astbpx\"" >> $arqconfas
  echo "export ast2mcommom=\"$ast2mcommom\"" >> $arqconfas
  echo "export asttargetid=\"$asttargetid\"" >> $arqconfas
  echo "export astpixsca=\"$astpixsca\"" >> $arqconfas
  echo "export astpixscaer=\"$astpixscaer\"" >> $arqconfas
  echo "export astsamepixe=\"$astsamepixe\"" >> $arqconfas
  echo "export astmaxcount=\"$astmaxcount\"" >> $arqconfas
  echo "export astmincount=\"$astmincount\"" >> $arqconfas
  echo "export astbackflat=\"$astbackflat\"" >> $arqconfas
  echo "export astbackfact=\"$astbackfact\"" >> $arqconfas
  echo "export astcroide=\"$astcroide\"" >> $arqconfas
  echo "export astpolin=\"$astpolin\"" >> $arqconfas
  echo "export astrad3=\"$astrad3\"" >> $arqconfas
  echo "export astrad5=\"$astrad5\"" >> $arqconfas
}
#
########################## Funcao PRAIA Big Table ################################################################################
#
PRAIA_big_table(){
  arqconfbt="$dirconf/bigtable.sh"
  . $arqconfbt
  praiabigtable="PRAIA_big_table_20_02"
  praiabigtabledat="PRAIA_big_table_20_02.dat"
  if [ "$1" = "-e" ];then
    if [ -s "$praiabigtabledat" ];then
      roda_PRAIA $praiabigtable $praiabigtabledat $(pwd);
    else
      echo "dat file not found"
    fi
  elif [ "$1" = "-r" ];then
    geradat_big_table;
    roda_PRAIA $praiabigtable $praiabigtabledat $(pwd);
  else
    cancel="0"
    until [ "$cancel" = "1" ];do
      bigtable=$( \
      yad --form \
        --center \
        --title="$bigttitle" \
        --width=650 \
        --height=500 \
        --window-icon=$imagem \
        --scroll \
        --field="$bigtoutput" "$btboutput" \
        --field="$bigttgoutput" "$btbtgoutput" \
        --field="$bigtsigthr" "$btbsigthr" \
        --field="$bigtsigfac" "$btbsigfac" \
        --field="$bigttimeint" "$btbtimeint" \
        --field="$bigttgname" "$btbtgname" \
        --field="$bigtpixscl" "$btbpixscl" \
        --field="$bigtinst" "$btbinst" \
        --field="$bigtobserver" "$btbobserver" \
        --field="$bigttreatm" "$btbtreatm" \
        --field="$bigtmodel" "$btbmodel" \
        --field="$bigtcatalo" "$btbcatalo" \
        --field="$bigtephem" "$btbephem" \
        --field="$bigtnotes" "$btbnotes" \
        --field="$bigtdirectory":DIR "$btbdirectory" \
        --field="$bigtsubdirectory":CHK "FALSE" \
        --field="$bigtchkdattxt":CHK "$btbchkdat" \
        --field="$bigtbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$bigtreadme" \
      )
      btboutput=$(echo "$bigtable" | cut -d"|" -f 1)
      btbtgoutput=$(echo "$bigtable" | cut -d"|" -f 2)
      btbsigthr=$(echo "$bigtable" | cut -d"|" -f 3)
      btbsigfac=$(echo "$bigtable" | cut -d"|" -f 4)
      btbtimeint=$(echo "$bigtable" | cut -d"|" -f 5)
      btbtgname=$(echo "$bigtable" | cut -d"|" -f 6)
      btbpixscl=$(echo "$bigtable" | cut -d"|" -f 7)
      btbinst=$(echo "$bigtable" | cut -d"|" -f 8)
      btbobserver=$(echo "$bigtable" | cut -d"|" -f 9)
      btbtreatm=$(echo "$bigtable" | cut -d"|" -f 10)
      btbmodel=$(echo "$bigtable" | cut -d"|" -f 11)
      btbcatalo=$(echo "$bigtable" | cut -d"|" -f 12)
      btbephem=$(echo "$bigtable" | cut -d"|" -f 13)
      btbnotes=$(echo "$bigtable" | cut -d"|" -f 14)
      btbdirectory=$(echo "$bigtable" | cut -d"|" -f 15)
      btbsubdirectory=$(echo "$bigtable" | cut -d"|" -f 16)
      btbchkdat=$(echo "$bigtable" | cut -d"|" -f 17)
      if [ -z "$bigtable" ];then cancel="1";fi
      if [ ! -z "$bigtable" ]; then  ## verifica se usuario nao clicou em cancel
        if [ "$btbsubdirectory" = "TRUE" ];then  ## usuario quer rodar em subdiretorios recursivamente que contenham fits
          busca="$btbtgoutput"
          lista_pastas $btbdirectory;
          btbdirectory=$(cat $dirpadrao/.lista_pastas)
          rm $dirpadrao/.lista_pastas
        fi
        (for d in $btbdirectory;do ## para toda pasta que existir fits
          cd $d
          if [ "$btbchkdat" = "FALSE" ];then  ## verifica se ja existe dat no diretorio
            geradat_big_table;
          fi
          roda_PRAIA $praiabigtable $praiabigtabledat $d;
        done) | yad --progress \
            --center \
            --title "Big Table" \
            --width="400" \
            --window-icon=$imagem \
            --progress-text="$bigtprogtxt" \
            --pulsate --auto-close
        salva_big_table;
      fi
    done
  fi
}
#
##### gerando arquivo dat ####
#
geradat_big_table(){
  btboutfil="$btbtgoutput""_filtered"
  btbouteli="$btbtgoutput""_eliminated"
  btbdate=$(date +"%d %b %Y")
  btbver="PRAIA V.$praiaversion"
  echo "$btboutput${intervalo4:${#btboutput}} Output Table file name" > $praiabigtabledat
  echo "$fimdedat" >> $praiabigtabledat
  echo "$btbtgoutput${intervalo4:${#btbtgoutput}} PRAIA output target file for averaging" >> $praiabigtabledat
  echo "$btboutfil${intervalo4:${#btboutfil}} output file: filtered PRAIA target file" >> $praiabigtabledat
  echo "$btbouteli${intervalo4:${#btbouteli}} output file: eliminated PRAIA target file" >> $praiabigtabledat
  echo "$btbsigthr${intervalo4:${#btbsigthr}} Position offset: sigma threshold (mas)" >> $praiabigtabledat
  echo "$btbsigfac${intervalo4:${#btbsigfac}} Outliers; sigma factor f: offset-mean > f*s" >> $praiabigtabledat
  echo "$btbtimeint${intervalo4:${#btbtimeint}} Grouping time interval (yr,mo,day,hh,mm,ss)" >> $praiabigtabledat
  echo "$btbtgname${intervalo4:${#btbtgname}} target name / IAU code" >> $praiabigtabledat
  echo "$btbpixscl${intervalo4:${#btbpixscl}} pixel scale (arcsec per pixel)" >> $praiabigtabledat
  echo "$btbinst${intervalo4:${#btbinst}} instrument / observatory / IAU code" >> $praiabigtabledat
  echo "$btbobserver${intervalo4:${#btbobserver}} observer" >> $praiabigtabledat
  echo "$btbtreatm${intervalo4:${#btbtreatm}} data treatment" >> $praiabigtabledat
  echo "$btbver${intervalo4:${#btbver}} package" >> $praiabigtabledat
  echo "$btbmodel${intervalo4:${#btbmodel}} model (x,y) <---> (RA,DEC) at tangent plane" >> $praiabigtabledat
  echo "$btbcatalo${intervalo4:${#btbcatalo}} reference catalogue" >> $praiabigtabledat
  echo "$btbephem${intervalo4:${#btbephem}} reference ephemeris or database" >> $praiabigtabledat
  echo "$btbdate${intervalo4:${#btbdate}} date of release in Table" >> $praiabigtabledat
  echo "$btbnotes${intervalo4:${#btbnotes}} general notes" >> $praiabigtabledat
  echo "$fimdedat" >> $praiabigtabledat
}
#
#### Salvando variaveis ####
#
salva_big_table(){
  echo "export  btboutput=\"$btboutput\"" > $arqconfbt
  echo "export  btbtgoutput=\"$btbtgoutput\"" >> $arqconfbt
  echo "export  btbsigthr=\"$btbsigthr\"" >> $arqconfbt
  echo "export  btbsigfac=\"$btbsigfac\"" >> $arqconfbt
  echo "export  btbtimeint=\"$btbtimeint\"" >> $arqconfbt
  echo "export  btbtgname=\"$btbtgname\"" >> $arqconfbt
  echo "export  btbpixscl=\"$btbpixscl\"" >> $arqconfbt
  echo "export  btbinst=\"$btbinst\"" >> $arqconfbt
  echo "export  btbobserver=\"$btbobserver\"" >> $arqconfbt
  echo "export  btbtreatm=\"$btbtreatm\"" >> $arqconfbt
  echo "export  btbmodel=\"$btbmodel\"" >> $arqconfbt
  echo "export  btbcatalo=\"$btbcatalo\"" >> $arqconfbt
  echo "export  btbephem=\"$btbephem\"" >> $arqconfbt
  echo "export  btbnotes=\"$btbnotes\"" >> $arqconfbt
}
#
########################## Funcao PRAIA Coronography ################################################################################
#
PRAIA_coronografy(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Global Reduction ################################################################################
#
PRAIA_global_reduction(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Header Extraction ################################################################################
#
PRAIA_header_extraction(){
  praiahextract="PRAIA_header_extraction_20_05"
  praiahextractdat="PRAIA_header_extraction_20_05.dat"
  cancel="0"
  until [ "$cancel" = "1" ];do
    header=$( \
    yad --form \
      --center \
      --title="$hextractitle" \
      --width=400 \
      --height=300 \
      --window-icon=$imagem \
      --field="$hextractoutput:" "output" \
      --field="$hextractfitsh:":CB "LNA CCD 301!ESO 1.2m Swiss CCD!ESO2p2/WFI mosaic!\
SOAR 4m/SOI CCDs!LNA SBIG-STL!CFHT MegaCam mosaic!Tubitak Observatory!Las Campanas 40cm Astrograph!\
LNA_S800!LNA IKON, IXON!Merlin/Raptor!Pic du Midi T1M (Astrometrica)!PRECAM - DES Project!Pic du Midi T1M - debug"  \
      --field="$hextradirectory":DIR "$hexdirectory" \
      --field="$hextrasubdirectory":CHK FALSE \
      --field="$hextrachkdattxt":CHK FALSE \
      --field="$hextrabtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$hextrareadme" \
    )
    if [ -z "$header" ];then cancel="1";fi
    if [ ! -z "$header" ]; then  ## verifica se usuario nao clicou em cancel
      hexoutput=$(echo "$header" | cut -d"|" -f 1)
      hexfitsheader=$(echo "$header" | cut -d"|" -f 2)
      hexdirectory=$(echo "$header" | cut -d"|" -f 3)
      hexsubdirectory=$(echo "$header" | cut -d"|" -f 4)
      hextrachkdat=$(echo "$header" | cut -d"|" -f 5)
      hexdirectory2="$hexdirectory"
      if [ "$hexsubdirectory" = "TRUE" ];then  ## usuario quer rodar em subdiretorios recursivamente que contenham fits
        busca="*.fits"
        lista_pastas $hexdirectory;
        hexdirectory2=$(cat $dirpadrao/.lista_pastas)
        rm $dirpadrao/.lista_pastas
      fi
      (for d in $hexdirectory2;do ## para toda pasta que existir fits
        cd $d
        ls *.fits > lista_fits
        if [ "$hextrachkdat" = "FALSE" ];then  ## verifica se ja existe dat no diretorio
          case $hexfitsheader in ## verifica tipo de header do fits
            "LNA CCD 301") hexfitshvalue="01";;
            "ESO 1.2m Swiss CCD") hexfitshvalue="02";;
            "ESO2p2/WFI mosaic") hexfitshvalue="03";;
            "SOAR 4m/SOI CCDs") hexfitshvalue="04";;
            "LNA SBIG-STL") hexfitshvalue="05";;
            "CFHT MegaCam mosaic") hexfitshvalue="06";;
            "Tubitak Observatory") hexfitshvalue="07";;
            "Las Campanas 40cm Astrograph") hexfitshvalue="08";;
            "LNA_S800") hexfitshvalue="09";;
            "LNA IKON, IXON") hexfitshvalue="10";;
            "Merlin/Raptor") hexfitshvalue="11";;
            "Pic du Midi T1M (Astrometrica)") hexfitshvalue="12";;
            "PRECAM - DES Project") hexfitshvalue="13";;
            "Pic du Midi T1M - debug") hexfitshvalue="14";;
          esac
#
##### gerando arquivo dat ####
#
          echo "lista_fits                                        |" > $praiahextractdat
          echo "$hexoutput ${intervalo:${#hexoutput}}" >> $praiahextractdat
          echo "$hexfitshvalue ${intervalo:${#hexfitshvalue}}" >> $praiahextractdat
          echo "$fimdedat" >> $praiahextractdat
        fi
#
#### Rodando PRAIA ####
#
        roda_PRAIA $praiahextract $praiahextractdat $d
      done) | yad --progress \
          --title "Header Extraction" \
          --center \
          --width="400" \
          --window-icon=$imagem \
          --progress-text="$hextraprogtxt" \
          --pulsate --auto-close
    fi
  done
}
#
########################## Funcao Adicionar Headers ################################################################################
#
add_header(){
  adicionahe=$(yad --form \
    --title "$addhetitle" \
    --center \
    --width=600 \
    --height=600 \
    --window-icon=$imagem \
    --scroll \
    --button="$sair":1 --button="$addhebtadd":2 \
    --field="$addhename" "" \
    --field="$addheax1" "NAXIS1" \
    --field="$addheax2" "NAXIS2" \
    --field="$addheobj" "OBJECT" \
    --field="$addhefil" "FILTER" \
    --field="$addhera" "RA" \
    --field="$addhedec" "DEC" \
    --field="$addhejd" "JD" \
    --field="$addhemjd" "JDM" \
    --field="$addheodt" "DATE-OBS" \
    --field="$addhesta" "UT" \
    --field="$addheedn" "TIME-END" \
    --field="$addheseexp" "TIME" \
    --field="$addheexp" "EXPTIME" \
    --field="*$addhetsk":NUM 1!1..4!1 \
    --field="*$addheccd":CB "Single CCD frame!Mosaic CCD frame" \
    --field="*$addherask":NUM 1!1..5!1 \
    --field="*$addhedesk":NUM 1!1..3!1 \
    --field="$addheraco":CB "HOUR!DEG!RAD" \
    --field="$addhedeco":CB "DEG!RAD!HOUR" \
    --field="$newhexbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$newhexreadme" \
  )
  if [ ! -z "$adicionahe" ];then
    addhexnm=$(echo "$adicionahe" | cut -d"|" -f 1)
    addhexa1=$(echo "$adicionahe" | cut -d"|" -f 2)
    addhexa2=$(echo "$adicionahe" | cut -d"|" -f 3)
    addhexob=$(echo "$adicionahe" | cut -d"|" -f 4)
    addhexfl=$(echo "$adicionahe" | cut -d"|" -f 5)
    addhexra=$(echo "$adicionahe" | cut -d"|" -f 6)
    addhexde=$(echo "$adicionahe" | cut -d"|" -f 7)
    addhexjd=$(echo "$adicionahe" | cut -d"|" -f 8)
    addhexmjd=$(echo "$adicionahe" | cut -d"|" -f 9)
    addhexodt=$(echo "$adicionahe" | cut -d"|" -f 10)
    addhexsta=$(echo "$adicionahe" | cut -d"|" -f 11)
    addhexend=$(echo "$adicionahe" | cut -d"|" -f 12)
    addhexsee=$(echo "$adicionahe" | cut -d"|" -f 13)
    addhexexp=$(echo "$adicionahe" | cut -d"|" -f 14)
    addhextsk=$(echo "$adicionahe" | cut -d"|" -f 15)
    read addhextsk <<< $(echo $addhextsk | cut -d . -f 1)
    addhexccd=$(echo "$adicionahe" | cut -d"|" -f 16)
    case $addhexccd in
      "Single CCD frame") addhexccdtyp="1";;
      "Mosaic CCD frame") addhexccdtyp="2";;
    esac
    addhexrask=$(echo "$adicionahe" | cut -d"|" -f 17)
    read addhexrask <<< $(echo $addhexrask | cut -d . -f 1)
    addhexdesk=$(echo "$adicionahe" | cut -d"|" -f 18)
    read addhexdesk <<< $(echo $addhexdesk | cut -d . -f 1)
    addhexraco=$(echo "$adicionahe" | cut -d"|" -f 19)
    case $addhexraco in
      "DEG") addhexratyp="1";;
      "RAD") addhexratyp="2";;
      "HOUR") addhexratyp="3";;
    esac
    addhexdeco=$(echo "$adicionahe" | cut -d"|" -f 20)
    case $addhexdeco in
      "DEG") addhexdetyp="1";;
      "RAD") addhexdetyp="2";;
      "HOUR") addhexdetyp="3";;
    esac
    echo "export nhea1=\"$addhexa1\"" > $dirheaders/"$addhexnm".sh
    echo "export nhea2=\"$addhexa2\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheob=\"$addhexob\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhefl=\"$addhexfl\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhera=\"$addhexra\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhede=\"$addhexde\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhejd=\"$addhexjd\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhemjd=\"$addhexmjd\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheodt=\"$addhexodt\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhesta=\"$addhexsta\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheend=\"$addhexend\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhesee=\"$addhexsee\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheexp=\"$addhexexp\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhetsk=\"$addhextsk\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheccd=\"$addhexccd\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheccdtyp=\"$addhexccdtyp\"" >> $dirheaders/"$addhexnm".sh
    echo "export nherask=\"$addhexrask\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhedesk=\"$addhexdesk\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheraco=\"$addhexraco\"" >> $dirheaders/"$addhexnm".sh
    echo "export nheratyp=\"$addhexratyp\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhedeco=\"$addhexdeco\"" >> $dirheaders/"$addhexnm".sh
    echo "export nhedetyp=\"$addhexdetyp\"" >> $dirheaders/"$addhexnm".sh
    echo -e "$addhexnm" >> $listaheaders
    var=$(sort $listaheaders)
    echo "$var" > $listaheaders
  fi
}
#
########################## Funcao Header Extraction NOVO ################################################################################
#
PRAIA_header_extraction_30(){
  arqconfnhe="$dirconf/newhex.sh"
  . $arqconfnhe
  praiahexnew="PRAIA_header_extraction_30_02"
  praiahexnewdat="PRAIA_header_extraction_30_02.dat"
  listaheaders="$dirheaders/.lista_headers"
  if [ "$1" = "-e" ];then
    if [ -s "$praiahexnewdat" ];then
      d=$(pwd)
      roda_PRAIA $praiahexnew $praiahexnewdat $d;
    else
      echo "file not found: $praiahexnewdat"
    fi
  elif [ "$1" = "-r" ];then
    . $dirheaders/"$nhefitsheader".sh
    geradat_header_extraction_30;
    roda_PRAIA $praiahexnew $praiahexnewdat $(pwd);
  else
    cancel="0"
    until [ "$cancel" = "1" ];do
      run="0"
      newhexoption=""
      while read nhnm;do
        newhexoption="$newhexoption$nhnm!"
      done < $listaheaders
      newhexoption=$(echo ${newhexoption%!*})
      nexhex=$(yad --form \
        --title "$newhextitle" \
        --center \
        --width=400 \
        --height=300 \
        --window-icon=$imagem \
        --button="$sair":1 --button="$newhexbtadd":2 --button="$newhexrodar":0 \
        --field="$newhexoutput:" "$nheoutput" \
        --field="$newhexfitsh":CB "${newhexoption//$nhefitsheader/^$nhefitsheader}" \
        --field="$newhextoff" "$nhetimeoff" \
        --field="$newhexdirectory":DIR "$nhedirectory" \
        --field="$newhexsubdirectory":CHK FALSE \
        --field="$newhexchkdattxt":CHK FALSE \
        --field="$newhexbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$newhexreadme" \
      )
      res=$?
      nheoutput=$(echo "$nexhex" | cut -d"|" -f 1)
      nhefitsheader=$(echo "$nexhex" | cut -d"|" -f 2)
      read nhefitsheader <<< $nhefitsheader
      nhetimeoff=$(echo "$nexhex" | cut -d"|" -f 3)
      nhedirectory=$(echo "$nexhex" | cut -d"|" -f 4)
      nhedirectory2="$nhedirectory"
      nhesubdirectory=$(echo "$nexhex" | cut -d"|" -f 5)
      nhechkdat=$(echo "$nexhex" | cut -d"|" -f 6)
      if [ "$nhesubdirectory" = "TRUE" ];then  ## usuario quer rodar em subdiretorios recursivamente que contenham fits
        busca="*.fits"
        lista_pastas $nhedirectory;
        nhedirectory2=$(cat $dirpadrao/.lista_pastas)
        rm $dirpadrao/.lista_pastas
      fi
      if [ "$nhechkdat" = "FALSE" ];then  ## verifica se ja existe dat no diretorio
        . $dirheaders/"$nhefitsheader".sh
      fi
      case $res in
        0) cancel="0"; run="1";;
        1) cancel="1";;
        2) cancel="0"; add_header;;
      esac
      if [ "$run" = "1" ];then
        (for d in $nhedirectory2;do ## para toda pasta que existir fits
          geradat_header_extraction_30;
          roda_PRAIA $praiahexnew $praiahexnewdat $d;
        done) | yad --progress \
            --title "Header Extraction" \
            --center \
            --width="400" \
            --window-icon=$imagem \
            --progress-text="$hextraprogtxt" \
            --pulsate --auto-close
        salva_header_extraction_30;
      fi
    done
  fi
}
#
##### gerando arquivo dat ####
#
geradat_header_extraction_30(){
  cd $d
  if [ ! -s lista_fits ];then
    ls *.fits > lista_fits
  fi
  echo "lista_fits                                        |" > $praiahexnewdat
  echo "$nheoutput${intervalo:${#nheoutput}}" >> $praiahexnewdat
  echo "$nhea1${intervalo2:${#nhea1}}" >> $praiahexnewdat
  echo "$nhea2${intervalo2:${#nhea2}}" >> $praiahexnewdat
  echo "$nheob${intervalo2:${#nheob}}" >> $praiahexnewdat
  echo "$nhefl${intervalo2:${#nhefl}}" >> $praiahexnewdat
  echo "$nhera${intervalo2:${#nhera}}" >> $praiahexnewdat
  echo "$nhede${intervalo2:${#nhede}}" >> $praiahexnewdat
  echo "$nhejd${intervalo2:${#nhejd}}" >> $praiahexnewdat
  echo "$nhemjd${intervalo2:${#nhemjd}}" >> $praiahexnewdat
  echo "$nheodt${intervalo2:${#nheodt}}" >> $praiahexnewdat
  echo "$nhesta${intervalo2:${#nhesta}}" >> $praiahexnewdat
  echo "$nheend${intervalo2:${#nheend}}" >> $praiahexnewdat
  echo "$nhesee${intervalo2:${#nhesee}}" >> $praiahexnewdat
  echo "$nheexp${intervalo2:${#nheexp}}" >> $praiahexnewdat
  echo "$nhetimeoff${intervalo:${#nhetimeoff}}" >> $praiahexnewdat
  echo "$nhetsk${intervalo3:${#nhetsk}}" >> $praiahexnewdat
  echo "$nheccdtyp${intervalo3:${#nheccdtyp}}" >> $praiahexnewdat
  echo "$nherask${intervalo3:${#nherask}}" >> $praiahexnewdat
  echo "$nhedesk${intervalo3:${#nhedesk}}" >> $praiahexnewdat
  echo "$nheratyp${intervalo3:${#nheratyp}}" >> $praiahexnewdat
  echo "$nhedetyp${intervalo3:${#nhedetyp}}" >> $praiahexnewdat
  echo "$fimdedat" >> $praiahexnewdat
}
#
#### Salvando variaveis ####
#
salva_header_extraction_30(){
  echo "export nheoutput=\"$nheoutput\"" > $arqconfnhe
  echo "export nhetimeoff=\"$nhetimeoff\"" >> $arqconfnhe
  echo "export nhefitsheader=\"$nhefitsheader\"" >> $arqconfnhe
}
#
########################## Funcao PRAIA Header Edit ################################################################################
#
PRAIA_header_edit(){
  arqconfhe="$dirconf/headeredit.sh"
  . $arqconfhe
  praiahedit="PRAIA_header_edit_20_01"
  praiaheditdat="PRAIA_header_edit_20_01.dat"
  if [ "$1" = "-e" ];then
    if [ -s "$praiaheditdat" ];then
      roda_PRAIA $praiahedit $praiaheditdat $(pwd);
    else
      echo "dat file not found"
    fi
  elif [ "$1" = "-r" ];then
    geradat_header_edit;
    roda_PRAIA $praiahedit $praiaheditdat $(pwd);
  else
    cancel="0"
    until [ "$cancel" = "1" ];do
#
#### lendo informacoes do usuario ####
#
      yad --plug=100 \
      --tabnum=1 \
      --form \
      --field="$heditoutputin:" "$hedoutin" \
      --field="$heditoutputou:" "$hedoutou" \
      --field="$heditaskname":CHK "$hedoutana" \
      --field="$heditaskfilt":CHK "$hedoutafl" \
      --field="$heditaskexpo":CHK "$hedoutaex" \
      --field="$heditaskano":CHK "$hedoutaan" \
      --field="$heditaskoff":CHK "$hedoutaof" \
      --field="$heditmulti":CHK "$hedoutamu" \
      --field="$heditdivide":CHK "$hedoutadi" \
      --field="$heditprecess":CHK "$hedoutapr" \
      --field="$heditaskephe":CHK "$hedoutaep" \
      --field="$heditdirectory":DIR "$hedoutadir" \
      --field="$heditdattxt":CHK FALSE \
      --field="$heditbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$heditreadme" &> .auxhed1 &
      pid1=$!
      yad --plug=100 \
      --tabnum=2 \
      --form \
      --text="$heditvatxt" \
      --field="$heditname:" "$hedoutnm" \
      --field="$heditfilt:" "$hedoutfi" \
      --field="$heditexpo:" "$hedoutep" \
      --field="$heditano:" "$hedoutyr" \
      --field="$heditoff:" "$hedoutof" \
      --field="$hedittargets:" "$hedouttg" \
      --field="$hedittargerr:" "$hedoutte" &> .auxhed2 &
      pid2=$!
      yad --notebook \
      --title="$hedittitle" \
      --center \
      --width=400 \
      --height=600 \
      --window-icon=$imagem \
      --key=100 --tab="$hedittab1" --tab="$hedittab2"
      while kill -0 $pid1 || kill -0 $pid2; do
        /bin/true
      done
      headeredit1=$(cat .auxhed1)
      headeredit2=$(cat .auxhed2)
      rm .auxhed1
      rm .auxhed2
      if [ ! -z "$headeredit1" ];then
        hedoutin=$(echo "$headeredit1" | cut -d"|" -f 1)
        hedoutou=$(echo "$headeredit1" | cut -d"|" -f 2)
        hedoutana=$(echo "$headeredit1" | cut -d"|" -f 3)
        hedoutafl=$(echo "$headeredit1" | cut -d"|" -f 4)
        hedoutaex=$(echo "$headeredit1" | cut -d"|" -f 5)
        hedoutaan=$(echo "$headeredit1" | cut -d"|" -f 6)
        hedoutaof=$(echo "$headeredit1" | cut -d"|" -f 7)
        hedoutamu=$(echo "$headeredit1" | cut -d"|" -f 8)
        hedoutadi=$(echo "$headeredit1" | cut -d"|" -f 9)
        hedoutapr=$(echo "$headeredit1" | cut -d"|" -f 10)
        hedoutaep=$(echo "$headeredit1" | cut -d"|" -f 11)
        hedoutadir=$(echo "$headeredit1" | cut -d"|" -f 12)
        hedoutdat=$(echo "$headeredit1" | cut -d"|" -f 13)
        if [ "$hedoutana" = "FALSE" ];then hedn=0;else hedn=1;fi
        if [ "$hedoutafl" = "FALSE" ];then hedf=0;else hedf=1;fi
        if [ "$hedoutaex" = "FALSE" ];then hede=0;else hede=1;fi
        if [ "$hedoutaan" = "FALSE" ];then heda=0;else heda=1;fi
        if [ "$hedoutaof" = "FALSE" ];then hedo=0;else hedo=1;fi
        if [ "$hedoutamu" = "FALSE" ];then hedm=0;else hedm=1;fi
        if [ "$hedoutadi" = "FALSE" ];then hedd=0;else hedd=1;fi
        if [ "$hedoutapr" = "FALSE" ];then hedp=0;else hedp=1;fi
        if [ "$hedoutaep" = "FALSE" ];then hedh=0;else hedh=1;fi
        hedoutnm=$(echo "$headeredit2" | cut -d"|" -f 1)
        hedoutfi=$(echo "$headeredit2" | cut -d"|" -f 2)
        hedoutep=$(echo "$headeredit2" | cut -d"|" -f 3)
        hedoutyr=$(echo "$headeredit2" | cut -d"|" -f 4)
        hedoutof=$(echo "$headeredit2" | cut -d"|" -f 5)
        hedouttg=$(echo "$headeredit2" | cut -d"|" -f 6)
        hedoutte=$(echo "$headeredit2" | cut -d"|" -f 7)
        cd $hedoutadir
        if [ "$hedoutdat" = "FALSE" ];then
          geradat_header_edit;
        fi
        (roda_PRAIA $praiahedit $praiaheditdat $hedoutadir;
        ) | yad --progress \
          --title "Header Edit" \
          --center \
          --width="400" \
          --window-icon=$imagem \
          --progress-text="$heditprogtxt" \
          --pulsate --auto-close
        salva_header_edit;
      fi
      if [ -z "$headeredit1" ];then cancel="1";fi
    done
  fi
}
#
##### gerando arquivo dat ####
#
geradat_header_edit(){
        echo "$hedoutin${intervalo:${#hedoutin}}" > $praiaheditdat
        echo "$hedouttg${intervalo:${#hedouttg}}" >> $praiaheditdat
        echo "$hedoutou${intervalo:${#hedoutou}}" >> $praiaheditdat
        z="$hedn $hedoutnm"; echo "$hedn $hedoutnm${intervalo:${#z}}" >> $praiaheditdat
        z="$hedf $hedoutfi"; echo "$hedf $hedoutfi${intervalo:${#z}}" >> $praiaheditdat
        z="$hede $hedoutep"; echo "$hede $hedoutep${intervalo:${#z}}" >> $praiaheditdat
        z="$heda $hedoutyr"; echo "$heda $hedoutyr${intervalo:${#z}}" >> $praiaheditdat
        z="$hedo $hedoutof"; echo "$hedo $hedoutof${intervalo:${#z}}" >> $praiaheditdat
        echo "$hedm${intervalo:${#hedm}}" >> $praiaheditdat
        echo "$hedd${intervalo:${#hedd}}" >> $praiaheditdat
        echo "$hedp${intervalo:${#hedp}}" >> $praiaheditdat
        z="$hedh $hedoutte"; echo "$hedh $hedoutte${intervalo:${#z}}" >> $praiaheditdat
        echo "$fimdedat" >> $praiaheditdat
}
#
#### Salvando variaveis ####
#
salva_header_edit(){
  echo "export hedoutin=\"$hedoutin\"" > $arqconfhe
  echo "export hedouttg=\"$hedouttg\"" >> $arqconfhe
  echo "export hedoutou=\"$hedoutou\"" >> $arqconfhe
  echo "export hedn=\"$hedn\"" >> $arqconfhe
  echo "export hedf=\"$hedf\"" >> $arqconfhe
  echo "export hede=\"$hede\"" >> $arqconfhe
  echo "export heda=\"$heda\"" >> $arqconfhe
  echo "export hedo=\"$hedo\"" >> $arqconfhe
  echo "export hedm=\"$hedm\"" >> $arqconfhe
  echo "export hedp=\"$hedp\"" >> $arqconfhe
  echo "export hedd=\"$hedd\"" >> $arqconfhe
  echo "export hedh=\"$hedh\"" >> $arqconfhe
  echo "export hedoutana=\"$hedoutana\"" >> $arqconfhe
  echo "export hedoutafl=\"$hedoutafl\"" >> $arqconfhe
  echo "export hedoutaex=\"$hedoutaex\"" >> $arqconfhe
  echo "export hedoutaan=\"$hedoutaan\"" >> $arqconfhe
  echo "export hedoutaof=\"$hedoutaof\"" >> $arqconfhe
  echo "export hedoutamu=\"$hedoutamu\"" >> $arqconfhe
  echo "export hedoutadi=\"$hedoutadi\"" >> $arqconfhe
  echo "export hedoutapr=\"$hedoutapr\"" >> $arqconfhe
  echo "export hedoutaep=\"$hedoutaep\"" >> $arqconfhe
  echo "export hedoutnm=\"$hedoutnm\"" >> $arqconfhe
  echo "export hedoutfi=\"$hedoutfi\"" >> $arqconfhe
  echo "export hedoutep=\"$hedoutep\"" >> $arqconfhe
  echo "export hedoutyr=\"$hedoutyr\"" >> $arqconfhe
  echo "export hedoutof=\"$hedoutof\"" >> $arqconfhe
  echo "export hedouttg=\"$hedouttg\"" >> $arqconfhe
  echo "export hedoutte=\"$hedoutte\"" >> $arqconfhe
}
#
########################## Funcao PRAIA Findscale ################################################################################
#
PRAIA_findscale(){
  arqconffs="$dirconf/findscale.sh"
#  . $arqconffs
  praiafindscale="PRAIA_findscale_20_04"
  praiafindscaledat="PRAIA_findscale_20_04.dat"
#
#### lendo informacoes do usuario ####
#
  if [ "$1" = "-e" ];then
    if [ -s "$praiafindscaledat" ];then
      roda_PRAIA $praiafindscale $praiafindscaledat $(pwd);
    else
      echo "dat file not found"
    fi
#  else
#    cancel="0"
#    until [ "$cancel" = "1" ];do
#      cancel="1"
#    done
  fi
#
##### gerando arquivo dat ####
#
}
#
########################## Funcao PRAIA Julian Date Gregorian ################################################################################
#
PRAIA_jd_gregorian(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA JPL Ephem Batch ################################################################################
#
PRAIA_JPL_ephem_batch(){
  arqconfeb="$dirconf/ephembatch.sh"
  . $arqconfeb
  praiajeb="PRAIA_JPL_ephem_batch_20_05"
  praiajebdat="PRAIA_JPL_ephem_batch_20_05.dat"
  if [ "$1" = "-e" ];then
    if [ -s "$praiajebdat" ];then
      roda_PRAIA $praiajeb $praiajebdat $(pwd);
    else
      echo "dat file not found"
    fi
  elif [ "$1" = "-r" ];then
    geradat_JPL_ephem_batch;
    roda_PRAIA $praiajeb $praiajebdat $(pwd);
  else
    cancel="0"
    until [ "$cancel" = "1" ];do
      jebvformatesc="$jebfohex!$jebfoxy!$jebfooff!$jebfompc!$jebfobig"
      yad --plug=200 \
      --tabnum=1 \
      --form \
      --field="$jebformat":CB "${jebvformatesc//$jebvformat/^$jebvformat}" \
      --field="$jebinfile" "$jebvfile" \
      --field="$jebobjnm" "$jebvnmobj" \
      --field="$jebobjcd" "$jebvcodeob" \
      --field="$jebckobj":CHK FALSE \
      --field="$jebjplcd" "$jebvcodelc" \
      --field="$jebckloc":CHK FALSE \
      --field="$jebemail" "$jebvmail" \
      --field="$jebdirectory":DIR "$jebvdir" \
      --field="$jebsubdirectory":CHK FALSE \
      --field="$jebdattxt":CHK FALSE \
      --field="$jebbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$jebreadme" &> .auxjeb1 &
      pid1=$!
      yad --plug=200 \
      --tabnum=2 \
      --form \
      --text="$jebcoortxt"  \
      --field="$jeblong" "$jebvlong" \
      --field="$jeblatit" "$jebvlat" \
      --field="$jebaltit" "$jebvalt" &> .auxjeb2 &
      pid2=$!
      yad --notebook \
      --center \
      --width=600 \
      --height=550 \
      --window-icon=$imagem \
      --title="$jebtitle" \
      --key=200 --tab="$jebtab1" --tab="$jebtab2"
      while kill -0 $pid1 || kill -0 $pid2; do
        /bin/true
      done
      jplephembatch1=$(cat .auxjeb1)
      jplephembatch2=$(cat .auxjeb2)
      rm .auxjeb1
      rm .auxjeb2
      if [ ! -z "$jplephembatch1" ];then   ## verifica se usuario nao clicou em cancel
        jebvformat=$(echo "$jplephembatch1" | cut -d"|" -f 1)
        jebvfile=$(echo "$jplephembatch1" | cut -d"|" -f 2)
        jebvnmobj=$(echo "$jplephembatch1" | cut -d"|" -f 3)
        jebvcodeob=$(echo "$jplephembatch1" | cut -d"|" -f 4)
        jebvckobj=$(echo "$jplephembatch1" | cut -d"|" -f 5)
        jebvcodelc=$(echo "$jplephembatch1" | cut -d"|" -f 6)
        jebvlocal=$(echo "$jplephembatch1" | cut -d"|" -f 7)
        jebvmail=$(echo "$jplephembatch1" | cut -d"|" -f 8)
        jebvdir=$(echo "$jplephembatch1" | cut -d"|" -f 9)
        jebvsubdir=$(echo "$jplephembatch1" | cut -d"|" -f 10)
        jebdat=$(echo "$jplephembatch1" | cut -d"|" -f 11)
        if [ "$jebvlocal" = "TRUE" ];then  ## usuario marcou pra definir coordenadas
          jebvcodelc="coord@399"
          jebvlong=$(echo "$jplephembatch2" | cut -d"|" -f 1)
          jebvlat=$(echo "$jplephembatch2" | cut -d"|" -f 2)
          jebvalt=$(echo "$jplephembatch2" | cut -d"|" -f 3)
        fi
        if [ "$jebvsubdir" = "TRUE" ];then  ## usuario quer buscar arquivo em subdiretorios
          busca="$jebvfile"
          lista_pastas $jebvdir;
          jebvdir2=$(cat $dirpadrao/.lista_pastas)
          rm $dirpadrao/.lista_pastas
        fi
        (if [ "$jebdat" = "FALSE" ];then
#          for d in $jebvdir2;do  ## escreve todos os arquivos em um só
#            cat $d/$jebvfile >> $jebvdir/$jebvfile
#          done
          case $jebvformat in
            "$jebfohex") jebfmt="1";;
            "$jebfoxy") jebfmt="2";;
            "$jebfooff") jebfmt="3";;
            "$jebfompc") jebfmt="4";;
            "$jebfobig") jebfmt="5";;
          esac
          if [ "$jebvckobj" = "TRUE" ];then jebvcodeob="$jebvnmobj";fi ## usuario nao sabe o codigo do obj. utiliza-se o nome
          cd $jebvdir
          geradat_JPL_ephem_batch;
        fi
        roda_PRAIA $praiajeb $praiajebdat $jebvdir;
        ) | yad --progress \
          --title "Ephem Batch" \
          --center \
          --width="400" \
          --window-icon=$imagem \
          --progress-text="$jebtxtpro" \
          --pulsate --auto-close
        salva_JPL_ephem_batch;
      fi
      if [ -z "$jplephembatch1" ];then cancel="1";fi
    done
  fi
}
#
#### gerando arquivo dat ####
#
geradat_JPL_ephem_batch(){
  echo "$jebfmt${intervalo:${#jebfmt}}" > $praiajebdat
  echo "$jebvfile${intervalo:${#jebvfile}}" >> $praiajebdat
  echo "$jebvnmobj${intervalo:${#jebvnmobj}}" >> $praiajebdat
  echo "$jebvcodelc${intervalo:${#jebvcodelc}}" >> $praiajebdat
  echo "$jebvlong${intervalo:${#jebvlong}}" >> $praiajebdat
  echo "$jebvlat${intervalo:${#jebvlat}}" >> $praiajebdat
  echo "$jebvalt${intervalo:${#jebvalt}}" >> $praiajebdat
  echo "$jebvcodeob${intervalo:${#jebvcodeob}}" >> $praiajebdat
  echo "$jebvmail${intervalo:${#jebvmail}}" >> $praiajebdat
  echo "$fimdedat" >> $praiajebdat
}
#
#### Salvando variaveis ####
#
salva_JPL_ephem_batch(){
  echo "export jebvformat=\"$jebvformat\"" > $arqconfeb
  echo "export jebfmt=\"$jebfmt\"" >> $arqconfeb
  echo "export jebvfile=\"$jebvfile\"" >> $arqconfeb
  echo "export jebvnmobj=\"$jebvnmobj\"" >> $arqconfeb
  echo "export jebvcodeob=\"$jebvcodeob\"" >> $arqconfeb
  echo "export jebvckobj=\"$jebvckobj\"" >> $arqconfeb
  echo "export jebvcodelc=\"$jebvcodelc\"" >> $arqconfeb
  echo "export jebvmail=\"$jebvmail\"" >> $arqconfeb
  echo "export jebvlong=\"$jebvlong\"" >> $arqconfeb
  echo "export jebvlat=\"$jebvlat\"" >> $arqconfeb
  echo "export jebvalt=\"$jebvalt\"" >> $arqconfeb
}
#
########################## Funcao NAIF ###########################################################################################
#
naif_ephem(){
  arqconfnaif="$dirconf/naif.sh"
  . $arqconfeb
  CONF=$( \
    yad --form \
    --center \
    --title="$naiftitle" \
    --width=400 \
    --height=200 \
    --window-icon=$imagem \
    --image=$imagem \
    --field="$confesidioma":CB "${confidioma//$idioma/^$idioma}" \
    --field="2MASS":DIR "$twomass" \
    --field="UCAC2":DIR "$ucac2" \
    --field="UCAC4":DIR "$ucac4" \
  )
}
#
########################## Funcao PRAIA Light Curve ################################################################################
#
PRAIA_light_curve_numerical_fit(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Occ Star Search ################################################################################
#
PRAIA_occ_star_search(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Proper Motions ################################################################################
#
PRAIA_proper_motions(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA RADEC DEC HEX ################################################################################
#
PRAIA_radec_dec_hex(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Statistics ################################################################################
#
PRAIA_statistics(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Targets JPL ################################################################################
#
PRAIA_targets_JPL(){
  arqconftjpl="$dirconf/targetsjpl.sh"
  . $arqconftjpl
  praiatgjpl="PRAIA_targets_JPL_20_01"
  praiatgjpldat="PRAIA_targets_JPL_20_01.dat"
  if [ "$1" = "-e" ];then  ####executa do terminal. ja possui dat
    if [ -s "$praiatgjpldat" ];then
      roda_PRAIA $praiatgjpl $praiatgjpldat $(pwd);
    else
      echo "dat file not found"
    fi
  elif [ "$1" = "-r" ];then  ####executa do terminal. usa ultimas configuracoes
    geradat_targets_JPL;
    roda_PRAIA $praiatgjpl $praiatgjpldat $(pwd);
  else                       ####abre interface grafica
    cancel="0"
    until [ "$cancel" = "1" ];do
      tgjplephemesc="NAIF!Horizons"
      targetsjpl=$( \
        yad --form \
        --center \
        --title="$tgjpltilte" \
        --width=400 \
        --height=300 \
        --window-icon=$imagem \
        --field="$tgjplnm" "$tgjplnome" \
        --field="$tgjpleph":CB "${tgjplephemesc//$tgjplephem/^$tgjplephem}" \
        --field="$tgjplineph" "$tgjplinput" \
        --field="$tgjploueph" "$tgjploutput" \
        --field="$tgjpldirectory":DIR "$tgjpldir" \
        --field="$tgjpldattxt":CHK FALSE \
        --field="$tgjplbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$tgjplreadme" \
      )
      if [ ! -z "$targetsjpl" ];then  ## verifica se usuario nao clicou em cancel
        tgjplnome=$(echo "$targetsjpl" | cut -d"|" -f 1)
        tgjplephem=$(echo "$targetsjpl" | cut -d"|" -f 2)
        tgjplinput=$(echo "$targetsjpl" | cut -d"|" -f 3)
        tgjploutput=$(echo "$targetsjpl" | cut -d"|" -f 4)
        tgjpldir=$(echo "$targetsjpl" | cut -d"|" -f 5)
        tgjpldat=$(echo "$targetsjpl" | cut -d"|" -f 6)
        case $tgjplephem in  ## ve qual efemeride escolhida
          "NAIF") tgjplfmt="1";;
          "Horizons") tgjplfmt="2";;
        esac
        cd $tgjpldir
        if [ "$tgjpldat" = "FALSE" ];then
          geradat_targets_JPL
        fi
        roda_PRAIA $praiatgjpl $praiatgjpldat $tgjpldir;
        salva_targets_JPL;
      fi
      if [ -z "$targetsjpl" ];then cancel="1";fi
    done
  fi
}
#
#### gerando arquivos dat ####
#
geradat_targets_JPL(){
  echo "$tgjplnome${intervalo:${#tgjplnome}}" > $praiatgjpldat
  echo "$tgjplfmt${intervalo:${#tgjplfmt}}" >> $praiatgjpldat
  echo "$tgjplinput${intervalo:${#tgjplinput}}" >> $praiatgjpldat
  echo "$tgjploutput${intervalo:${#tgjploutput}}" >> $praiatgjpldat
  echo "$fimdedat" >> $praiatgjpldat
}
#
#### Salvando variaveis ####
#
salva_targets_JPL(){
  echo "export tgjplnome=\"$tgjplnome\"" > $arqconftjpl
  echo "export tgjplephem=\"$tgjplephem\"" >> $arqconftjpl
  echo "export tgjplfmt=\"$tgjplfmt\"" >> $arqconftjpl
  echo "export tgjplinput=\"$tgjplinput\"" >> $arqconftjpl
  echo "export tgjploutput=\"$tgjploutput\"" >> $arqconftjpl
}
#
########################## Funcao PRAIA Target Seach ################################################################################
#
PRAIA_targets_search(){
  arqconfts="$dirconf/targetsearch.sh"
 . $arqconfts
  praiatgsrc="PRAIA_targets_search_20_02"
  praiatgsrcdat="PRAIA_targets_search_20_02.dat"
  if [ "$1" = "-e" ];then
    if [ -s "$praiatgsrcdat" ];then
      roda_PRAIA $praiatgsrc $praiatgsrcdat $(pwd);
    else
      echo "dat file not found"
    fi
  elif [ "$1" = "-r" ];then
    geradat_target_search;
    roda_PRAIA $praiatgsrc $praiatgsrcdat $(pwd);
  else
    cancel=0
    until [ "$cancel" = "1" ];do
      tgsrcvfmtesc="$tgsrcfoxy!$tgsrcfooff!$tgsrcfompc"
      targetsearch=$( \
        yad --form \
        --center \
        --title="$tgsrctitle" \
        --width=400 \
        --height=600 \
        --window-icon=$imagem \
        --field="$tgsrcformat":CB "${tgsrcvfmtesc//$tgsrcvfmt/^$tgsrcvfmt}" \
        --field="$tgsrcxy2mass":CHK "$tgsrcxy2m" \
        --field="$tgsrcxyucac2":CHK "$tgsrcxyu2" \
        --field="$tgsrcxyucac4":CHK "$tgsrcxyu4" \
        --field="$tgsrcxyuser":CHK "$tgsrcxyus" \
        --field="$tgsrcinxy" "$tgsrcxy" \
        --field="$tgsrcchktg":CHK "$tgsrctgck" \
        --field="$tgsrcintg" "$tgsrctg" \
        --field="$tgsrcout" "$tgsrcoutput" \
        --field="$tgsrcaerr" "$tgsrcerr" \
        --field="$tgsrcdirectory":DIR "$tgsrcdir" \
        --field="$tgsrccksubdir":CHK FALSE \
        --field="$tgsrcdattxt":CHK FALSE \
        --field="$tgsrcbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$tgsrcreadme" \
      )
      if [ ! -z "$targetsearch" ];then  ## verifica se usuario nao clicou em cancel
        tgsrcvfmt=$(echo "$targetsearch" | cut -d"|" -f 1)
        tgsrcxy2m=$(echo "$targetsearch" | cut -d"|" -f 2)
        tgsrcxyu2=$(echo "$targetsearch" | cut -d"|" -f 3)
        tgsrcxyu4=$(echo "$targetsearch" | cut -d"|" -f 4)
        tgsrcxyus=$(echo "$targetsearch" | cut -d"|" -f 5)
        tgsrcxy=$(echo "$targetsearch" | cut -d"|" -f 6)
        tgsrctgck=$(echo "$targetsearch" | cut -d"|" -f 7)
        tgsrctg=$(echo "$targetsearch" | cut -d"|" -f 8)
        tgsrcoutput=$(echo "$targetsearch" | cut -d"|" -f 9)
        tgsrcerr=$(echo "$targetsearch" | cut -d"|" -f 10)
        tgsrcdir=$(echo "$targetsearch" | cut -d"|" -f 11)
        tgsrcsubdir=$(echo "$targetsearch" | cut -d"|" -f 12)
        tgsrcdat=$(echo "$targetsearch" | cut -d"|" -f 13)
        case $tgsrcvfmt in  ## ve qual formato escolhida
          "$tgsrcfoxy") tgsrcfmt="1";;
          "$tgsrcfooff") tgsrcfmt="2";;
          "$tgsrcfompc") tgsrcfmt="3";;
        esac
        (cd $tgsrcdir
        if [ "$tgsrcdat" = "FALSE" ];then
          rm $praiatgsrcdat
          if [ "$tgsrctgck" = "TRUE" ];then
            alvos=$(ls target_*)
          else
            alvos="$tgsrctg"
          fi
          for alvo in $alvos;do
            geradat_target_search
          done
        fi
        roda_PRAIA $praiatgsrc $praiatgsrcdat $tgsrcdir;
        ) | yad --progress \
          --title "Target Search" \
          --center \
          --width="400" \
          --window-icon=$imagem \
          --progress-text="$tgsrcpro" \
          --pulsate --auto-close
        salva_target_search;
      fi
      if [ -z "$targetsearch" ];then cancel="1";fi
    done
  fi
}
#
##### gerando arquivo dat Target Seach####
#
geradat_target_search(){
  if [ "$tgsrcxy2m" = "TRUE" ];then
    if [ -f "list.2mass.xy" ];then
      ls *.2mass.red.xy > list.2mass.xy
    fi
    if [ "$tgsrctgck" = "TRUE" ];then
      tgoutput="2mass_${alvo#target_}"
    else
      tgoutput="2mass_$tgsrcoutput"
    fi
    echo "$tgsrcfmt${intervalo:${#tgsrcfmt}}" >> $praiatgsrcdat
    echo "list.2mass.xy                                     |" >> $praiatgsrcdat
    echo "$alvo${intervalo:${#alvo}}" >> $praiatgsrcdat
    echo "$tgoutput${intervalo:${#tgoutput}}" >> $praiatgsrcdat
    echo "$tgsrcerr${intervalo:${#tgsrcerr}}" >> $praiatgsrcdat
    echo "$fimdedat" >> $praiatgsrcdat
  fi
  if [ "$tgsrcxyu2" = "TRUE" ];then
    if [ -f "list.ucac2.xy" ];then
      ls *.ucac2.red.xy > list.ucac2.xy
    fi
    if [ "$tgsrctgck" = "TRUE" ];then
      tgoutput="ucac2_${alvo#target_}"
    else
      tgoutput="ucac2_$tgsrcoutput"
    fi
    echo "$tgsrcfmt${intervalo:${#tgsrcfmt}}" >> $praiatgsrcdat
    echo "list.ucac2.xy                                     |" >> $praiatgsrcdat
    echo "$alvo${intervalo:${#alvo}}" >> $praiatgsrcdat
    echo "$tgoutput${intervalo:${#tgoutput}}" >> $praiatgsrcdat
    echo "$tgsrcerr${intervalo:${#tgsrcerr}}" >> $praiatgsrcdat
    echo "$fimdedat" >> $praiatgsrcdat
  fi
  if [ "$tgsrcxyu4" = "TRUE" ];then
    if [ -f "list.ucac4.xy" ];then
      ls *.ucac4.red.xy > list.ucac4.xy
    fi
    if [ "$tgsrctgck" = "TRUE" ];then
      tgoutput="ucac4_${alvo#target_}"
    else
      tgoutput="ucac4_$tgsrcoutput"
    fi
    echo "$tgsrcfmt${intervalo:${#tgsrcfmt}}" >> $praiatgsrcdat
    echo "list.ucac4.xy                                     |" >> $praiatgsrcdat
    echo "$alvo${intervalo:${#alvo}}" >> $praiatgsrcdat
    echo "$tgoutput${intervalo:${#tgoutput}}" >> $praiatgsrcdat
    echo "$tgsrcerr${intervalo:${#tgsrcerr}}" >> $praiatgsrcdat
    echo "$fimdedat" >> $praiatgsrcdat
  fi
  if [ "$tgsrcxyus" = "TRUE" ];then
    if [ -f "list.wfi.xy" ];then
      ls *.wfi.red.xy > list.wfi.xy
    fi
    if [ "$tgsrctgck" = "TRUE" ];then
      tgoutput="wfi_${alvo#target_}"
    else
      tgoutput="wfi_$tgsrcoutput"
    fi
    echo "$tgsrcfmt${intervalo:${#tgsrcfmt}}" >> $praiatgsrcdat
    echo "list.wfi.xy                                       |" >> $praiatgsrcdat
    echo "$alvo${intervalo:${#alvo}}" >> $praiatgsrcdat
    echo "$tgoutput${intervalo:${#tgoutput}}" >> $praiatgsrcdat
    echo "$tgsrcerr${intervalo:${#tgsrcerr}}" >> $praiatgsrcdat
    echo "$fimdedat" >> $praiatgsrcdat
  fi
  if [ "$tgsrcxy" != "" ];then
    if [ "$tgsrctgck" = "TRUE" ];then
      tgoutput="search_results_${alvo#target_}"
    else
      tgoutput="search_results_$tgsrcoutput"
    fi
    echo "$tgsrcfmt${intervalo:${#tgsrcfmt}}" >> $praiatgsrcdat
    echo "$tgsrcxy${intervalo:${#tgsrcxy}}" >> $praiatgsrcdat
    echo "$alvo${intervalo:${#alvo}}" >> $praiatgsrcdat
    echo "$tgoutput${intervalo:${#tgoutput}}" >> $praiatgsrcdat
    echo "$tgsrcerr${intervalo:${#tgsrcerr}}" >> $praiatgsrcdat
    echo "$fimdedat" >> $praiatgsrcdat
  fi
}
#
#### Salvando variaveis ####
#
salva_target_search(){
  echo "export tgsrcfmt=\"$tgsrcfmt\"" > $arqconfts
  echo "export tgsrcvfmt=\"$tgsrcvfmt\"" >> $arqconfts
  echo "export tgsrcxy2m=\"$tgsrcxy2m\"" >> $arqconfts
  echo "export tgsrcxyu2=\"$tgsrcxyu2\"" >> $arqconfts
  echo "export tgsrcxyu4=\"$tgsrcxyu4\"" >> $arqconfts
  echo "export tgsrcxyus=\"$tgsrcxyus\"" >> $arqconfts
  echo "export tgsrcxy=\"$tgsrcxy\"" >> $arqconfts
  echo "export tgsrctgck=\"$tgsrctgck\"" >> $arqconfts
  echo "export tgsrctg=\"$tgsrctg\"" >> $arqconfts
  echo "export tgsrcoutput=\"$tgsrcoutput\"" >> $arqconfts
  echo "export tgsrcerr=\"$tgsrcerr\"" >> $arqconfts
}
#
########################## Funcao PRAIA Trajectory ################################################################################
#
PRAIA_trajectory(){
x=1
#
#### lendo informacoes do usuario ####
#
#
##### gerando arquivo dat ####
#
#
#### Rodando PRAIA ####
#
}
#
########################## Funcao PRAIA Trajectory ################################################################################
#
PRAIA_uniform(){
  arqconfuni="$dirconf/uniform.sh"
  . $arqconfuni
  praiauniform="PRAIA_uniform_astrometry_20_02"
  praiauniformdat="PRAIA_uniform_astrometry_20_02.dat"
  cancel="0"
#
#### lendo informacoes do usuario ####
#
  until [ "$cancel" = "1" ];do
    if [ "$1" != "-e" ];then
      unicatesc="2MASS!UCAC2!UCAC4!$unifuscat"
      uniform=$( \
      yad --form \
        --center \
        --title="$uniftitle" \
        --width=400 \
        --height=400 \
        --window-icon=$imagem \
        --field="$unifchcat":CB "${unicatesc//$unicat/^$unicat}" \
        --field="$unifperfra":SCL $uniper \
        --field="$unifcutoff" "$unicut" \
        --field="$unifpolin":NUM $unipolin!0..3!1 \
        --field="$unifrad3":CHK "$unir3" \
        --field="$unifrad5":CHK "$unir5" \
        --field="$unifdir":DIR "$unidir" \
        --field="$unifckdat":CHK FALSE \
        --field="$unifbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/manual_unif_$manidioma" \
      )
      if [ -z "$uniform" ];then cancel="1";fi
      if [ ! -z "$uniform" ];then
        unicat=$(echo "$uniform" | cut -d"|" -f 1)
        uniper=$(echo "$uniform" | cut -d"|" -f 2)
        unicut=$(echo "$uniform" | cut -d"|" -f 3)
        unipolin=$(echo "$uniform" | cut -d"|" -f 4)
        unipolin=$(echo $unipolin | cut -d . -f 1 )
        unir3=$(echo "$uniform" | cut -d"|" -f 5)
        unir5=$(echo "$uniform" | cut -d"|" -f 6)
        unidir=$(echo "$uniform" | cut -d"|" -f 7)
        unidat=$(echo "$uniform" | cut -d"|" -f 8)
        case $unicat in  ## ve qual efemeride escolhida
          "2MASS") unisig="2mass";;
          "UCAC2") unisig="ucac2";;
          "UCAC4") unisig="ucac4";;
          "$unifuscat") unisig="wfi";;
        esac
        if [ "$unir3" = "TRUE" ];then unir3p="3"; else unir3p="0"; fi
        if [ "$unir5" = "TRUE" ];then unir5p="5"; else unir5p="0"; fi
        cd $unidir
        unixylist="list.$unisig.xy"
        unixylist2="list2.$unisig.xy"
        ref="0"
        until [ "$ref" = "1" ];do
          unixyext="$unisig.uni.xy"
          rm *.$unixyext list.uni.$unisig.xy
          if [ ! -f "$unixylist" ];then
            ls *.$unisig.red.xy > $unixylist
          fi
          if [ ! -s "$unixylist2" ];then
            cat $unixylist > $unixylist2
          fi
          if [ ! -f "$unixylist2" ];then
            while read line;do
              echo "TRUE $line" >> .auxuni2
            done < $unixylist
          else
            rm .auxuni2
            while read line;do
              uniaux=$(grep "$line" $unixylist2)
              if [ ! -z $uniaux ];then
                echo "TRUE $line" >> .auxuni2
              else
                echo "FALSE $line" >> .auxuni2
              fi
            done < $unixylist
          fi
          unilist=$(cat .auxuni2)
#
##### gerando arquivo dat ####
#
          (if [ "$unidat" = "FALSE" ];then
            echo "$unixylist2${intervalo:${#unixylist2}}" > $praiauniformdat
            echo "$unixyext${intervalo:${#unixyext}}" >> $praiauniformdat
            echo "$uniper${intervalo:${#uniper}}" >> $praiauniformdat
            echo "$unicut${intervalo:${#unicut}}" >> $praiauniformdat
            echo "$unipolin${intervalo:${#unipolin}}" >> $praiauniformdat
            echo "$unir3p${intervalo:${#unir3p}}" >> $praiauniformdat
            echo "$unir5p${intervalo:${#unir5p}}" >> $praiauniformdat
            echo "$fimdedat" >> $praiauniformdat
          fi
#
#### Rodando PRAIA ####
#
          cp $dirfontes/$praiauniform $unidir
          ./$praiauniform < $praiauniformdat > .auxuni
          rm $praiauniform) | yad --progress \
            --title "Uniform Astrometry" \
            --center \
            --width="400" \
            --window-icon=$imagem \
            --progress-text="$unifpro" \
            --pulsate --auto-close
          ls *.$unixyext > list.uni.$unisig.xy
#
#### Salvando variaveis ####
#
          echo "export uniper=\"$uniper\"" > $arqconfuni
          echo "export unicat=\"$unicat\"" >> $arqconfuni
          echo "export unicut=\"$unicut\"" >> $arqconfuni
          echo "export unipolin=\"$unipolin\"" >> $arqconfuni
          echo "export unir3=\"$unir3\"" >> $arqconfuni
          echo "export unir5=\"$unir5\"" >> $arqconfuni
#
#### Vendo Resultados ####
#
          yad --plug=400 \
            --tabnum=1 \
            --text-info \
            --wrap \
            --filename=".auxuni" \
            --fontname="LMMonoSlant10" &> .auxuni4 &
          pid1=$!
          yad --plug=400 \
            --tabnum=2 \
            --list \
            --checklist \
            --separator="" \
            --print-column=2 \
            --no-headers \
            --column="M" \
            --column="Arquivo XY" \
            $unilist > .auxuni3 &
          pid2=$!
          yad --notebook \
            --center \
            --width=650 \
            --height=600 \
            --window-icon=$imagem \
            --title="$uniftitle" \
            --button="$unifbtref":0 --button="$unifbtfim":1 \
            --key=400 --tab="$uniftab1" --tab="$uniftab2"
          res=$?
#          while kill -0 $pid1 || kill -0 $pid2; do
#            /bin/true
#          done
#
#### Verificando re-reduçao ####
#
          case $res in
            "1") ref="1";;
            "0") ref="0"; rm *.$unixyext; unidat="TRUE";;
          esac
          rm .auxuni4
          if [ "$ref" = "0" ];then
            uniflist=$(cat .auxuni3)
            rm $unixylist2
            for i in $uniflist;do
              echo "$i" >> $unixylist2
            done
          fi
          rm .auxuni3
        done
      fi
#
#### Executando do terminal ####
#
    else
      if [ -s "$praiauniformdat" ];then
        cancel="1"
        roda_PRAIA $praiauniform $praiauniformdat $(pwd);
      else
        echo "dat file not found"
        cancel="1"
      fi
    fi
  done
}
#
#### Rodando PRAIA ####
#
roda_PRAIA(){
  cp $dirfontes/$1 $3
  ./$1 < $2
  rm $1
}
#
########################## Funcao Adicionar Script ################################################################################
#
add_script(){
  adiciona=$(yad --form \
    --title "$addscripttitle" \
    --height=325 \
    --window-icon=$imagem \
    --button="$sair":1 --button="$addscriptbtadd":2 \
    --text="$addscripttext" \
    --field="$addscriptfile":FL "" \
    --field="$addscriptname" "" \
    --field="$addscriptcom" "./script ou bash script" \
    --field="$addscriptinf":TXT "$addscriptrdm" \
    --field="$scriptbtnreadme":BTN "yad --text-info --center --width=800 --height=400 --window-icon=$imagem --fontname=LMMonoSlant10 --filename=$dirmanual/$scriptreadme" \
  )
  addscfl=$(echo "$adiciona" | cut -d"|" -f 1)
  addscnm=$(echo "$adiciona" | cut -d"|" -f 2)
  addsccm=$(echo "$adiciona" | cut -d"|" -f 3)
  addscinfo=$(echo "$adiciona" | cut -d"|" -f 4)
  if [ ! -z "$adiciona" ];then
    rm $manualscript
    cp $addscfl $dirscript
    addfile=$(echo ${addscfl##*/})
    echo -e "$addscnm\t$addfile\t$addsccm" >> $listascript
    echo -e "$addscinfo" > $dirscript/"$addscnm".nfo
    var=$(sort -n --key=1 $listascript)
    echo "$var" > $listascript
    while read line;do
      addnm=$(echo "$line" | cut -f 1)
      echo -e "\n--- $addnm ---------------------------------" >> $manualscript
      cat $dirscript/"$addnm".nfo >> $manualscript
    done < $listascript
  fi
}
#
########################## Funcao Rodar Proprio Programa ################################################################################
#
Rodar_script(){
#
#### lendo informacoes do usuario ####
#
  arqconfscr="$dirconf/script.sh"
  . $arqconfscr
  listascript="$dirscript/.lista_scripts"
  manualscript="$dirmanual/manual_script"
  cancel="0"
  until [ "$cancel" = "1" ];do
    run="0"
    while read line;do
      scnm=$(echo "$line" | cut -f 1)
      echo "$scnm!" >> .auxadd
    done < $listascript
    scriptoption=$(cat .auxadd)
    scriptoption=$(echo ${scriptoption%!*})
    rm .auxadd
    script=$(yad --form \
      --title "$scripttitle" \
      --center
      --height=200 \
      --window-icon=$imagem \
      --button="$sair":1 --button="$scriptbtadd":2 --button="$scriptrodar":0 \
      --text="$scripttext" \
      --field="$scriptchtxt":CB "${scriptoption//$scriptch/^$scriptch}" \
      --field="$scriptdirectory":DIR "$scriptdir" \
      --field="$scriptprobar":CHK FALSE \
      --field="$scriptbtnreadme":BTN "yad --text-info --center --width=600 --height=400 --window-icon=$imagem --filename=$manualscript" \
    )
    res=$?
    scriptch=$(echo "$script" | cut -d"|" -f 1)
    read scriptch <<< $scriptch
    scriptdir=$(echo "$script" | cut -d"|" -f 2)
    scriptpro=$(echo "$script" | cut -d"|" -f 3)
    case $res in
      0) cancel="0"; run="1";;
      1) cancel="1";;
      2) cancel="0"; add_script;;
    esac
#
#### Rodando Script ####
#
    if [ "$run" = "1" ];then
      scarq=$(grep "$scriptch" $listascript | cut -f 2)
      sccomand=$(grep "$scriptch" $listascript | cut -f 3)
      cp $dirscript/"$scarq" $scriptdir
      cd $scriptdir
      if [ "$scriptpro" = "FALSE" ];then
        $sccomand
      else
        ($sccomand) | yad --progress \
      --title "$scripttitle" \
      --center \
      --width="400" \
      --window-icon=$imagem \
      --progress-text="$scriptruntitle $scriptch" \
      --pulsate --auto-close --auto-kill
      fi
      echo "export scriptch=\"$scriptch\"" > $arqconfscr
    fi
  done
}
#
########################## Funcao Configuration ################################################################################
#
configuration(){
  confidioma="English!Français!Português"
  CONF=$( \
    yad --form \
    --center \
    --title="$conftitle" \
    --width=400 \
    --height=200 \
    --window-icon=$imagem \
    --image=$imagem \
    --field="$confesidioma":CB "${confidioma//$idioma/^$idioma}" \
    --field="2MASS":DIR "$twomass" \
    --field="UCAC2":DIR "$ucac2" \
    --field="UCAC4":DIR "$ucac4" \
  )
  if [ ! -z "$CONF" ];then
    idioma=$(echo "$CONF" | cut -d"|" -f 1)
    twomass=$(echo "$CONF" | cut -d"|" -f 2)
    ucac2=$(echo "$CONF" | cut -d"|" -f 3)
    ucac4=$(echo "$CONF" | cut -d"|" -f 4)
    case $idioma in
      "Português") manidioma="pt";;
      "English") manidioma="en";;
      "Français") manidioma="fr";;
    esac
    . $dirlanguage/language_$manidioma.sh;
    echo "export manidioma=\"$manidioma\"" > $arqconf
    echo "export idioma=\"$idioma\"" >> $arqconf
    echo "export twomass=\"$twomass\"" >> $arqconf
    echo "export ucac2=\"$ucac2\"" >> $arqconf
    echo "export ucac4=\"$ucac4\"" >> $arqconf
    echo "export usercat=\"$usercat\"" >> $arqconf
  fi
}
#export -f configuration
#
########################## Funcao Busca Subdiretorios ################################################################################
#
lista_pastas(){
  cd $1
  if [ "$(ls $busca)" ];then
    pwd >> $dirpadrao/.lista_pastas
  fi
  i=$(ls -d */)
  for j in $i;do
    lista_pastas $j
  done
  cd ..
}
#
########################## Corpo do Programa ##########################################
#
. $arqconf
. $dirlanguage/language_$manidioma.sh
case $1 in
  "") inicio;;
  "--astrometry") PRAIA_astrometry $2;;
  "--big-table") PRAIA_big_table $2;;
#  "--coronografy") PRAIA_coronografy $2; echo "function not yet installed";;
#  "--global-red") PRAIA_global_reduction $2; echo "function not yet installed";;
  "--hextraction") PRAIA_header_extraction_30 $2;;
  "--hedit") PRAIA_header_edit $2;;
  "--findscale") PRAIA_findscale $2;;
#  "--jd-gregorian") PRAIA_jd_gregorian $2 echo "function not yet installed";;
  "--eph-batch") PRAIA_JPL_ephem_batch $2;;
#  "--light-curve") PRAIA_light_curve_numerical_fit $2; echo "function not yet installed";;
#  "--occ-star") PRAIA_occ_star_search $2; echo "function not yet installed";;
#  "--proper-motion") PRAIA_proper_motions $2; echo "function not yet installed";;
#  "--radec") PRAIA_radec_dec_hex $2; echo "function not yet installed";;
#  "--statistics") PRAIA_statistics $2; echo "function not yet installed";;
  "--tgt-jpl") PRAIA_targets_JPL $2;;
  "--tgt-search") PRAIA_targets_search $2;;
#  "--trajectory") PRAIA_trajectory $2; echo "function not yet installed";;
  "--uniform") PRAIA_uniform $2;;
  "--script") Rodar_script $2;;
  "--config") configuration;;
#  "--options") cat $dirmanual/options;;
  "--version") echo "$praiaversion";;
esac
