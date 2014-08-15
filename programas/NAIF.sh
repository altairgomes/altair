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
