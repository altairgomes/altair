while read linha
do
    echo "${linha:68:16} JD" >> julian_date_OHP
done < header_extraction.dat
