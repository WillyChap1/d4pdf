#!/bin/bash

levs='925 850 700 600 500 400 300 1000'

for kk in $levs
do

    cd /Users/wchapman/SIO/research/100Amip/d4pdf/Qlev
    echo $kk

    for (( jj=1951; jj<=2011; jj++ ))
    do
        echo another year!                                                                                                                                                                  
	for i in $( ls Q"$kk"_*.nc ); do
	    FILEBASENAME=$(echo $i | cut -d. -f1)
	    SUBSTRING=$(echo $FILEBASENAME| cut -d'_' -f 2)

	    
	    
    #grab time domain 
	    ncks -d time,"$jj-01-01 00:00","$jj-12-31 18:00" "Q$kk"_$SUBSTRING.nc "Q$kk"____"$jj"0101_"$jj"1231_$SUBSTRING.nc
   
	done

#concatenate the files                                                                                                                            
	ncecat -O -u ENS Q"$kk"____"$jj"0101_"$jj"1231_??.nc Q"$kk"____"$jj"0101_"$jj"1231.nc

#remove the files                                                                                                                                
	rm -f Q"$kk"____"$jj"0101_"$jj"1231_*.nc
     

	((j++)) #increment year. 
    done
done 

cd /Users/wchapman/Git/d4pdf

