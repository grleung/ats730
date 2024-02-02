#!/bin/bash

for i in {0..70}; do
    rvp=$(echo "scale=4;0.007 + ($i*0.0001)" |bc)

      # Generate the NAMELIST file
      { echo "&output"
        echo " base_out = .False."
        echo " base_outpath = ''"
        echo " parcel_out = .True."
        echo " parcel_outpath = 'output-rv/rvp"$i"_parceltraj.txt'"
        echo "/"
        echo "&parcel"
        echo " rvpsurf = $rvp"
        echo "/"
        echo "&grid"
        echo " nz = 40"
        echo " dz0= 700"
        echo "/"
        
      } > Namelist

      # Now run your program
      ./MAC.exe

done

