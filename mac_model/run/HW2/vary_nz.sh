#!/bin/bash

for dz in 700 350 175 1400 2800; do

      # Generate the NAMELIST file
      { echo "&output"
        echo " base_out = .False."
        echo " base_outpath = ''"
        echo " parcel_out = .True."
        echo " parcel_outpath = 'output-nz/dz"$dz"_parceltraj.txt'"
        echo "/"
        echo "&parcel"
        echo " rvpsurf = 11.5e-3"
        echo "/"
        echo "&grid"
        echo " nz = $(((27300/dz) + 1))"
        echo " dz0= $dz"
        echo "/"
        
      } > Namelist

      # Now run your program
      ./MAC.exe

done
