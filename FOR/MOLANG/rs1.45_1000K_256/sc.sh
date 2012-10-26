#!/bin/tcsh

../exe.x < in

foreach file ( angle*.hist )


  ../png/matrix.py $file output.rix
  ../png/png.py output.rix output.png
  mv output.png "$file".png

end
