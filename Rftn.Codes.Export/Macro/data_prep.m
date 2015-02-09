sc mkdir Trash
sc mkdir GoodOnes
qdp off
ygrid on
do file wild *BHZ
 setbb vert $file
 setbb east '( CHANGE 'BHZ' 'BHE' %vert )'
 setbb north '( CHANGE 'BHZ' 'BHN' %vert )'
 r %vert %east %north
 synch
 w over
 r %vert
 rmean
 rtr
 ppk
 setbb t0 &1,t0
 r %vert %east %north
 ch t0 %t0
 w over
 cut t0 -60 t0 +90
 r %vert %east %north
 rmean
 rtr
 taper w 0.2
 div 6.307450 6.465442 6.478608
 w over
 p1
 setbb resp (REPLY "Enter t to trash the file")
 if %resp eq "t" then
   sc mv %vert Trash
   sc mv %east Trash
   sc mv %north Trash
 else
   sc mv %vert GoodOnes
   sc mv %east GoodOnes
   sc mv %north GoodOnes
 
 endif
 cut off
enddo
