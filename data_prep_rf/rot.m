$KEYS dirot te
do file wild *.BHZ.M.SAC
  setbb vert $file
  setbb east '(CHANGE 'BHZ' 'BHE' %vert )'
  setbb north '( CHANGE 'BHZ' 'BHN' %vert )'
  setbb cmpr '( CHANGE 'BHZ.M.SAC' 'r' %vert )'
  setbb cmpt '( CHANGE 'BHZ.M.SAC' 't' %vert )'
  setbb cmpz '( CHANGE 'BHZ.M.SAC' 'z' %vert )'
  cut 0 $te
  r %north %east
  rmean
  rtr
  rotate to GCP
  w $dirot$/%cmpr $dirot$/%cmpt
  r %vert
  rmean
  rtr
  w $dirot$/%cmpz
enddo
