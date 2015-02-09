sc mkdir SAC_rot
do file wild  *.BHN.M.SAC
 setbb north $file
 setbb east '( CHANGE 'BHN' 'BHE' %north )'
 setbb vert '( CHANGE 'BHN' 'BHZ' %north )'
 setbb cmpr '( CHANGE 'BHN.M.SAC' 'r' %north )'
 setbb cmpt '( CHANGE 'BHN.M.SAC' 't' %north )'
 setbb cmpz '( CHANGE 'BHN.M.SAC' 'z' %north )'
 cut 0 140
 r  %north %east
 decimate 2
 rotate to GCP
 w  SAC_rot/%cmpr SAC_rot/%cmpt
 r  %vert
 decimate 2
 w  SAC_rot/%cmpz
enddo

