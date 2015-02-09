sc mkdir SAC_rot
do file wild  *.BH1.M.SAC
 setbb north $file
 setbb east '( CHANGE 'BH1' 'BH2' %north )'
 setbb vert '( CHANGE 'BH1' 'BHZ' %north )'
 setbb cmpr '( CHANGE 'BH1.M.SAC' 'r' %north )'
 setbb cmpt '( CHANGE 'BH1.M.SAC' 't' %north )'
 setbb cmpz '( CHANGE 'BH1.M.SAC' 'z' %north )'
 cut 0 149
 r  %north %east
 rotate to GCP
 w SAC_rot/%cmpr SAC_rot/%cmpt
 r %vert
 w SAC_rot/%cmpz
enddo

