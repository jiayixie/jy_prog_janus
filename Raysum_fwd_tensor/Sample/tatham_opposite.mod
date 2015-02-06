# Layers are listed from top to bottom. The bottom layer is
# assumed to be a half-space. Interface strike and dip apply
# to the upper interface of the layer.
#
# Format:
#     Column   Contents
#        1    Thickness (m)
#        2    Density (kg/m^3)
#        3    Average P-wave velocity (m/s)
#        4    Average S-wave velocity (m/s)
#        5    Isotropic-layer flag (1:isotropic, 0:anisotropic)
#        6    %P anisotropy
#        7    %S anisotropy (if 5 and 6 are zero, isotropic layer)
#        8    Trend of fast axis (degrees)
#        9    Plunge of fast axis (degrees)
#       10    Interface strike (degrees)
#       11    Interface dip (degrees)  
#     Note that the percentages of anisotropy are peak-to-peak
#     (the expressions used are from Farra et al. (1991))
#
# seis-tensor is a hack and needs 3 aniso layers, no more.
# Layers: crust, anisotropic wedge, isotropic half-space.
#thick rho  alph beta iso %P  %S  tr   pl  st di
20000  2800 6500 3810  1   0   0   0     0   0  0
 5000  2800 6500 3810  0  -6  -6  240   20   0   0
 5000  2800 6500 3810  0  -6  -6  160   55  70  55
 5000  2800 6500 3810  0  -6  -6  240   20  70  55
