from math import sqrt, cos, sin, atan2, degrees
import numpy as np

p = 6.8 # micron from doc
Distance = [int(i+1) for i in range(31)]
Distortion  = [-0.022,-0.175,-0.591,-1.397,-2.72,-4.68,-7.396,-10.979,-15.534,-21.158,\
-27.941,-35.963,-45.295,-55.999,-68.124,-81.708,-96.779,-113.352,-131.427,-150.996,-172.036,\
-194.509,-218.368,-243.552,-269.987,-297.588,-326.258,-355.892,-386.37,-417.569,-449.354]


def FindPhotogCoord( ij):
    i, j = ij
    i_PPAC = p / 1000 * ( i - 3616.46 )
    j_PPAC = p / 1000 * ( j - 2723.55 )
    
#--------- distortion-------
    r = sqrt( i_PPAC**2 + j_PPAC**2 )
    theta = atan2( j_PPAC, i_PPAC )
    print('theta = {:.5} degs. '.format( degrees(theta) ))
    delta_r = np.interp( r, Distance, Distortion )
    r = r - delta_r/1000
    print('r = {:.4} mm. '.format( r )  )

#----------find Xp Yp-------
    Xp = r * cos(theta)
    Yp = r * sin(theta)
    return ( Xp, Yp )

XY = { 'BottomLeft': ( 5516.18, 2689.11 ),
       'TopLeft': ( 5516, 1848 ),
       'TopRight': ( 5902, 1858 ),
       'BottomRight': ( 5904, 2700 )}

for key, value in XY.items() :
    print(key , '|  X : {} pixel , Y : {} pixel'.format( *value) )
    Xp, Yp = FindPhotogCoord( value )
    print(('Xp = {:.4} mm. | Yp = {:.4} mm.').format(Xp, Yp), '\n')

