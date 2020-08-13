import numpy as np


AffMat62 = np.array([ [  0.015011580041652846,  3.5451578812471363e-07,  -113.49223266903768],
                    [  -1.1070871431591042e-07,  -0.015011955912427185,  113.50497970297891]  ] )

AffMat63 = np.array([ [  0.0150119124,  2.21365058e-08,   -113.494109],
                    [  -4.43075446e-07,  -1.50119559e-02,  1.13505616e+02]  ] )


ImageCoord = {
    '1001' : ( 15088, 7561 )  ,
    '1002' : ( 34, 7561 ) ,
    '1003' : ( 7561, 32 )
       }

for pnt,coord in ImageCoord.items():
    xy = np.array( [ [coord[0]] , 
                     [coord[1]] , 
                     [1 ] 
                     ] )
    photo_coord = np.matmul( AffMat62 , xy )
    print( '{} :\n {}'.format( pnt,photo_coord) )

print('-----end------')
