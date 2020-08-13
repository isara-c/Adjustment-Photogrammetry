import numpy as np
import pickle

AffMat62 = np.array([ [  1.5011859226653e-02, -2.8365288601049e-7, -113.50077766574],
                    [  -1.8716778648295e-7, -1.5012983008061e-02, 113.50467204681]  ] )

AffMat63 = np.array([ [  1.5011224390783e-02, 1.4843241403917e-07, -113.49751765623],
                    [  -2.9018885425355e-07, -1.501169666325e-02, 0.11349384986359e+03]  ] )
#----Tle
TiePoints62 = { 't1' : ( 8494.16, 11028.81),
            't2' : ( 10081.89, 9717.18),
            't3' : ( 11374.63, 10243.39),
            't4' : ( 10324.39, 7792.17),
            't5' : ( 9578.75, 7417.41),
            't6' : ( 9351.51, 5321.74)
            }

TiePoints63 = { 't1' : ( 2402.89, 11222.63),
            't2' : ( 3927.12, 9919.85),
            't3' : ( 5176.54, 10461.05),
            't4' : ( 4328.46, 7983.93),
            't5' : ( 3600.49, 7600.35),
            't6' : ( 3130.99, 5502.19)
            }


for pnt,coord in TiePoints62.items():
    xy = np.array( [ coord[0] , 
                     coord[1] , 
                     1  ] )
    
    photo_coord = np.matmul( AffMat62 , xy )
    print( '{} : {:.3f} {:.3f}'.format(pnt, photo_coord[0], photo_coord[1] ))
    
print('\n')
for pnt,coord in TiePoints63.items():
    xy = np.array( [ coord[0] , 
                     coord[1] , 
                     1 ] )
    
    photo_coord = np.matmul( AffMat63 , xy )
    print( '{} : {:.3f} {:.3f}'.format(pnt, photo_coord[0], photo_coord[1] ))
