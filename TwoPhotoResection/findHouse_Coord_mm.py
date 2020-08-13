import numpy as np
import pickle

with open('Affine.pickle', 'rb') as f: # from older Lab Affine para
    Affine =  pickle.load(f) 

house_coord_mm = {}
house_coord = { '621' : ( 7609.95, 9104.50),
            '622' : ( 7675.58, 9088.25),
            '623' : ( 7688.45, 9141.88),
            '624' : ( 7623.08, 9158.25),
            '625' : ( 7697.83, 9162.38),
            '631' : ( 1518.41, 9281.66),
            '632' : ( 1585.17, 9265.79),
            '633' : ( 1598.05, 9317.54),
            '634' : ( 1529.92, 9335.41),
            '635' : ( 1608.05, 9339.29)
            }

for pnt,coord in house_coord.items():
    xy = np.array( [ [coord[0]] , 
                     [coord[1]] , 
                     [1 ] 
                     ] )
    for i in range(4):
        photo_coord = np.matmul( Affine['AffMat62'] , xy )

    else:
        photo_coord = np.matmul( Affine['AffMat63'] , xy )
        
    house_coord_mm[pnt] = ( photo_coord[0][0], photo_coord[1][0] )

    

print('-----end------')
with open('house_coord_mm.pickle', 'wb') as f: # from older Lab Affine para
    pickle.dump( house_coord_mm, f, pickle.HIGHEST_PROTOCOL )
1529.92, 9335.4
