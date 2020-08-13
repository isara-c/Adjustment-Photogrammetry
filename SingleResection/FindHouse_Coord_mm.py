import numpy as np
import pickle

with open('Affine.pickle', 'rb') as f: # from older Lab Affine para
    Affine =  pickle.load(f) 

house_coord = { '621' : ( 8101, 7515),
            '622' : ( 8118, 7555),
            '623' : ( 8141, 7544),
            '624' : ( 8124, 7506),
            '625' : ( 8142, 7555),
            '631' : ( 2063, 7689),
            '632' : ( 2079, 7727),
            '633' : ( 2103, 7718),
            '634' : ( 2089, 7681),
            '635' : ( 2109, 7719)
            }

print('\n','-----house_coord-----')
for pnt,coord in house_coord.items():
    xy = np.array( [ [coord[0]] , 
                     [coord[1]] , 
                     [1 ] 
                     ] )
    photo_coord = np.matmul( Affine['AffMat63'] , xy )
    print( '{} :\n {}'.format( pnt,photo_coord) )

print('-----end------')
