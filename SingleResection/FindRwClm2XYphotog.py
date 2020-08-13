import numpy as np
import pickle

with open('Affine.pickle', 'rb') as f: # from older Lab Affine para
    Affine =  pickle.load(f) 

ImageCoord62 = { '40301' : (6647, 2668), # Colums, Row at GCPs  
        '40401' : (9840, 2376),   
        '40501' : (13847, 2155),   
        '30301' : (6563, 8907),   
        '30401' : (10098, 7945),   
        '20301' : (7534, 13325),   
        '20401' : (10161, 13288),
        '20501' : (14027, 13356)
       }

ImageCoord63 = { '40301' : (578, 2822),   
        '40401' : (3819, 2551),
        '40501' : (7940, 2354),   
        '30301' : (458, 9074),   
        '30401' : (4110, 8134),   
        '20301' : (1358, 13528),   
        '20401' : (3902, 13517),  
        '20501' : (7682, 13624)
       }

print('----- Pic 62-----')
for pnt,coord in ImageCoord62.items(): 
    xy = np.array( [ [coord[0]] , 
                     [coord[1]] , 
                     [1 ] 
                     ] )

    # transform Coord of IMG to Photo with Affine
    photo_coord = np.matmul( Affine['AffMat62'] , xy )
    
    print( '{} :\n {}'.format( pnt,photo_coord) )

print('\n','-----Pic 63-----')
for pnt,coord in ImageCoord63.items():
    xy = np.array( [ [coord[0]] , 
                     [coord[1]] , 
                     [1 ] 
                     ] )
    photo_coord = np.matmul( Affine['AffMat63'] , xy )
    print( '{} :\n {}'.format( pnt,photo_coord) )

print('-----end------')
