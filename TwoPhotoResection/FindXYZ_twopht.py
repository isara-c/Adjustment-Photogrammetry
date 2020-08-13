import numpy as np
from sympy import *
from math import radians, pi
import pickle


import matplotlib.pyplot as plt

with open('EO.pickle', 'rb') as f:
    EO = pickle.load(f)

with open('house_coord_mm.pickle', 'rb') as f:
    house_coord_mm =  pickle.load(f) 


def ParamterAndResidal( A, L ):
    
    E = np.array ( [ [10e-9, 0, 0, 0],
                 [ 0, 10e-9, 0, 0],
                 [ 0, 0, 10e-9, 0],
                 [ 0, 0, 0, 10e-9] ] )

    Q = ( ( 1/10e-5 )**2  ) * E
    P = np.linalg.inv( Q )
    
    AtP = np.matmul( A.T, P )
    
    N = np.matmul( AtP, A ).astype(float)
    U = np.matmul( AtP, L ).astype(float)


    X_ = -np.matmul( np.linalg.inv( N ), U )
    V = np.matmul( A, X_ ) + L
    return ( X_, V, N)


f = 154.006

XpYp_dict, F = house_coord_mm , []

for no in ['62', '63']:
    X, Y, Z  = symbols('X, Y, Z')
    photo = 'photo'+no

    X0, Y0, Z0, ome, phi, kap = EO['EO'+no]

    a11 = cos( phi )*cos( kap )
    a12 = -cos( phi )*sin( kap )
    a13 = sin( phi )
    a21 = cos( ome )*sin( kap )+sin( ome )*sin( phi )*cos( kap )
    a22 = cos( ome )*cos( kap )-sin( ome )*sin( phi )*sin( kap )
    a23 = -sin( ome )*cos( phi )
    a31 = sin( ome )*sin( kap )-cos( ome )*sin( phi )*cos( kap )
    a32 = sin( ome )*cos( kap )+cos( ome )*sin( phi )*sin( kap )
    a33 = cos( ome )*cos( phi )
    
    x = (-f) * ((a11*(X-X0))+(a12*(Y-Y0))+(a13*(Z-Z0)))/(((a31*(X-X0))+(a32*(Y-Y0))+(a33*(Z-Z0))));
    y = (-f) * ((a21*(X-X0))+(a22*(Y-Y0))+(a23*(Z-Z0)))/(((a31*(X-X0))+(a32*(Y-Y0))+(a33*(Z-Z0))));

    F.append( [x,y] )


print('\n', 'project name : Two Photo Resection Image 62 & 63', '\n' )
print('Focal Length : {} mm. '.format(f))
print('All point counted : {} '.format( '4' )  )
print('Unknowns : {} '.format( '3' ), '\n' )
print( '================== Coord House62&63 in mm. ====================')
for key, value in XpYp_dict.items():
    print( 'Coor Photo{} of pnt{} |   Xp : {:8.3f} mm.    Yp :{:8.3f} mm.'.format( key[:-1], key[-1], value[0], value[1] ))
print('\n')



xyz0 = [4000, 2400, 200]
result_house = []
result_V = []
result_Si = []

for i in range(1,6):
    
    Fx62, Fy62, Fx63, Fy63 = F[0][0], F[0][1], F[1][0], F[1][1]
    Lb = np.array( [ XpYp_dict['62'+str(i)][0], XpYp_dict['62'+str(i)][1], XpYp_dict['63'+str(i)][0], XpYp_dict['63'+str(i)][1] ] )
    L0 = np.array( [ Fx62, Fy62, Fx63, Fy63])
    l = L0 - Lb

    error = 1;
    while error >= 1e-5 :
        L = np.transpose(  np.array( [  l[0].subs(list(zip( [X, Y, Z],xyz0 ))),
                                        l[1].subs(list(zip( [X, Y, Z],xyz0 ))),
                                        l[2].subs(list(zip( [X, Y, Z],xyz0 ))),
                                        l[3].subs(list(zip( [X, Y, Z],xyz0 ))) ]) )

        
        J = []
        J.append([          diff( Fx62, X ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fx62, Y ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fx62, Z ).subs(list(zip( [X, Y, Z],xyz0 ))) ] )
        
        J.append([          diff( Fy62, X ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fy62, Y ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fy62, Z ).subs(list(zip( [X, Y, Z],xyz0 )))] )
                           
        J.append( [         diff( Fx63, X ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fx63, Y ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fx63, Z ).subs(list(zip( [X, Y, Z],xyz0 ))) ] )
        
        J.append([          diff( Fy63, X ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fy63, Y ).subs(list(zip( [X, Y, Z],xyz0 ))),
                            diff( Fy63, Z ).subs(list(zip( [X, Y, Z],xyz0 )))] )
        A = np.array( J )
     
        X_, V, N  = ParamterAndResidal( A, L )

        
        Xa = X_ +  np.array( xyz0 )
        
        error = abs( Xa[1] - xyz0[1] );
        xyz0 = [ Xa[0],Xa[1], Xa[2] ]

    House_Ground = xyz0
    xyz0 = [3000, 3000, 200]

    VtV = np.matmul( V.T, V )
    S0 = sqrt( VtV )
    Qxx = ( np.linalg.inv( N ) )
    
    Si_X, Si_Y, Si_Z  = S0 * np.sqrt(Qxx[0][0]), S0 * np.sqrt(Qxx[1][1]), S0 * np.sqrt(Qxx[2][2])

    result_house.append(House_Ground)
    result_V.append(V.tolist())
    result_Si.append( [Si_X, Si_Y, Si_Z] )

print( '============== Ground Coord & Adjusted Value ==================')
for i, value in enumerate( result_house ):
    print('House{} : X = {:8.4f} m.  |  std. : {:>3}{:.3f} m. '.format(i+1, value[0], u"\u00B1", result_Si[i][0]))
    print('House{} : Y = {:8.4f} m.  |  std. : {:>3}{:.3f} m. '.format(i+1, value[1], u"\u00B1", result_Si[i][1]))
    print('House{} : Z = {:8.4f}  m.  |  std. : {:>3}{:.3f} m. '.format(i+1, value[2], u"\u00B1", result_Si[i][2]), '\n')

print('\n', '================ Measurement Residues ===================')
for j, value in enumerate( result_V ):
    print ('photo_{} : pnt {}| vx :  {:+.4f} mm.  vy :  {:+.4f} mm.  | vx :  {:+.1f} pix  vy :  {:+.1f} pix'.format( \
        '62' , j,  value[0] , value[1], value[0]* 15232 /230 , value[1]* 15232 /230))

print('\n')
for k, value in enumerate( result_V ):
    print ('photo_{} : pnt {} | vx :  {:+.4f} mm.  vy :  {:+.4f} mm.  | vx :  {:+.1f} pix  vy :  {:+.1f} pix'.format( \
        '63' , k,  value[2] , value[3], value[2]* 15232 /230 , value[3]* 15232 /230)  )
    
    


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X_house = [ result_house[0][0], result_house[1][0], result_house[2][0], result_house[3][0], result_house[0][0]]
Y_house = [ result_house[0][1], result_house[1][1], result_house[2][1], result_house[3][1], result_house[0][1]]
Z_house = np.array([
                [result_house[0][2], result_house[1][2], result_house[2][2], result_house[3][2], result_house[0][2] ],
                [result_house[0][2], result_house[1][2], result_house[2][2], result_house[3][2], result_house[0][2] ]] )

ax.plot_wireframe(X_house,Y_house, Z_house, color = 'black')

Z_ground = np.array([
                [result_house[-1][2]],
                [result_house[-1][2]] ])
ax.scatter(result_house[4][0], result_house[4][1],Z_ground )
plt.show()






























