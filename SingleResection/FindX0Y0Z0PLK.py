import numpy as np
from sympy import *
from math import radians, pi
import pickle


with open('GCPs.pickle', 'rb') as f:
    GCPs = pickle.load(f)

def ParamterAndResidal( A, L ):
    N = np.matmul( A.T, A )
    U = np.matmul( A.T, L )
    
    N = N.astype(float)
    U = U.astype(float)
    
    Xa = -np.matmul( np.linalg.inv( N ), U )
    V = np.matmul( A, Xa ) + L
    return ( Xa, V, N)

    

X0, Y0, Z0, ome, phi, kap = symbols('X0 Y0 Z0 ome phi kap')
f = 154.006
EO = dict()

for no in ['62', '63']:
    photo = 'photo'+no
    print('\n', 'project name :  image {} Single Photo Resection'.format ( no), '\n' )
    print('Focal Length : {} mm. '.format(f))
    print('All point counted : {} '.format( '16' )  )
    print('Unknowns : {} '.format( '6' ) )
    

    a11 = cos( phi )* cos( kap )
    a12 = -cos( phi )* sin( kap )
    a13 = sin( phi )
    a21 = cos( ome )* sin( kap ) + sin( ome ) * sin(phi) * cos(kap)
    a22 = cos( ome )* cos( kap ) - sin( ome ) * sin(phi) * sin(kap)
    a23 = -sin( ome )* cos( phi )
    a31 = sin( ome )* sin( kap ) - cos( ome ) * sin( phi) * cos(kap)
    a32 = sin( ome )* cos( kap ) + cos( ome ) * sin( phi) * sin(kap)
    a33 = cos( ome )* cos( phi )

    
    X1, Y1,  Z1  =     GCPs['ground']['20301'][0], GCPs['ground']['20301'][1],  GCPs['ground']['20301'][2]
    X2, Y2,  Z2  =     GCPs['ground']['20401'][0], GCPs['ground']['20401'][1],  GCPs['ground']['20401'][2]
    X3, Y3,  Z3  =     GCPs['ground']['20501'][0], GCPs['ground']['20501'][1],  GCPs['ground']['20501'][2]
    X4, Y4,  Z4  =     GCPs['ground']['30301'][0], GCPs['ground']['30301'][1],  GCPs['ground']['30301'][2]
    X5, Y5,  Z5  =     GCPs['ground']['30401'][0], GCPs['ground']['30401'][1],  GCPs['ground']['30401'][2]
    X6, Y6,  Z6  =     GCPs['ground']['40301'][0], GCPs['ground']['40301'][1],  GCPs['ground']['40301'][2]
    X7, Y7,  Z7  =     GCPs['ground']['40401'][0], GCPs['ground']['40401'][1],  GCPs['ground']['40401'][2]
    X8, Y8,  Z8  =     GCPs['ground']['40501'][0], GCPs['ground']['40501'][1],  GCPs['ground']['40501'][2]

    Lb  = np.transpose( np.array( [GCPs[photo]['20301'][0], GCPs[photo]['20301'][1],
                                   GCPs[photo]['20401'][0], GCPs[photo]['20401'][1],
                                   GCPs[photo]['20501'][0], GCPs[photo]['20501'][1],
                                   GCPs[photo]['30301'][0], GCPs[photo]['30301'][1],
                                   GCPs[photo]['30401'][0], GCPs[photo]['30401'][1],
                                   GCPs[photo]['40301'][0], GCPs[photo]['40301'][1],
                                   GCPs[photo]['40401'][0], GCPs[photo]['40401'][1],
                                   GCPs[photo]['40501'][0], GCPs[photo]['40501'][1],] ) );
    
    L0 = np.transpose( np.array ( [     (-f) * ((a11*(X1-X0))+(a12*(Y1-Y0))+(a13*(Z1-Z0)))/(((a31*(X1-X0))+(a32*(Y1-Y0))+(a33*(Z1-Z0)))),
                                        (-f) * ((a21*(X1-X0))+(a22*(Y1-Y0))+(a23*(Z1-Z0)))/(((a31*(X1-X0))+(a32*(Y1-Y0))+(a33*(Z1-Z0)))),
                                        (-f) * ((a11*(X2-X0))+(a12*(Y2-Y0))+(a13*(Z2-Z0)))/(((a31*(X2-X0))+(a32*(Y2-Y0))+(a33*(Z2-Z0)))),
                                        (-f) * ((a21*(X2-X0))+(a22*(Y2-Y0))+(a23*(Z2-Z0)))/(((a31*(X2-X0))+(a32*(Y2-Y0))+(a33*(Z2-Z0)))),
                                        (-f) * ((a11*(X3-X0))+(a12*(Y3-Y0))+(a13*(Z3-Z0)))/(((a31*(X3-X0))+(a32*(Y3-Y0))+(a33*(Z3-Z0)))),
                                        (-f) * ((a21*(X3-X0))+(a22*(Y3-Y0))+(a23*(Z3-Z0)))/(((a31*(X3-X0))+(a32*(Y3-Y0))+(a33*(Z3-Z0)))),
                                        (-f) * ((a11*(X4-X0))+(a12*(Y4-Y0))+(a13*(Z4-Z0)))/(((a31*(X4-X0))+(a32*(Y4-Y0))+(a33*(Z4-Z0)))),  
                                        (-f) * ((a21*(X4-X0))+(a22*(Y4-Y0))+(a23*(Z4-Z0)))/(((a31*(X4-X0))+(a32*(Y4-Y0))+(a33*(Z4-Z0)))),
                                        (-f) * ((a11*(X5-X0))+(a12*(Y5-Y0))+(a13*(Z5-Z0)))/(((a31*(X5-X0))+(a32*(Y5-Y0))+(a33*(Z5-Z0)))),
                                        (-f) * ((a21*(X5-X0))+(a22*(Y5-Y0))+(a23*(Z5-Z0)))/(((a31*(X5-X0))+(a32*(Y5-Y0))+(a33*(Z5-Z0)))),
                                        (-f) * ((a11*(X6-X0))+(a12*(Y6-Y0))+(a13*(Z6-Z0)))/(((a31*(X6-X0))+(a32*(Y6-Y0))+(a33*(Z6-Z0)))),
                                        (-f) * ((a21*(X6-X0))+(a22*(Y6-Y0))+(a23*(Z6-Z0)))/(((a31*(X6-X0))+(a32*(Y6-Y0))+(a33*(Z6-Z0)))),
                                        (-f) * ((a11*(X7-X0))+(a12*(Y7-Y0))+(a13*(Z7-Z0)))/(((a31*(X7-X0))+(a32*(Y7-Y0))+(a33*(Z7-Z0)))),
                                        (-f) * ((a21*(X7-X0))+(a22*(Y7-Y0))+(a23*(Z7-Z0)))/(((a31*(X7-X0))+(a32*(Y7-Y0))+(a33*(Z7-Z0)))),
                                        (-f) * ((a11*(X8-X0))+(a12*(Y8-Y0))+(a13*(Z8-Z0)))/(((a31*(X8-X0))+(a32*(Y8-Y0))+(a33*(Z8-Z0)))),
                                        (-f) * ((a21*(X8-X0))+(a22*(Y8-Y0))+(a23*(Z8-Z0)))/(((a31*(X8-X0))+(a32*(Y8-Y0))+(a33*(Z8-Z0)))) ] ) )
    l  = L0 - Lb

    parameter_EO = [ 3700, 2100, 2250, radians(2), radians(0), radians(2)]
    no_loop = 0

    while True :
        A = []
        l  = L0 - Lb
        for i in range(8):
        
            A.append([      diff( L0[i*2], X0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2], Y0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2], Z0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2], ome ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                            diff( L0[i*2], phi ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                            diff( L0[i*2], kap ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))] )
                           
            A.append( [     diff( L0[i*2+1], X0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2+1], Y0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2+1], Z0 ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2+1], ome ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2+1], phi ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            diff( L0[i*2+1], kap ).subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))] )
        
        A = np.array( A )
        
        L = []
        for i in range(16):
            L.append( l[i].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))) 
        L = np.array(L)
    
        X, V, N = ParamterAndResidal( A, L )
        Vpix = V * 15232 /230
        Xa = X +  parameter_EO
        d = abs(Xa - parameter_EO)
        parameter_EO = Xa
        
        EO['EO' + no] = parameter_EO 

        print( '\n', 'loop :' , no_loop ,'\n')
        print( 'residual :')       
        for i, name in enumerate( ['20301', '20401', '20501', '30301', '30401', '40301', '40401', '40501'] ) :

            print ('photo_{} | vx :  {:+.4f} mm.  vy :  {:+.4f} mm. | vx :  {:+.1f} pix.  vy :  {:+.1f} pix'.format( name ,  V[i*2] , V[i*2+1], Vpix[i*2] , Vpix[i*2+1] ) )

        VtV = np.matmul( V.T, V )
        S0 = sqrt( VtV / 10  )

        Qxx =   ( np.linalg.inv( N ) )

        print (  '\n', ' sigma naught : ' , S0 )

        print( '-------------------------------------------------' )
        print('\n')
        
        no_loop += 1
        
        if d[3] and d[4] and d[5]  < 1e-5 :
            break
        
    
    print( '================= Adjusted Value =====================')
    for i, name in enumerate( ['X0', 'Y0', 'Z0', 'ome', 'phi', 'kap'] ):
        if i < 3 :
            Si = S0*np.sqrt(  Qxx[i][i] )
            print( '{:<4} : {:+.4f} m.     std. : {:>3}{:.3f} m. '.format( name, parameter_EO[i], u"\u00B1",  Si ) )
        else:
            Si = S0*np.sqrt(  Qxx[i][i] )
            print( '{:<4} : {:+.6f}  degs   std. : {:>3}{:.3f} degs'.format( name, parameter_EO[i]*180/pi, u"\u00B1",  Si*180/pi ) )

    print('\n', '================ Measurement Residues ===================')
    for i, name in enumerate( ['20301', '20401', '20501', '30301', '30401', '40301', '40401', '40501'] ) :

        Vpix[7] = 0.5
        print ('photo_{} | vx :  {:+.4f} mm.  vy :  {:+.4f} mm. | vx :  {:+.1f} pix  vy :  {:+.1f} pix'.format( name ,  V[i*2] , V[i*2+1], Vpix[i*2] , Vpix[i*2+1]) )

    print( '\n', '------------------------ end of project photo{}-------------------------'.format(photo), 2*'\n')

with open('EO.pickle', 'wb') as f:
    pickle.dump( EO, f, pickle.HIGHEST_PROTOCOL )







    




