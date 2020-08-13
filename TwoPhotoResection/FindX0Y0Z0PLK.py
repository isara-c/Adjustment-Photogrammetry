import numpy as np
import sympy
from sympy import *
from math import radians, pi
import pickle


with open('GCPs.pickle', 'rb') as f:
    GCPs = pickle.load(f)

def ParamterAndResidal( A, L ):
    N = np.matmul( np.transpose( A ), A )
    U = np.matmul( np.transpose( A ), L )
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
    

    a11 = cos( phi )*cos( kap )
    a12 = -cos( phi )*sin( kap )
    a13 = sin( phi )
    a21 = cos( ome )*sin( kap )+sin( ome )*sin(phi)*cos(kap)
    a22 = cos( ome )*cos( kap )-sin( ome )*sin(phi)*sin(kap)
    a23 = -sin( ome )*cos( phi )
    a31 = sin( ome )*sin( kap )-cos( ome )*sin( phi)*cos(kap)
    a32 = sin( ome )*cos( kap )+cos( ome )*sin( phi)*sin(kap)
    a33 = cos( ome )*cos( phi )
    
    X1, Y1, x_1, y_1, Z1  =     GCPs['ground']['20301'][0], GCPs['ground']['20301'][1], GCPs[photo]['20301'][0], GCPs[photo]['20301'][1], GCPs['ground']['20301'][2]
    X2, Y2, x_2, y_2, Z2  =     GCPs['ground']['20401'][0], GCPs['ground']['20401'][1], GCPs[photo]['20401'][0], GCPs[photo]['20401'][1], GCPs['ground']['20401'][2]
    X3, Y3, x_3, y_3, Z3  =     GCPs['ground']['20501'][0], GCPs['ground']['20501'][1], GCPs[photo]['20501'][0], GCPs[photo]['20501'][1], GCPs['ground']['20501'][2]
    X4, Y4, x_4, y_4, Z4  =     GCPs['ground']['30301'][0], GCPs['ground']['30301'][1], GCPs[photo]['30301'][0], GCPs[photo]['30301'][1], GCPs['ground']['30301'][2]
    X5, Y5, x_5, y_5, Z5  =     GCPs['ground']['30401'][0], GCPs['ground']['30401'][1], GCPs[photo]['30401'][0], GCPs[photo]['30401'][1], GCPs['ground']['30401'][2]
    X6, Y6, x_6, y_6, Z6  =     GCPs['ground']['40301'][0], GCPs['ground']['40301'][1], GCPs[photo]['40301'][0], GCPs[photo]['40301'][1], GCPs['ground']['40301'][2]
    X7, Y7, x_7, y_7, Z7  =     GCPs['ground']['40401'][0], GCPs['ground']['40401'][1], GCPs[photo]['40401'][0], GCPs[photo]['40401'][1], GCPs['ground']['40401'][2]
    X8, Y8, x_8, y_8, Z8  =     GCPs['ground']['40501'][0], GCPs['ground']['40501'][1], GCPs[photo]['40501'][0], GCPs[photo]['40501'][1], GCPs['ground']['40501'][2]


    x1 = (-f) * ((a11*(X1-X0))+(a12*(Y1-Y0))+(a13*(Z1-Z0)))/(((a31*(X1-X0))+(a32*(Y1-Y0))+(a33*(Z1-Z0))));
    y1 = (-f) * ((a21*(X1-X0))+(a22*(Y1-Y0))+(a23*(Z1-Z0)))/(((a31*(X1-X0))+(a32*(Y1-Y0))+(a33*(Z1-Z0))));
    x2 = (-f) * ((a11*(X2-X0))+(a12*(Y2-Y0))+(a13*(Z2-Z0)))/(((a31*(X2-X0))+(a32*(Y2-Y0))+(a33*(Z2-Z0))));
    y2 = (-f) * ((a21*(X2-X0))+(a22*(Y2-Y0))+(a23*(Z2-Z0)))/(((a31*(X2-X0))+(a32*(Y2-Y0))+(a33*(Z2-Z0))));
    x3 = (-f) * ((a11*(X3-X0))+(a12*(Y3-Y0))+(a13*(Z3-Z0)))/(((a31*(X3-X0))+(a32*(Y3-Y0))+(a33*(Z3-Z0))));
    y3 = (-f) * ((a21*(X3-X0))+(a22*(Y3-Y0))+(a23*(Z3-Z0)))/(((a31*(X3-X0))+(a32*(Y3-Y0))+(a33*(Z3-Z0))));
    x4 = (-f) * ((a11*(X4-X0))+(a12*(Y4-Y0))+(a13*(Z4-Z0)))/(((a31*(X4-X0))+(a32*(Y4-Y0))+(a33*(Z4-Z0))));  
    y4 = (-f) * ((a21*(X4-X0))+(a22*(Y4-Y0))+(a23*(Z4-Z0)))/(((a31*(X4-X0))+(a32*(Y4-Y0))+(a33*(Z4-Z0))));
    x5 = (-f) * ((a11*(X5-X0))+(a12*(Y5-Y0))+(a13*(Z5-Z0)))/(((a31*(X5-X0))+(a32*(Y5-Y0))+(a33*(Z5-Z0))));
    y5 = (-f) * ((a21*(X5-X0))+(a22*(Y5-Y0))+(a23*(Z5-Z0)))/(((a31*(X5-X0))+(a32*(Y5-Y0))+(a33*(Z5-Z0))));
    x6 = (-f) * ((a11*(X6-X0))+(a12*(Y6-Y0))+(a13*(Z6-Z0)))/(((a31*(X6-X0))+(a32*(Y6-Y0))+(a33*(Z6-Z0))));
    y6 = (-f) * ((a21*(X6-X0))+(a22*(Y6-Y0))+(a23*(Z6-Z0)))/(((a31*(X6-X0))+(a32*(Y6-Y0))+(a33*(Z6-Z0))));
    x7 = (-f) * ((a11*(X7-X0))+(a12*(Y7-Y0))+(a13*(Z7-Z0)))/(((a31*(X7-X0))+(a32*(Y7-Y0))+(a33*(Z7-Z0))));
    y7 = (-f) * ((a21*(X7-X0))+(a22*(Y7-Y0))+(a23*(Z7-Z0)))/(((a31*(X7-X0))+(a32*(Y7-Y0))+(a33*(Z7-Z0))));
    x8 = (-f) * ((a11*(X8-X0))+(a12*(Y8-Y0))+(a13*(Z8-Z0)))/(((a31*(X8-X0))+(a32*(Y8-Y0))+(a33*(Z8-Z0))));
    y8 = (-f) * ((a21*(X8-X0))+(a22*(Y8-Y0))+(a23*(Z8-Z0)))/(((a31*(X8-X0))+(a32*(Y8-Y0))+(a33*(Z8-Z0))));

    x = [x1, x2, x3, x4, x5, x6, x7, x8]
    y = [y1, y2, y3, y4, y5, y6, y7, y8]

    L  = np.transpose( np.array( [x_1, y_1, x_2,y_2, x_3, y_3, x_4, y_4, x_5, y_5, x_6, y_6, x_7, y_7, x_8, y_8] ) );
    L0 = np.transpose( np.array( [x1, y1, x2, y2, x3, y3, x4, y4, x5, y5, x6, y6, x7, y7, x8, y8] ) );
    l  = L0 - L

    parameter_EO = [ 3700, 2100, 2250, radians(2), radians(0), radians(2)]

    while True :
        A = []
        l  = L0 - L
        for i in range(8):
            Fx_x0   =   diff( x[i], X0 )
            Fx_y0   =   diff( x[i], Y0 )
            Fx_z0   =   diff( x[i], Z0 )
            Fx_ome  =   diff( x[i], ome)
            Fx_phi  =   diff( x[i], phi)
            Fx_kap  =   diff( x[i], kap)
            Fy_x0   =   diff( y[i], X0 )
            Fy_y0   =   diff( y[i], Y0 )
            Fy_z0   =   diff( y[i], Z0 )
            Fy_ome  =   diff( y[i], ome )
            Fy_phi  =   diff( y[i], phi )
            Fy_kap  =   diff( y[i], kap )
        
            A.append([Fx_x0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fx_y0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fx_z0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fx_ome.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                            Fx_phi.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                            Fx_kap.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))] )
                           
            A.append( [Fy_x0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fy_y0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fy_z0.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fy_ome.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fy_phi.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))) ,
                            Fy_kap.subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))] )
        
        A = np.array( A )
        l = np.transpose( np.array( [l[0].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[1].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[2].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[3].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[4].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[5].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[6].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[7].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[8].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[9].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[10].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[11].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[12].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[13].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[14].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO))),
                                l[15].subs(list(zip([X0,Y0,Z0,ome,phi,kap],parameter_EO)))] ) )
    
        X, V, N = ParamterAndResidal( A, l )
        Xa = X +  parameter_EO
        d = abs(Xa - parameter_EO)
        parameter_EO = Xa
        
        EO['EO' + no] = parameter_EO 
    
        if d[3] and d[4] and d[5]  < 0.00001 :
            break
    VtV = np.matmul( V.T, V )
    S0 = sqrt( VtV / 10  )

    Qxx =   ( np.linalg.inv( N ) )
    
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

        print ('photo_{} | vx :  {:+.4f} mm.  vy :  {:+.4f} mm.'.format( name ,  V[i*2] , V[i*2+1]) )

    print( '-------------------------------------------------' )

with open('EO.pickle', 'wb') as f:
    pickle.dump( EO, f, pickle.HIGHEST_PROTOCOL )


