import numpy as np

def ParamterAndResidal( A, L ):
    N = np.matmul( np.transpose( A ), A )
    U = np.matmul( np.transpose( A ), L )
    Xa = np.matmul( -1*np.linalg.inv( N ), U )
    V = np.matmul( A, Xa ) + L
    return ( Xa, V )

def ScreenCoord2PhotCoord( parameter ):
    a, b, c, d, e, f = parameter
    Pa = e / ( a*e - b*d )
    Pb = -b  / ( a*e - b*d )
    Pc = -( c*e - b*f ) / ( a*e - b*d )
    Pd = -d / ( a*e - b*d )
    Pe = a / ( a*e - b*d )
    Pf = ( c*d - a*f ) / ( a*e - b*d )
    return ( Pa, Pb, Pc, Pd, Pe, Pf )


A = np.array( [ [113.003,       0.006,       1,       0,           0,           0],
     [0,             0,           0,       113.003,     0.006,       1],
     [-112.987,      0.006,       1,       0,           0,           0],
     [0,             0,           0,       -112.987,    0.006,       1],
     [0.012,         113.023,     1,       0,           0,           0],
     [0,             0,           0,       0.012,       113.023,     1],
     [0.012,	    -113.004,     1,       0,           0,           0],
     [0,             0,           0,       0.012,       -113.004,    1],
     [113.007,	    113.005,      1,       0,           0,           0],
     [0,             0,           0,       113.007,     113.005,     1],
     [-112.989,	    -112.994,     1,       0,           0,           0],
     [0,             0,           0,       -112.989,    -112.994,    1],
     [-112.998,	    113.006,      1,       0,           0,           0],
     [0,             0,           0,       -112.998,    113.006,     1],
     [112.999,	    -112.998,     1,       0,           0,           0],
     [0,             0,           0,       112.999,     -112.998,    1] ] )


L_62 = -1*np.array( [ [   15088],	
        [7561],
        [34], 
        [7561],	
        [7561],	
        [32],
        [7560],	
        [15088],
        [15088],
        [33],	
        [33],	
        [15088],
        [33],
        [33],
        [15088],
        [15088] ] )

L_63 = -1*np.array( [ [   15088],	
        [7560],
        [34], 
        [7561],	
        [7561],	
        [32],
        [7560],	
        [15088],
        [15088],
        [34],	
        [33],	
        [15088],
        [33],
        [33],
        [15088],
        [15088] ] )

Xa_62, V_62 = ParamterAndResidal( A, L_62 )
Xa_63, V_63 = ParamterAndResidal( A, L_63 )

print( 'parameter of img 62', Xa_62, '\n',  'residal of img 62', V_62 ,'\n')
print(  '\n', 'parameter of img 63', Xa_63, '\n',  'residal of img 63', V_63 ,'\n' ) 

print( 'parameter of img 62 for convert to photo coord')
for item in ScreenCoord2PhotCoord( Xa_62 ) :print( float(item) )
print(  '\n', 'parameter of img 63 for convert to photo coord')
for item in ScreenCoord2PhotCoord( Xa_63 ) : print( float(item) )




























