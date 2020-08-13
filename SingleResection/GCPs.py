import pickle

photo62 = {'20301' : ( -0.405, -86.545 ),
        '20401' : ( 39.031,  -85.990 ),
        '20501' : ( 97.067, -87.011 ),  
        '30301' : ( -14.980, -20.217 ),
        '30401' : ( 38.087,  -5.775 ),
        '40301' : ( -13.718, 73.449 ),
        '40401' : ( 34.215,  77.831 ),
        '40501' : ( 94.368,  81.149 ),
         }

photo63 = {'20301' : ( -93.110,  -89.585 ),
        '20401' : ( -54.922,  -89.420 ),
        '20501' : ( 1.821,   -91.027 ),  
        '30301' : (  -106.621, -22.722 ),
        '30401' : ( -51.800,  -8.612 ),
        '40301' : ( -104.821, 71.131 ),
        '40401' : ( -56.169, 75.197 ),
        '40501' : ( 5.692,   78.154 ),
         }

ground = { '20301' : (3674.942,    1050.052,    206.364),
        '20401' : (4182.090,    1053.451,   240.941),
        '20501' : (4906.054,    1030.005,    279.686),
        '30301' : (3517.415,    1921.055,    228.798),
        '30401' : (4228.470,    2082.750,    207.739),
        '40301' : (3581.964,    3142.402,    285.339),
        '40401' : (4213.013,    3182.528,    277.559),
        '40501' : (5027.192,    3214.013,    250.182),
         }

GCPs = { 'photo62' : photo62 ,
         'photo63' : photo63 ,
           'ground' : ground }

with open('GCPs.pickle', 'wb') as f:
    pickle.dump( GCPs, f, pickle.HIGHEST_PROTOCOL )
