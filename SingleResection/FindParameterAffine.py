import numpy as np
import pickle

AffMat62 = np.array([ [  1.5011859226653e-02, -2.8365288601049e-7, -113.50077766574],
                    [  -1.8716778648295e-7, -1.5012983008061e-02, 113.50467204681]  ] )

AffMat63 = np.array([ [  1.5011224390783e-02, 1.4843241403917e-07, -113.49751765623],
                    [  -2.9018885425355e-07, -1.501169666325e-02, 0.11349384986359e+03]  ] )

Affine = { 'AffMat62' : AffMat62 ,
           'AffMat63' : AffMat63 }

with open('Affine.pickle', 'wb') as f:
    pickle.dump( Affine, f, pickle.HIGHEST_PROTOCOL )
