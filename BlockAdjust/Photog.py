#
#
#  Photog : Photogrammetry object class
#  8 Mar 2017 : Phisan Santitamnot
#
#
import numpy as np
import math

PI2 = 2.*math.pi
def chk_ang( ang_rad ):
    if ang_rad<-PI2 or ang_rad>+PI2:
        print('***WARNING** chk_ang {}'.format( ang_rad ) )
        #raise ValueError('***ERROR*** Rotation.__init__() angle {} must be radian and within +/-2.*PI'.format(ang_rad) )
    return

class Photo():
    def __init__(self, Name, FocLen, POS, ROT, SEQUENCE='PIX4D'):
        self.Name = Name
        (self.Ome, self.Phi, self.Kap) = ROT
        map( chk_ang, ROT )
        self.FocLen = FocLen
        (self.x0 , self.y0)           = 0.0, 0.0
        (self.X0 , self.Y0 , self.Z0) = POS
        self.SEQUENCE =SEQUENCE
        return

    #####################################################################
    def World2Photo(self, XYZ ):
        (X,Y,Z) = XYZ
        X0,Y0,Z0    = self.X0, self.Y0, self.Z0
        f,x0,y0     = self.FocLen, self.x0, self.y0
        (r11,r12,r13,r21,r22,r23,r31,r32,r33) =\
                self.RotationMatrix((self.Ome, self.Phi, self.Kap))
        nom_xp = r11*(X-X0)+r21*(Y-Y0)+r31*(Z-Z0)
        nom_yp = r12*(X-X0)+r22*(Y-Y0)+r32*(Z-Z0)
        denom  = r13*(X-X0)+r23*(Y-Y0)+r33*(Z-Z0)
        xp = x0 - f*nom_xp/denom
        yp = y0 - f*nom_yp/denom
        return (xp,yp)

    ####################################################################
    def __str__(self):
        """ return string of omega,phi,kappa"""
        pos_ori = 'XYZ=({:.3f},{:.3f},{:.3f}) opk=({:.8f},{:.8f},{:.8f})' .format (
                self.X0 , self.Y0 , self.Z0,
                math.degrees(self.Ome),math.degrees(self.Phi),math.degrees(self.Kap) )
        return pos_ori

    ####################################################################
    def RotationMatrix(self, OmePhiKap):
        (Ome,Phi,Kap) = OmePhiKap
        coso = math.cos(Ome)
        sino = math.sin(Ome)
        cosp = math.cos(Phi)
        sinp = math.sin(Phi)
        cosk = math.cos(Kap)
        sink = math.sin(Kap)
        # ERDAS
        if self.SEQUENCE=='ERDAS':
            r11 = cosp * cosk - sinp * sino * sink
            r12 = -sinp * coso
            r13 = cosp * sink + sinp * sino * cosk
            r21 = sinp * cosk + cosp * sino * sink
            r22 = cosp * coso
            r23 = sinp * sink - cosp * sino * cosk
            r31 = -coso * sink
            r32 =  sino
            r33 =  coso * cosk
        # Bonn Dewitt
        elif self.SEQUENCE=='BONN':
            r11 = cosp * cosk;
            r12 = sino * sinp * cosk + coso * sink;
            r13 = -coso * sinp * cosk + sino * sink;
            r21 = -cosp * sink;
            r22 = -sino * sinp * sink + coso * cosk;
            r23 = coso * sinp * sink + sino * cosk;
            r31 = sinp;
            r32 = -sino * cosp;
            r33 = coso * cosp;
        # Pix4D, Digital Photogrammetry (P.Santitamnont)
        else:
            r11 = cosk*cosp
            r12 = -sink*cosp
            r13 = sinp
            r21 = cosk*sino*sinp+sink*coso
            r22 = cosk*coso-sink*sino*sinp
            r23 = -sino*cosp
            r31 = sink*sino-cosk*coso*sinp
            r32 = sink*coso*sinp+cosk*sino
            r33 = coso*cosp
        return (r11,r12,r13,r21,r22,r23,r31,r32,r33)

################################################################################
if __name__ == "__main__":
    # photo 62
    Solu = {     'XL'   :  3709.1004, 'YL'   :  2100.6322, 'ZL'   :  2258.6875,
    'Ome'  :  math.radians( 2.1947),  'Phi'  :  math.radians( -0.4220),  'Kap'  :  math.radians( -2.2214) }
    PhotMeasFile = r'D:\Phisan\Photogrammetry\PhotogrammetryII\Lab-2016\Lab_4_data\InputResect62_19_20.dat'

    ##  photo 63
    #Solu = {     'XL'   :  4907.812, 'YL'   :  2088.579, 'ZL'   :  2257.266,
    #'Ome'  :  math.radians(2.42343 ), 'Phi'  :  math.radians(-0.31138 ), 'Kap'  :  math.radians(-1.78346) }
    #PhotMeasFile = r'D:\Phisan\Photogrammetry\PhotogrammetryII\Lab-2016\Lab_4_data\InputResect63_19_20.dat'

    with open(PhotMeasFile,'r') as fd:
        MeasLines = fd.readlines()
    FocLen = float( MeasLines[0] )  # f in mm.

    photo_meas = np.genfromtxt(PhotMeasFile, skip_header=1, delimiter=None,
            dtype=[('f0','S5'),('f1','<f8'),('f2','<f8'),('f3','<f8'),('f4','<f8'),('f5','<f8')])
    print('photo measurement : {}\n'.format( photo_meas ) ) 


    photo = Photo("62", FocLen, (Solu['XL'], Solu['YL'], Solu['ZL'] ),
                    (Solu['Ome'],Solu['Phi'],Solu['Kap'] ) )

    for i in photo_meas:
        name,xp,yp,X,Y,Z = i
        xp_,yp_ = photo.World2Photo( (X,Y,Z) )
        dxp,dyp = xp_-xp,yp_-yp
        print('{:6s}{:10.3f}{:10.3f}{:10.3f}{:10.3f}{:10.4f}{:10.4f}'.\
                format(name,xp_,yp_,xp,yp,dxp,dyp) )


    print( 'Photo = {}'.format( str(photo) ) )
    print( 'Position (internal) = {} : {} : {} '.format( photo.X0,photo.Y0,photo.Z0 ) )
    print( 'Rotation (internal) = {} : {} : {} '.format( photo.Ome,photo.Phi,photo.Kap ) )

