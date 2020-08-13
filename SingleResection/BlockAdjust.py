#
#
PROG= \
"""BlockAdjust : Bundle Block Adjustment, reads a YAML photogrammetric 
   block file and computes all exterior orientations (EOP) of the photos.
   Optional EOPs will be read and applied to specified photo.
   If measurements on tie points are existed, it will estimate 3-D coordinate 
   of tie points too.
   Author  : Phisan Santitamnont (phisan.chula@gmail.com)
   History : 23 Mar 2017 : Initialed , read INI format
             23 Feb 2018 : version 0.7, read YAML format instead 
             1  Mar 2018 : version 0.75, YAML attribute could be separated by 
                                         SPC,COMMA,COLON,SEMI-COLON all handled
                                         by SplitFloat() function.  
             23 Jan 2019 : version 0.8, porting to Python3 with Python2 compatible
             """
#
#
import numpy as np
import math,re,sys
import six
import yaml, yamlordereddictloader
from collections import OrderedDict,defaultdict,Counter
from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import argparse

from Photog import *

def SplitFloat( float_str, CAST=None ):
    """ split multiple floating point numbers separated by space, comma, colon, semicolon.
        If additional casting is need , user can specify casting function via CAST variable
    """
    str_list = re.findall(r'[^,:;\s]+', float_str)
    flt_list = list( map( float, str_list) )  # python3 compatible
    if CAST:
        return list( map(CAST,flt_list) )
    else:
        return flt_list

####################################################################################
class BundleBlock:
    EOP_SUFFIX = ('_X0','_Y0','_Z0','_Om','_Ph','_Ka')
    PNT_SUFFIX = ('_X','_Y','_Z')
    def __init__(self, BlockFile):
        self.ReadYAML( BlockFile )
        self.AnalyzeBlock()
        self.DoAdjustment()
        pass

    ################################################################################
    def ReadYAML(self, BlockFile):
        """ read YAML and restructure its hiearchy such a way that all attributes are
         encode in ordinates string (x,y) or (x,y,z) will be splitted into list and 
         casted to floating numbers using SplitFloat() function."""
        stream = open( BlockFile, "r")
        self.BLOCK_PHOTO = yaml.load(stream, Loader=yamlordereddictloader.Loader)
        #print '==== BlockData ===',self.BLOCK_PHOTO
        # 1) check if GCP and EOP exist ..., 
        if ( 'GCP' not in self.BLOCK_PHOTO.keys() ) or self.BLOCK_PHOTO['GCP'] is None :
            self.BLOCK_PHOTO['GCP'] = dict()   #  all sections will be dict() strictly!
        else:
            for gcp_name,str_xyz in six.iteritems( self.BLOCK_PHOTO['GCP'] ):
                self.BLOCK_PHOTO['GCP'][gcp_name] = SplitFloat( str_xyz )
        if ( 'EOP' not in self.BLOCK_PHOTO.keys() ) or self.BLOCK_PHOTO['EOP'] is None :
            self.BLOCK_PHOTO['EOP'] = dict()   #  all sections will be dict() strictly!
        else:
            for img_name, EOPs in six.iteritems( self.BLOCK_PHOTO['EOP'] ):
                eop = SplitFloat( EOPs['XYZ'] ) + SplitFloat( EOPs['OPK'], CAST=math.radians )
                self.BLOCK_PHOTO['EOP'][img_name] = eop
        return

    ################################################################################
    def FindPointsOnImages(self, Points, Images):
        Point_on_Image = defaultdict( list )
        if Points is None: 
            return Point_on_Image
        else:
            for pnt in Points.keys(): 
                for img in Images.keys():
                    if pnt in Images[img].keys():
                        Point_on_Image[pnt].append( img )        
        return Point_on_Image
        
    ################################################################################
    def AnalyzeBlock(self):        
        gcps = self.BLOCK_PHOTO['GCP']
        # find GCP actually on at least one image
        self.GCP_ON_IMAGE = self.FindPointsOnImages( gcps, self.BLOCK_PHOTO['Image'] )
        if len( self.GCP_ON_IMAGE ) == 0:
            print('Number of GCPs on images : None!')
        else:
            print('Number of GCPs : {:3d}-listed  {:3d}-used'.format( len(gcps), len(self.GCP_ON_IMAGE) ) )
            ugcp = list( set(gcps.keys())  -  set(self.GCP_ON_IMAGE.keys()) )
            print('Unused GCPs    : {}'.format( str( ugcp ) ) )

        # 2) merge all poions on every images
        MeasPntList = list()
        for img, img_pnt in six.iteritems( self.BLOCK_PHOTO['Image'] ):
            #print img, img_pnt.keys()
            MeasPntList += img_pnt.keys()

        # 3) count frequency of point names
        PntFreq = Counter( MeasPntList )  # unique point name and its frequency
        self.TPNT_LIST = list( set(PntFreq.keys() ) - set(self.GCP_ON_IMAGE.keys()) )
        print( 'Number of tie-points : {}'.format(len(self.TPNT_LIST) ) )
        for pnt in self.TPNT_LIST:
            if PntFreq[pnt] <= 1:
                print('Tie Point : "{}" measured on ONE image only.'.format( pnt ) )
        return
    
    ################################################################################
    # define objective function: returns the array to be minimized
    def ObjectiveBundle(self, UnkPar):
        """ colinearity equation model """
        BP = self.BLOCK_PHOTO
        FocLen =  BP['Project']['FocLen']  # mm
        residu = list()
        for photo in BP['Image'].keys():
            for (pnt,xpyp) in six.iteritems( BP['Image'][photo] ):
                xp,yp = SplitFloat( xpyp )
                if photo in BP['EOP']:
                    EOP = BP['EOP'][photo]
                else:
                    EOP = [ UnkPar[photo+suffix].value for suffix in self.EOP_SUFFIX ]
                photo_obj = Photo(photo, FocLen, EOP[:3], EOP[3:], SEQUENCE='PIX4D' )

                if pnt in BP['GCP']:
                    XYZ = BP['GCP'][pnt]
                elif pnt in self.TPNT_LIST:
                    XYZ = [ UnkPar[pxyz].value for pxyz in self.DecorTP(pnt) ]
                else:
                    raise '***ERROR**** "{}" is neight GPS nor TP'.format(pnt)

                xp_,yp_ = photo_obj.World2Photo( XYZ )
                residu.extend( [xp_-xp,yp_-yp] )
        return residu      # model-data
    
    ######################################################################################
    def DecorTP(self, tp):
        tp_xyz = list()
        for suffix in self.PNT_SUFFIX:
            tp_xyz.append( '_' + str(tp) +suffix )
        return tp_xyz
        
    ######################################################################################
    def DoAdjustment(self):
        # create a set of Unknowns/Parameters
        Blk_Pho = self.BLOCK_PHOTO
        Params = Parameters()
        for tp_name in self.TPNT_LIST:
            tp_xyz = self.DecorTP( tp_name )
            for tp, approx in zip(tp_xyz, [1000., 1000., 100.] ):
                Params.add( tp, value=approx,  vary=True ) 

        # create unknow parameters for all EOPs
        for image_name in Blk_Pho['Image'].keys():
            if image_name not in Blk_Pho['EOP']:
                for i,approx in enumerate( [1000.,1000.,3000.,0.,0.,0.] ):
                    Params.add( image_name+ self.EOP_SUFFIX[i] , value=approx, vary=True )
                    
        # do fit, here with leastsq model
        minner = Minimizer( self.ObjectiveBundle, Params )
        self.RESULT = minner.minimize( method='leastsq')
        #self.RESULT = minner.minimize( method='cg')  # wrong solution !
    
        # write error report
        #report_fit(self.RESULT)
        return

    ################################################################################
    def Is_OMK(self,str_var):
        for kw in self.EOP_SUFFIX[3:]:
            #print('Is_OMK {}  {}'.format(kw,str_var))
            if kw in str_var[-3:]:
                return True
            else:
                continue
            return False

    ######################################################################################
    def Report(self):
        block = self.BLOCK_PHOTO
        result = self.RESULT
        print('Focal Length : {} mm.'.format( block['Project']['FocLen'] ) )
        print('Degree of Freedom : {}'.format(result.nfree))
        print('Reduced Chi-square : {:.6f} micron'.format(1000.0*result.redchi))
        ScanRes = block['Project']['ScanRes']
        print('Scanned Resolution : {} micron'.format( ScanRes ) )
        print('\n================ Adjusted Value ====================')
        for i,symbol in enumerate( list(result.params) ):
            par = result.params[symbol]
            if self.Is_OMK( par.name):
                print('{:15s} {:10.6f} deg     std.+/-{:5.3f} deg.'.format(
                    par.name, math.degrees(par.value), math.degrees( par.stderr) ) )
            else:
                print('{:15s} {:10.4f} m          std.+/-{:5.3f} m.'.format( par.name, par.value, par.stderr) )
            if (i%3)==2:   # every 3 lines, make a white line
                print('')
            
        print('\n================= Measurement Residues ====================')
        print('{:^10s}{:^10s}{:^10s}{:^10s}{:^10s}{:^10s}'.format(
              'photo','pnt','vx[mm]','vy[mm]','vx[pix]','vy[pix]'))
        
        for i, photo_name in enumerate( block['Image'].keys() ):
            for pnt, xpyp in six.iteritems( block['Image'][photo_name] ):
                vx_mm , vy_mm = result.residual[i], result.residual[i+1]
                vx_px, vy_px = 1000.*vx_mm/ScanRes, 1000.*vy_mm/ScanRes
                print('{:10s}{:10s}{:9.4f}{:9.4f}{:9.1f}{:9.1f}'.format(
                    str(photo_name), str(pnt), vx_mm , vy_mm, vx_px, vy_px ) )
                i = i+2
            print('')
        return

def SELF_TEST():
    f1 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\SinglResect_62.yml'
    f2 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\SinglResect_63.yml'
    f3 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\TwoPhotoInter_62_63.yml'
    f4 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\Block_62_63.yml'
    f5 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\Block_62_63_TiePnt.yml'
    f6 = r'D:\Phisan\Dev\PhotoComVis\DataYAML\Block_62_63_nogcp.yml'
    
    for icase, FileName in enumerate( [f1,f2,f3,f4,f5] ):
    #for icase, FileName in enumerate( [f1,f2,f3,f4] ):
    #for icase, FileName in enumerate( [f3,] ):
        print( '{} Case : {} {}'.format( 35*'%', icase+1, 35*'%' ) )
        print('Reading Block File : "{}"'.format(FileName))
        bb = BundleBlock( FileName )
        print('Project Name : {}\n'.format( bb.BLOCK_PHOTO['Project']['Name']) )

        bb.Report()
    return

################################################################################
if __name__ == "__main__":
    print(70*'=') ;  print(PROG) ; print( 70*'=' )
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--test", help="system test on variety of datasets",
                        action="store_true")
    parser.add_argument("-f", "--infile",   type=argparse.FileType('r'),
                        help="YAML block data file" )  
    args = parser.parse_args()
    
    if args.test:
        SELF_TEST()
        sys.exit(1)
        
    if args.infile:
        print('Reading Block File : "{}"'.format( args.infile.name ))
        bb = BundleBlock( args.infile.name )
        print('Project Name : {}\n'.format( bb.BLOCK_PHOTO['Project']['Name']) )
        bb.Report()

    print(70*'=')
    
