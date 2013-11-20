from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggs(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.muAsPOI = False
        self.muFloating = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False
        
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "muAsPOI":
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            if po == "muFloating":
                print "Will consider the signal strength as a floating parameter (as a paraeter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True
            if po.startswith("verbose"):
                self.verbose = True
            if po.startswith("map="):
                self.muFloating = True
                (maplist,poi) = po.replace("map=","").split(":")
                maps = maplist.split(",")
                poiname = re.sub("\[.*","",poi)
                if poiname not in self.pois:
                    if self.verbose: print "Will create a var ",poiname," with factory ",poi
                    self.pois[poiname] = poi
                if self.verbose: print "Mapping ",poiname," to ",maps," patterns"
                self.poiMap.append((poiname,maps))
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("CMS_zz4l_fg4[0,0,1]")
        poi = "CMS_zz4l_fg4"

        if self.muFloating:
            if self.pois:
                for ph,pf in self.pois.items():
                    self.modelBuilder.doVar(pf)
                    if self.muAsPOI:
                        print 'Treating %(pn)s as a POI' %locals()
                        poi += ','+pn
            else:
                self.modelBuilder.doVar("r[1,0,4]")
                if self.muAsPOI:
                    print 'Treating r as a POI'
                    poi += ",r"
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
