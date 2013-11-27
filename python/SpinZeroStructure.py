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

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        
        if not self.DC.isSignal[process]: return 1   

        if self.pois:
            target = "%(bin)s/%(process)s" % locals()
            scale = 1
            for p,l in self.poiMap:
                for el in l:
                    if re.martch(el, target):
                        scale = p + self.my_norm
            print "Will scale ", target, " by ",scale
            return scale
        else:
            print "Process ",process, " will get norm ", self.my_norm
            return self.my_norm

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "muAsPOI":
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            if po == "muFloating":
                print "Will consider the signal strength as a floating parameter (as a paraeter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True
            if po.startswith("map="):
                self.muFloating = True
                (maplist,poi) = po.replace("map=","").split(":")
                maps = maplist.split(",")
                poiname = re.sub("\[.*","", poi)
                if poiname not in self.pois:
                    if self.verbose: print "Will create a var ",poiname," with factory ",poi
                    self.pois[poiname] = poi
                if self.verbose:  print "Mapping ",poiname," to ",maps," patterns"
                self.poiMap.append((poiname, maps))
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("CMS_zz4l_fg4[0,0,1]")
        poi = "CMS_zz4l_fg4"

        if self.muFloating:
            if self.pois:
                for pn,pf in self.pois.items():
                    self.modelBuilder.doVar(pf)
                    if self.muAsPOI:
                        print 'Treating %(pn)s as a POI' % locals()
                        poi += ','+pn
                    self.modelBuilder.factory_('%(pn)s_tiems_temp("@0*(1.-0.826667*@1)",%(pn)s,CMS_zz4l_fg4)' % locals())
                self.my_norm = '_times_temp'
            else:
                self.modelBuilder.doVar("r[1,0,4]")
                self.modelBuilder.factory_('expr::r_times_temp("@0*(1.-0.826667*@1)",r,CMS_zz4l_fg4)')
                self.my_norm = 'r_times_temp'
                if self.muAsPOI:
                    print 'Treating r as a POI'
                    poi += ",r"
        else:
            self.modelBuilder.factory_('expr::One_times_temp(\"(1.-0.826667*@0)\",CMS_zz4l_fg4)')
            self.my_norm = 'One_times_temp'
        
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
