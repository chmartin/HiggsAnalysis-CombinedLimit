from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the SpinZero structure

class SpinZeroHiggs(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.muAsPOI = False
        self.muFloating = False
        self.poiMap = []
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        "Split in production and decay, and call getHiggsSignalYieldScale; return 1 for backgrounds "
        
        if not self.DC.isSignal[process]: return 1

        if self.muFloating:
            print "Process ", process, "will get norm r"
            return "r"

        else:
            print "Process ", process, " will get norm 1"
            return 1

    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "muAsPOI":
                print "Will consider the signal strength as a parameter of interest"
                self.muAsPOI = True
                self.muFloating = True
            elif po == "muFloating":
                print "Will consider the signal strength as a floating parameter (as a paraeter of interest if --PO muAsPOI is specified, as a nuisance otherwise)"
                self.muFloating = True
            else:
                print "Option not understood"
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("CMS_zz4l_fg4[0,0,1]")
        poi = "CMS_zz4l_fg4"

        if self.muFloating:
            self.modelBuilder.doVar("r[1,0,4]")
            if self.muAsPOI:
                print 'Treating r as a POI'
                poi += ",r"
        
        self.modelBuilder.doSet("POI",poi)
        
spinZeroHiggs = SpinZeroHiggs()
