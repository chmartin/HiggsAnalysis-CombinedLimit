from HiggsAnalysis.CombinedLimit.PhysicsModel import *

### This is the base python class to study the Higgs width

class Higgswidth(PhysicsModel):
    def __init__(self):
        self.mHRange = []
        self.GGsmfixed = False
        self.is2l2nu = False
        self.poiMap = []
        self.pois = {}
        self.verbose = False
    def setModelBuilder(self, modelBuilder):
        PhysicsModel.setModelBuilder(self,modelBuilder)
        self.modelBuilder.doModelBOnly = False

    def getYieldScale(self,bin,process):
        if self.is2l2nu:
            name = "%s_%s_func"%(process,bin)
            if process == "ggH_s":
                self.modelBuilder.factory_("expr::%s(\"@0*@1-sqrt(@0*@1)\",CMS_zz4l_mu,CMS_zz4l_GGsm)"%(name))
                return name
            elif process == "ggH_b":
                self.modelBuilder.factory_("expr::%s(\"1-sqrt(@0*@1)\",CMS_zz4l_mu,CMS_zz4l_GGsm)"%(name))
                return name
            elif process == "ggH_sbi":
                self.modelBuilder.factory_("expr::%s(\"sqrt(@0*@1)\",CMS_zz4l_mu,CMS_zz4l_GGsm)"%(name))
                return name
            else:
                return 1
        else:
            return 1
            
    def setPhysicsOptions(self,physOptions):
        for po in physOptions:
            if po == "GGsmfixed":
                print "Will fix CMS_zz4l_GGsm to 1 and float mu"
                self.GGsmfixed = True
            if po == "is2l2nu":
                print "Will consider cards in 2l2nu style (separated S, B, S+B+I)"
                self.is2l2nu = True
            
    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        if self.is2l2nu:
            self.modelBuilder.doVar("CMS_zz4l_GGsm[1.,0.,50.]")
            self.modelBuilder.doVar("CMS_zz4l_mu[1.,0.,4]")
            if self.GGsmfixed:
                print "Fixing CMS_zz4l_GGsm"
                poi = "CMS_zz4l_mu"
            else:
                poi = "CMS_zz4l_GGsm"
        if self.GGsmfixed:
            self.modelBuilder.out.var("CMS_zz4l_GGsm[1.]")
            self.modelBuilder.out.var("CMS_zz4l_mu")
            print "Fixing CMS_zz4l_GGsm"
            poi = "CMS_zz4l_mu"
        else:
            self.modelBuilder.out.var("CMS_zz4l_GGsm[25.]")
            self.modelBuilder.out.var("CMS_zz4l_mu[1.]")
            poi = "CMS_zz4l_GGsm"
        
        self.modelBuilder.doSet("POI",poi)
        
higgswidth = Higgswidth()
