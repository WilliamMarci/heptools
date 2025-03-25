lumi_units = {
    'nb-1':1e-3,
    'pb-1':1,
    'fb-1':1e3,
    'ab-1':1e6}
xsec_units = {
    'nb':1e3,
    'pb':1,
    'fb':1e-3,
    'ab':1e-6}
######################
class HEPExp:
    def __init__(self):
        self.lumi_num=0
        self.lumi_unit=0
        self.xsec_num=0
        self.xsec_unit=0
        self.pa_exp=0
        self.event=0
        self.exps=0
        self.title='default'
    def set_title(self,title):
        self.title=title
    def set_lumi(self,lumi_num,lumi_unit):
        self.lumi_num=lumi_num
        self.lumi_unit=lumi_unit
    def set_xsec(self,xsec_num,xsec_unit):
        self.xsec_num=xsec_num
        self.xsec_unit=xsec_unit
    def set_LHC(self, run='Run2'):
        if run=='Run2':
            self.set_lumi(140, 'fb-1')
        elif run=='Run3':
            self.set_lumi(300, 'fb-1')
        elif run=='HL-LHC':
            self.set_lumi(3000, 'fb-1')
        else:
            raise ValueError('Please input Run2, Run3, or HL-LHC')
    def set_simulation(self,pa_exp):
        self.pa_exp=pa_exp
    def set_crosssection(self,xsec_num,xsec_unit):
        self.xsec_num=xsec_num
        self.xsec_unit=xsec_unit
    def calculate(self):
        lumi=lumi_units[self.lumi_unit]*self.lumi_num
        xsec=xsec_units[self.xsec_unit]*self.xsec_num
        event=xsec*lumi
        exps=self.pa_exp/event
        self.event=event
        self.exps=exps
    def info(self):
        self.calculate()
        if self.title=='default':
            print("--------------------------------")
        else:
            print("-----",self.title,"-----")
        print("Lumi: ",self.lumi_num, self.lumi_unit)
        print("Cross-Section: ",self.xsec_num, self.xsec_unit)
        print("Events: ",self.event)
        print("Experiment: ",self.exps)
        print("--------------------------------")
    def get_event(self):
        self.calculate()
        return self.event
    def get_exp(self):
        self.calculate()
        return self.exps

