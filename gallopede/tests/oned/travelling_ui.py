#!/usr/bin/env pythonw
import matplotlib
matplotlib.interactive(True)
#matplotlib.use('WX')
import wx
from string import *
from scipy import *
import trav_eqns_3layer
import signal
import trav_eqns_3layer_rigid
import interrupt 
class Form(wx.Panel):
    ''' The Form class is a wx.Panel that creates a bunch of controls
        and handlers for callbacks. Doing the layout of the controls is 
        the responsibility of subclasses (by means of the doLayout()
        method). '''

    def __init__(self,type,layer_no, *args, **kwargs):
        self.type=type
        self.layer_no=layer_no
        wx.Panel.__init__(self,*args, **kwargs)
        self.createControls(self,*args, **kwargs)
        self.bindEvents(self,*args, **kwargs)
        self.doLayout(*args,**kwargs)
        
    def createControls(self,*args, **kwargs):
        self.runButton = wx.Button(self, label="Run")
        self.plotButton = wx.Button(self, label="Plots")
        self.bifButton = wx.Button(self, label="Bifurcation")
        self.potButton = wx.Button(self, label="Potential")
        self.exitButton = wx.Button(self, label="Exit")
        if self.type=="layer":
            self.layerLabel = wx.StaticText(self, label="Number of layers")
            self.layerTextCtrl = wx.TextCtrl(self, value="%i"%self.layer_no)
            self.speedLabel = wx.StaticText(self, label="Wave group speed")
            self.speedTextCtrl = wx.TextCtrl(self, value="1.3")
            self.t0Label = wx.StaticText(self, label="x_0")
            self.t0TextCtrl = wx.TextCtrl(self, value="-50.0")
            self.tfLabel = wx.StaticText(self, label="x_f")
            self.tfTextCtrl = wx.TextCtrl(self, value="2000.0")
            self.dtLabel = wx.StaticText(self, label="Delta x")
            self.dtTextCtrl = wx.TextCtrl(self, value="1.0")
            self.ptLabel = wx.StaticText(self, label="Output x")
            self.ptTextCtrl = wx.TextCtrl(self, value="5.0")
            self.epsilonLabel = wx.StaticText(self, label="Epsilon")
            self.epsilonTextCtrl = wx.TextCtrl(self, value="1.0e-5")
            self.thLabel= wx.StaticText(self, label="Theta")
            self.thTextCtrl= wx.TextCtrl(self, value="0.0")
            self.rigidLabel = wx.StaticText(self, label="Top surface")
            self.rigidButton= wx.CheckBox(self, -1, label="Rigid Lid")
            self.shootLabel = wx.StaticText(self, label="Shooting")
            self.shootButton= wx.CheckBox(self, -1, label="On")
            
        if self.type=="rhodata":
            self.rhoLabel=[]
            self.rhoTextCtrl=[]
            for i in range(10):
                self.rhoLabel.append(wx.StaticText(self, 
                                                   label=('rho_%i'%(i+1))))
                self.rhoTextCtrl.append(wx.TextCtrl(self,
                                                    id=(1001+i),
                                                    value=('%f'%(1023+3*i))))
            for i in range(self.layer_no,10):
                self.rhoTextCtrl[i].Enable(False)

        if self.type=="ddata":
            self.dLabel=[]
            self.dTextCtrl=[]
            self.dLabel.append(wx.StaticText(self, 
                                             label=('D_%i'%(1))))
            self.dTextCtrl.append(wx.TextCtrl(self,
                                              id=(1001),
                                              value=('%f'%(50.0))))
            for i in range(1,10):
                self.dLabel.append(wx.StaticText(self, 
                                                   label=('D_%i'%(i+1))))
                self.dTextCtrl.append(wx.TextCtrl(self,
                                                    id=(1001+i),
                                                    value=('%f'%(350.0))))
            for i in range(self.layer_no,10):
                self.dTextCtrl[i].Enable(False)



    def bindEvents(self,*args, **kwargs):
        pass

    def doLayout(self,*args,**kwargs):
        ''' Layout the controls that were created by createControls(). 
            Form.doLayout() will raise a NotImplementedError because it 
            is the responsibility of subclasses to layout the controls. '''
        raise NotImplementedError

    # Callback methods:

    def onSave(self, event):
         self.__log('Save_data')


    def onLayerEntered(self, event):
        self.__log('User entered name: %s'%event.GetString())
        if len(event.GetString())>0:
            self.layer_no=atoi(event.GetString())

    def onSpeedEntered(self, event):
        self.__log('User entered name: %s'%event.GetString())
        if len(event.GetString())>0:
            self.speed=atof(event.GetString())

    # Helper method(s):

    def __log(self, message):
        ''' Private method to append a string to the logger text
            control. '''
        print message

class FormWithSizer1(Form):
    def doLayout(self,type):
        ''' Layout the controls by means of sizers. '''

        # A horizontal BoxSizer will contain the GridSizer (on the left)
        # and the logger text control (on the right):
        boxSizerbut = wx.FlexGridSizer(rows=1, cols=5, vgap=10, hgap=10)
        # A GridSizer will contain the other controls:
        gridSizer = wx.FlexGridSizer(rows=10, cols=2, vgap=10, hgap=10)
        boxSizer = wx.FlexGridSizer(rows=2, cols=1, vgap=10, hgap=10)

        # Prepare some reusable arguments for calling sizer.Add():
        expandOption = dict(flag=wx.EXPAND)
        noOptions = dict()
        emptySpace = ((0, 0), noOptions)

        for control, options in \
                [(self.runButton, dict(flag=wx.ALIGN_CENTER)),
                 (self.plotButton, dict(flag=wx.ALIGN_CENTER)),
                 (self.bifButton, dict(flag=wx.ALIGN_CENTER)),
                  (self.potButton, dict(flag=wx.ALIGN_CENTER)),
                 (self.exitButton, dict(flag=wx.ALIGN_CENTER))]:
              boxSizerbut.Add(control, **options)

        # Add the controls to the sizers:
        if self.type == 'layer':
            for control, options in \
                    [(self.layerLabel, noOptions),
                     (self.layerTextCtrl, expandOption),
                     (self.speedLabel, noOptions),
                     (self.speedTextCtrl, expandOption),
                     (self.t0Label, noOptions),
                     (self.t0TextCtrl, expandOption),
                     (self.tfLabel, noOptions),
                     (self.tfTextCtrl, expandOption),
                     (self.dtLabel, noOptions),
                     (self.dtTextCtrl, expandOption),
                     (self.ptLabel, noOptions),
                     (self.ptTextCtrl, expandOption),
                     (self.epsilonLabel, noOptions),
                     (self.epsilonTextCtrl, expandOption),
                     (self.thLabel, noOptions),
                     (self.thTextCtrl, expandOption),
                     (self.rigidLabel, noOptions),
                     (self.rigidButton, dict(flag=wx.ALIGN_CENTER)),
                     (self.shootLabel, noOptions),
                     (self.shootButton, dict(flag=wx.ALIGN_CENTER))]:
                gridSizer.Add(control, **options)
                
        if self.type == 'rhodata':
            mlist=[]
            for i in range(10):
                mlist.append((self.rhoLabel[i], noOptions))
                mlist.append((self.rhoTextCtrl[i], expandOption))
            for control, options in mlist:
                gridSizer.Add(control, **options)

        if self.type == 'ddata':
            mlist=[]
            for i in range(10):
                mlist.append((self.dLabel[i], noOptions))
                mlist.append((self.dTextCtrl[i], expandOption))
            for control, options in mlist:
                gridSizer.Add(control, **options)

        for control, options in \
                [(gridSizer, dict(border=5, flag=wx.ALL)),
                 (boxSizerbut, dict(border=5, flag=wx.ALL))]:
            boxSizer.Add(control, **options)

        self.SetSizerAndFit(boxSizer)

class FrameWithForms(wx.Frame):
    def __init__(self, *args, **kwargs):
        super(FrameWithForms, self).__init__(*args, **kwargs)
        self.notebook = wx.Notebook(self)
        self.layer_no=2
        self.th=zeros(19)
        self.rho=zeros(10)
        self.d=zeros(10)
        self.forms=[]
        self.D=zeros([0,0])
        self.q=zeros([0,0])
        self.t=zeros([0,0])
        self.forms.append(FormWithSizer1('layer',
                                    self.layer_no,self.notebook,))
        self.forms.append( FormWithSizer1('rhodata',
                                    self.layer_no,self.notebook,))
        self.forms.append(FormWithSizer1('ddata',self.layer_no,self.notebook,))
        self.notebook.AddPage(self.forms[0], 'Globals')
        self.notebook.AddPage(self.forms[1], 'Densities')
        self.notebook.AddPage(self.forms[2], 'Heights')

        self.c=atof(self.forms[0].speedTextCtrl.GetValue())
        self.t0=atof(self.forms[0].t0TextCtrl.GetValue())
        self.tf=atof(self.forms[0].tfTextCtrl.GetValue())
        self.dt=atof(self.forms[0].dtTextCtrl.GetValue())
        self.pt=atof(self.forms[0].ptTextCtrl.GetValue())
        self.epsilon=atof(self.forms[0].epsilonTextCtrl.GetValue())
        for i in range(10):
            self.rho[i]=atof(self.forms[1].rhoTextCtrl[i].GetValue())
            self.d[i]=atof(self.forms[2].dTextCtrl[i].GetValue())
        self.bindEvents()
        # We just set the frame to the right size manually. This is feasible
        # for the frame since the frame contains just one component. If the
        # frame had contained more than one component, we would use sizers
        # of course, as demonstrated in the FormWithSizer class above.
        self.SetClientSize(self.notebook.GetBestSize())
        
    def bindEvents(self):
        for i in range(3):
            self.Bind(wx.EVT_BUTTON,self.onExit,self.forms[i].exitButton)
            self.Bind(wx.EVT_BUTTON,self.onRun,self.forms[i].runButton)
            self.Bind(wx.EVT_BUTTON,self.onPlot,self.forms[i].plotButton)
            self.Bind(wx.EVT_BUTTON,self.onBif,self.forms[i].bifButton)
            self.Bind(wx.EVT_BUTTON,self.onPot,self.forms[i].potButton)
        self.Bind(wx.EVT_TEXT,self.onLayerEntered,self.forms[0].layerTextCtrl)

    def onExit(self,e):
        print "Hello"
        self.Close(True)
    
    def onBif(self,e):
        print ('Make bifuctation diagram')
        self.reloadModules()

        self.c=atof(self.forms[0].speedTextCtrl.GetValue())
        for i in range(self.layer_no):
            self.rho[i]=atof(self.forms[1].rhoTextCtrl[i].GetValue())
            self.d[i]=atof(self.forms[2].dTextCtrl[i].GetValue())

        print "passed"
        print "N=", self.layer_no
        print "c=", self.c
        print "rho=", self.rho[0:self.layer_no]
        print "D=", self.d[0:self.layer_no]

        if self.forms[0].rigidButton.GetValue():
            trav_eqns_3layer_rigid.bifurcationDiagram(self.c,
                                                      self.rho[0:self.layer_no],
                                                      self.d[0:self.layer_no])
        else :
            trav_eqns_3layer.bifurcationDiagram(self.c,
                                                self.rho[0:self.layer_no],
                                                self.d[0:self.layer_no])
            
    def onPot(self,e):
        print ('Make potential diagrams')
        self.reloadModules()
        self.c=atof(self.forms[0].speedTextCtrl.GetValue())
        for i in range(self.layer_no):
            self.rho[i]=atof(self.forms[1].rhoTextCtrl[i].GetValue())
            self.d[i]=atof(self.forms[2].dTextCtrl[i].GetValue())

        print "passed"
        print "N=", self.layer_no
        print "c=", self.c
        print "rho=", self.rho[0:self.layer_no]
        print "D=", self.d[0:self.layer_no]

        
        trav_eqns_3layer.potMap(self.c,
                                            self.rho[0:self.layer_no],
                                            self.d[0:self.layer_no],False)

    def onRun(self,e):
        signal.signal(signal.SIGINT,goodbyeCruelWorld)
        print signal.getsignal(signal.SIGINT)
        self.reloadModules()
        print ('Run shooting code')
        self.c=atof(self.forms[0].speedTextCtrl.GetValue())
        self.t0=atof(self.forms[0].t0TextCtrl.GetValue())
        self.tf=atof(self.forms[0].tfTextCtrl.GetValue())
        self.dt=atof(self.forms[0].dtTextCtrl.GetValue())
        self.pt=atof(self.forms[0].ptTextCtrl.GetValue())
        self.epsilon=atof(self.forms[0].epsilonTextCtrl.GetValue())
        self.th[0]=atof(self.forms[0].thTextCtrl.GetValue())
        for i in range(self.layer_no):
            self.rho[i]=atof(self.forms[1].rhoTextCtrl[i].GetValue())
            self.d[i]=atof(self.forms[2].dTextCtrl[i].GetValue())

        print self.c, self.t0, self.tf, self.dt, self.pt, self.epsilon

        if self.forms[0].rigidButton.GetValue():
            trav_eqns_3layer_rigid.trav_maker(self.c,
                                              self.t0,
                                              self.tf,
                                              self.dt,
                                              self.pt,
                                              self.epsilon,
                                              self.rho[0:self.layer_no],
                                              self.d[0:self.layer_no],
                                              self.th[0:2*self.layer_no-1])
        else :
            if self.forms[0].shootButton.GetValue():
                self.th=trav_eqns_3layer.shooting(self.c,
                                        self.t0,
                                        self.tf,
                                        self.dt,
                                        self.pt,
                                        self.epsilon,
                                        self.rho[0:self.layer_no],
                                        self.d[0:self.layer_no],
                                                self.th[0:2*self.layer_no-1])
            (self.D,self.q,self.t)=trav_eqns_3layer.trav_maker(self.c,
                                        self.t0,
                                        self.tf,
                                        self.dt,
                                        self.pt,
                                        self.epsilon,
                                        self.rho[0:self.layer_no],
                                        self.d[0:self.layer_no],
                                        self.th[0:2*self.layer_no-1])

        return
                      
    def onPlot(self,e):
        self.reloadModules()
        N=self.layer_no
        print ('Run shooting code')
        print self.c, self.t0, self.tf, self.dt, self.pt, self.epsilon
        trav_eqns_3layer.pic_maker(self.D,self.q,self.t,self.c,
                                   self.d[0:N],self.rho[0:N])
        return

    def reloadModules(self):
        reload(trav_eqns_3layer)
        reload(trav_eqns_3layer_rigid)
        return

    def onLayerEntered(self,e):
        print ('User entered layer no. %s'%e.GetString())
        if len(e.GetString())>0:
            self.layer_no=atoi(e.GetString())
#            self.form3.createControls('rhodata',self.layer_no,self.notebook)
#            self.form3.doLayout('rhodata')
            for i in range(self.layer_no):
                self.forms[1].rhoTextCtrl[i].Enable(True)
                self.forms[2].dTextCtrl[i].Enable(True)
            for i in range(self.layer_no,10):
                self.forms[1].rhoTextCtrl[i].Enable(False)
                self.forms[2].dTextCtrl[i].Enable(False)
            
def goodbyeCruelWorld(signum,frame):
    raise interrupt.MyError(signum)


signal.signal(signal.SIGINT,goodbyeCruelWorld)
if __name__ == '__main__':
   
    print signal.getsignal(signal.SIGINT)
    app = wx.App(0)
    frame = FrameWithForms(None, title='Travelling Wave UI')
    frame.Show()
    app.MainLoop()
