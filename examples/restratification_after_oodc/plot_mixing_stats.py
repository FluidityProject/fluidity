import fluidity_tools
import pylab
from matplotlib.pyplot import *
subplot(121)


labels=['0<= T <0.1', '0.1<= T <0.2', '0.2<= T <0.3', '0.3<= T <0.4', '0.4<= T <0.5', '0.5<= T <0.6', '0.6<= T <0.7', '0.7<= T <0.8', '0.8<= T <0.9', '0.9<= T < 1.0', '1.0<= T < 1.1', '1.1<= T <1.2', '1.2<= T <1.3', '1.3<= T <1.4', '1.4<= T <1.5', '1.5<= T <1.6', '1.6<= T <1.8']
lines=['-','--',  ':',  '-','--',  ':','-','--',  ':', '-','--',  ':','-','--',  ':','-','--',  ':',]
colours=['black', 'black','black','blue','blue','blue' , 'purple',   'purple',  'purple', 'red' ,  'red' , 'red' , 'green',  'green', 'green', 'grey', 'grey','grey']

for i in range (0,(len(labels)-1)):
  bins=fluidity_tools.stat_parser('restratification_after_oodc.stat')['BoussinesqFluid']['Temperature']['mixing_bins%mixingf'][i]
  t=fluidity_tools.stat_parser('restratification_after_oodc.stat')['ElapsedTime']['value']
  bins=fluidity_tools.stat_parser('restratification_after_oodc.stat')['BoussinesqFluid']['Temperature']['mixing_bins%mixingf'][i]
  if (i==16):
    bins2=fluidity_tools.stat_parser('restratification_after_oodc.stat')['BoussinesqFluid']['Temperature']['mixing_bins%mixingf'][i+1]
    pylab.plot(t/3600/24,(bins+bins2)/(1.9634e14),label=labels[i],linewidth=3.0,linestyle=lines[i], color=colours[i])
  else:
    pylab.plot(t/3600/24,bins/(1.9634e14),label=labels[i],linewidth=3.0,linestyle=lines[i], color=colours[i])
  #pylab.legend(ncol=2,markerscale=0.01,columnspacing=0.1,prop=12)
legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
# setting font size
fontsize=14
ax = pylab.gca()
for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
for label in [ax.get_xaxis().get_label(), ax.get_yaxis().get_label()]:
    label.set_fontsize(fontsize)
for label in ax.get_axes().get_xticklabels(): label.set_fontsize(fontsize)
for label in ax.get_axes().get_yticklabels(): label.set_fontsize(fontsize)
pylab.rc("axes", linewidth=3.0)
label = pylab.gca().get_yaxis().get_offset_text()
label.set_fontsize(fontsize)
label = pylab.gca().get_xaxis().get_offset_text()
label.set_fontsize(fontsize)
#ax.yaxis.major.formatter.set_powerlimits((-2,2))
#ax.xaxis.major.formatter.set_powerlimits((-2,2))
pylab.xlabel('Time (days)',fontsize=fontsize)
pylab.ylabel('Proportion within the given temperature range',fontsize=fontsize)
pylab.savefig('mixing_stats.png')
#pylab.show()


