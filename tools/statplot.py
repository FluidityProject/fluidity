#!/usr/bin/env python3

import gi
gi.require_version('Gdk', '3.0')
gi.require_version('Gtk', '3.0')
from gi.repository import Gdk, Gtk  # noqa: E402

import argparse  # noqa: E402
from lxml import etree  # noqa: E402
from matplotlib import rc  # noqa: E402
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)  # noqa: E402
from matplotlib.backends.backend_gtk3 import (
    NavigationToolbar2GTK3 as NavigationToolbar)  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import matplotlib.ticker as tck  # noqa: E402
import numpy  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402

rc('text', usetex=False)

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='''GUI for Fluidity .stat outputs

When the window is displayed, press the following keys for additional options:
- x, y -> Change the scale of the x, y axis (automatically chooses between
          linear, log and symlog).
- l    -> Change the representation of the data (solid line or markers).
- r    -> Reload the .stat file. Only relevant for simulations that are still
          running.
- a, o -> Grab focus onto the abscissa, ordinate text boxes.
- Esc  -> Release focus on text boxes.
- q    -> Exit the GUI.''')
parser.add_argument('statfile', nargs='+', help='Path(s) to the .stat file(s)')
args = parser.parse_args(sys.argv[1:])


class StatplotWindow(Gtk.Window):
    def __init__(self, statfile):
        Gtk.Window.__init__(self)
        self.connect('key-release-event', self.KeyPressed)
        self.set_border_width(8)
        self.set_default_size(1600, 900)
        self.set_position(Gtk.WindowPosition.CENTER)
        self.set_title(statfile[-1])
        self.statfile = statfile
        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        vbox.set_homogeneous(False)
        self.add(vbox)
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
        hbox.set_homogeneous(True)
        vbox.pack_end(hbox, False, False, 0)
        self.entries, self.values = self.ReadData()
        store = Gtk.ListStore.new([str])
        for entry in self.entries:
            store.append([entry])
        self.xCombo = Gtk.ComboBox.new_with_model_and_entry(store)
        self.yCombo = Gtk.ComboBox.new_with_model_and_entry(store)
        self.ComboParams(self.xCombo)
        self.ComboParams(self.yCombo)
        self.InitCombo('load')
        self.InitCompletion(self.xCombo)
        self.InitCompletion(self.yCombo)
        self.xCon = self.xCombo.connect('changed', self.ComboChanged)
        self.xCombo.connect('key-release-event', self.ReleaseFocus)
        self.yCon = self.yCombo.connect('changed', self.ComboChanged)
        self.yCombo.connect('key-release-event', self.ReleaseFocus)
        hbox.pack_start(self.xCombo, True, True, 0)
        hbox.pack_start(self.yCombo, True, True, 0)
        self.fig = Figure(figsize=(14, 7))
        self.ax = self.fig.gca()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_can_focus(True)
        vbox.pack_start(self.canvas, True, True, 0)
        self.PlotType = 'line'
        self.PlotData('create', self.PlotType, 'linear', 'linear')
        self.toolbar = NavigationToolbar(self.canvas, self)
        vbox.pack_start(self.toolbar, False, False, 0)

    def ComboChanged(self, combo_box):
        activeText = combo_box.get_child().get_text()
        if activeText in self.entries:
            combo_box.set_active(self.entries.index(activeText))
            self.PlotData('update', self.PlotType, 'linear',
                          self.ax.get_yscale())
            self.canvas.grab_focus()

    def ComboParams(self, comboBox):
        comboBox.set_entry_text_column(0)
        comboBox.set_margin_left(100)
        comboBox.set_margin_right(100)
        comboBox.set_wrap_width(3)

    def FormatAxis(self, ax2fmt, scale):
        if ax2fmt == 'x':
            setScale = self.ax.set_xscale
            curAxis = self.ax.xaxis
            curData = self.xData
        elif ax2fmt == 'y':
            setScale = self.ax.set_yscale
            curAxis = self.ax.yaxis
            curData = self.yData
        setScale(scale)
        if scale == 'linear':
            self.ax.ticklabel_format(style='sci', axis=ax2fmt,
                                     scilimits=(0, 0), useMathText=True)
            curAxis.set_minor_locator(tck.AutoMinorLocator())
            self.ax.relim()
            self.ax.autoscale(True, ax2fmt, None)
        elif scale == 'log':
            curAxis.set_minor_locator(
                tck.LogLocator(subs=numpy.arange(2, 10)))
            logFmt = tck.LogFormatterSciNotation(base=10,
                                                 labelOnlyBase=False,
                                                 minor_thresholds=(4, 1))
            curAxis.set_minor_formatter(logFmt)
            self.ax.relim()
            self.ax.autoscale(True, ax2fmt, None)
        elif scale == 'symlog':
            axMin = min(abs(curData[curData != 0]))
            axMax = max(abs(curData))
            axRange = numpy.log10(axMax / axMin)
            if ax2fmt == 'x':
                setScale(
                    'symlog', basex=10, subsx=numpy.arange(2, 10),
                    linthreshx=axMin * 10 ** (axRange / 2))
            elif ax2fmt == 'y':
                setScale(
                    'symlog', basey=10, subsy=numpy.arange(2, 10),
                    linthreshy=axMin * 10 ** (axRange / 2))
            # Thomas Duvernay, 06/01/19
            # There seems to be a bug with the labelling of the 0 tick
            # when a 'symlog' is used as an axis scale. It looks like
            # it is considered as a minor tick.
            symLogLoc = tck.SymmetricalLogLocator(
                subs=numpy.arange(2, 10),
                linthresh=axMin * 10 ** (axRange / 2), base=10)
            curAxis.set_minor_locator(symLogLoc)
            logFmt = tck.LogFormatterSciNotation(
                base=10, labelOnlyBase=False, minor_thresholds=(4, 1),
                linthresh=axMin * 10 ** (axRange / 2))
            curAxis.set_minor_formatter(logFmt)
        self.ax.set_xlabel(self.xCombo.get_child().get_text(),
                           fontweight='bold', fontsize=20)
        self.ax.set_ylabel(self.yCombo.get_child().get_text(),
                           fontweight='bold', fontsize=20)
        self.ax.tick_params(which='major', length=7, labelsize=16, width=2)
        self.ax.tick_params(which='minor', length=4, labelsize=10, width=2,
                            colors='xkcd:scarlet', labelrotation=45)
        self.ax.xaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                            color='xkcd:black')
        self.ax.yaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                            color='xkcd:black')
        self.fig.set_tight_layout(True)

    def InitCombo(self, action, prevIterX=None, prevIterY=None):
        if action == 'load':
            if 'ElapsedTime' in self.entries:
                iterX = self.xCombo.get_model().get_iter(
                    self.entries.index('ElapsedTime'))
                self.xCombo.set_active_iter(iterX)
                self.yCombo.set_active(0)
            else:
                self.xCombo.set_active(0)
                self.yCombo.set_active(1)
        elif action == 'reload':
            self.xCombo.set_active_iter(prevIterX)
            self.yCombo.set_active_iter(prevIterY)

    def InitCompletion(self, comboBox):
        completion = Gtk.EntryCompletion.new()
        completion.set_text_column(0)
        completion.set_inline_completion(True)
        completion.set_inline_selection(False)
        completion.set_model(comboBox.get_model())
        comboBox.get_child().set_completion(completion)

    def KeyPressed(self, widget, event):
        key = event.string
        if (self.xCombo.get_child().has_focus()
                or self.yCombo.get_child().has_focus()):
            pass
        elif key == 'r':
            self.RefreshData()
            self.PlotData('update', self.PlotType, self.ax.get_xscale(),
                          self.ax.get_yscale())
        elif key == 'q':
            self.destroy()
        elif key == 'x' or key == 'y':
            if key == 'x':
                self.get_scale = self.ax.get_xscale
                curData = self.xData
            elif key == 'y':
                self.get_scale = self.ax.get_yscale
                curData = self.yData
            if (self.get_scale() == 'linear'
                    and (curData == 0).all()):
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the ' + key + ' axis '
                              + 'is null.', stacklevel=2)
                scale = 'linear'
            elif (self.get_scale() == 'linear'
                    and max(abs(curData)) / min(abs(curData)) < 10):
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the ' + key + ' axis '
                              + 'has a range of variation smaller than one '
                              + 'order of magnitude.', stacklevel=2)
                scale = 'linear'
            elif self.get_scale() == 'linear':
                axMin = min(abs(curData[curData != 0]))
                axMax = max(abs(curData))
                if axMin != axMax and min(curData) < 0:
                    scale = 'symlog'
                elif axMin != axMax:
                    scale = 'log'
                else:
                    warnings.warn('Change to logarithmic scale denied: the '
                                  + 'selected variable for the ' + key + ' '
                                  + 'axis is a constant.', stacklevel=2)
                    scale = 'linear'
            elif self.get_scale() in ['log', 'symlog']:
                scale = 'linear'
            self.FormatAxis(key, scale)
            self.fig.canvas.draw()
        elif key == 'l':
            self.PlotType = 'line' if self.PlotType == 'marker' else 'marker'
            self.PlotData('update', self.PlotType, self.ax.get_xscale(),
                          self.ax.get_yscale())
        elif key == 'a':
            self.xCombo.grab_focus()
        elif key == 'o':
            self.yCombo.grab_focus()

    def PlotData(self, action, type, xscale, yscale):
        self.xData = self.values[..., self.xCombo.get_active()]
        self.yData = self.values[..., self.yCombo.get_active()]
        if (numpy.unique(self.xData).size == 1
                and numpy.unique(self.yData).size == 1):
            type = 'marker'
        if action == 'create':
            self.statplot, = self.ax.plot(self.xData, self.yData, linewidth=2,
                                          color='xkcd:light purple')
        elif action == 'update':
            self.statplot.set_xdata(self.xData)
            self.statplot.set_ydata(self.yData)
            if type == 'line':
                self.statplot.set_color('xkcd:light purple')
                self.statplot.set_linestyle('solid')
                self.statplot.set_linewidth(2)
                self.statplot.set_marker('None')
            elif type == 'marker':
                self.statplot.set_linestyle('None')
                self.statplot.set_marker('d')
                self.statplot.set_markeredgecolor('xkcd:black')
                self.statplot.set_markerfacecolor('xkcd:light purple')
                self.statplot.set_markersize(7)
                self.statplot.set_markeredgewidth(0.3)
            self.toolbar.update()
        self.FormatAxis('x', xscale)
        self.FormatAxis('y', yscale)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def ReadData(self):
        def GatherEntries(stat2read):
            print('Reading ' + stat2read)
            entries = []
            with open(stat2read, 'r') as fid:
                for num, line in enumerate(fid):
                    if line.startswith('<field'):
                        info = etree.fromstring(line)
                        if info.attrib['statistic'] == 'value':
                            entries.append(info.attrib['name'])
                        elif all([x in info.attrib for x in
                                  ['components', 'material_phase']]):
                            for i in range(int(info.attrib['components'])):
                                entries.append('{e[material_phase]}%{e[name]}'
                                               '%{e[statistic]}%i'
                                               .format(e=info.attrib, i=i))
                        elif 'components' in info.attrib:
                            for i in range(int(info.get('components'))):
                                entries.append('{e[name]}%{e[statistic]}%i'
                                               .format(e=info.attrib, i=i))
                        elif 'material_phase' in info.attrib:
                            entries.append('{e[material_phase]}%{e[name]}'
                                           '%{e[statistic]}'
                                           .format(e=info.attrib))
                        else:
                            entries.append('{e[name]}%{e[statistic]}'
                                           .format(e=info.attrib))
                    elif line.startswith('</header>'):
                        break
            return numpy.asarray(entries), num

        if self.statfile[0].endswith('stat'):
            entries, num = GatherEntries(self.statfile[0])
            values = numpy.genfromtxt(self.statfile[0], skip_header=num + 1)
            values = values[:, numpy.argsort(entries)]
            entries = entries[numpy.argsort(entries)]
            if len(self.statfile) > 1:
                for item in self.statfile[1:]:
                    entriesTemp, num = GatherEntries(item)
                    if numpy.array_equal(
                            entries, entriesTemp[numpy.argsort(entriesTemp)]):
                        valuesTemp = numpy.genfromtxt(item,
                                                      skip_header=num + 1)
                        valuesTemp = valuesTemp.T[numpy.argsort(entriesTemp)].T
                        values = numpy.vstack((values, valuesTemp))
                    else:
                        sys.exit('Statfiles entries do not match')
        elif self.statfile[0].endswith('detectors'):
            entries = GatherEntries(self.statfile[0])[0]
            values = numpy.fromfile(self.statfile[0] + '.dat',
                                    dtype=numpy.float64)
            ncols = entries.size
            nrows = values.size // ncols
            values = values[:nrows * ncols].reshape(nrows, ncols)
        return list(entries), values

    def RefreshData(self):
        with (self.xCombo.handler_block(self.xCon),
              self.yCombo.handler_block(self.yCon)):
            self.xCombo.remove_all()
            self.yCombo.remove_all()
            self.entries, self.values = self.ReadData()
            self.InitCombo('reload', prevIterX=self.xCombo.get_active_iter(),
                           prevIterY=self.yCombo.get_active_iter())

    def ReleaseFocus(self, widget, event):
        if event.keyval == Gdk.KEY_Escape:
            self.canvas.grab_focus()


window = StatplotWindow(args.statfile)
window.connect('delete-event', Gtk.main_quit)
window.connect('destroy', Gtk.main_quit)
window.show_all()
Gtk.main()
