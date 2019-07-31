#!/usr/bin/env python3

import gi
gi.require_version('Gdk', '3.0')
gi.require_version('Gtk', '3.0')
from gi.repository import Gdk, Gtk  # noqa: E402

import argparse  # noqa: E402
import fluidity.diagnostics.fluiditytools as fluidity_tools  # noqa: E402
from matplotlib import rc  # noqa: E402
from matplotlib.backends.backend_gtk3agg import (
    FigureCanvasGTK3Agg as FigureCanvas)  # noqa: E402
from matplotlib.backends.backend_gtk3 import (
    NavigationToolbar2GTK3 as NavigationToolbar)  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
import matplotlib.ticker as tck  # noqa: E402
import numpy  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
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
- Esc  -> Release focus on text boxes.
- q    -> Exit the GUI.''')
parser.add_argument('statfile', nargs='+', help='path to the .stat file')
args = parser.parse_args(sys.argv[1:])


class StatplotWindow(Gtk.Window):
    def __init__(self, statfile):
        Gtk.Window.__init__(self)
        self.connect('key-press-event', self.KeyPressed)
        self.set_border_width(8)
        self.set_default_size(1600, 900)
        self.set_position(Gtk.WindowPosition.CENTER)
        self.set_title(statfile[-1])
        self.statfile = statfile
        self.vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=8)
        self.vbox.set_homogeneous(False)
        self.add(self.vbox)
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=8)
        hbox.set_homogeneous(True)
        self.vbox.pack_end(hbox, False, False, 0)
        self.entries = self.ReadData(statfile)
        self.xCombo = Gtk.ComboBoxText.new_with_entry()
        self.xCombo.set_wrap_width(3)
        self.yCombo = Gtk.ComboBoxText.new_with_entry()
        self.yCombo.set_wrap_width(3)
        self.PopulateCombo('load')
        self.InitCompletion(self.xCombo)
        self.InitCompletion(self.yCombo)
        self.xCon = self.xCombo.connect('changed', self.ComboChangedX)
        self.xCombo.connect('key-release-event', self.ReleaseFocus)
        self.yCon = self.yCombo.connect('changed', self.ComboChangedY)
        self.yCombo.connect('key-release-event', self.ReleaseFocus)
        hbox.pack_start(self.xCombo, True, True, 0)
        hbox.pack_start(self.yCombo, True, True, 0)
        self.fig = Figure(figsize=(14, 7))
        self.ax = self.fig.gca()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.set_can_focus(True)
        self.vbox.pack_start(self.canvas, True, True, 0)
        self.PlotType = 'line'
        self.PlotData('create', self.PlotType, 'linear', 'linear')
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.vbox.pack_start(self.toolbar, False, False, 0)

    def ComboChangedX(self, widget):
        if self.xCombo.get_active_text() in sorted(self.entries.Paths()):
            self.PlotData('update', self.PlotType, 'linear',
                          self.ax.get_yscale())
            self.canvas.grab_focus()

    def ComboChangedY(self, widget):
        if self.yCombo.get_active_text() in sorted(self.entries.Paths()):
            self.PlotData('update', self.PlotType, self.ax.get_xscale(),
                          'linear')
            self.canvas.grab_focus()

    def FormatAxis(self, ax2fmt, scale):
        if ax2fmt == 'x':
            self.set_scale = self.ax.set_xscale
            self.axis = self.ax.xaxis
            self.Data = self.xData
        elif ax2fmt == 'y':
            self.set_scale = self.ax.set_yscale
            self.axis = self.ax.yaxis
            self.Data = self.yData
        self.set_scale(scale)
        if scale == 'linear':
            self.ax.ticklabel_format(style='sci', axis=ax2fmt,
                                     scilimits=(0, 0), useMathText=True)
            self.axis.set_minor_locator(tck.AutoMinorLocator())
            self.ax.relim()
            self.ax.autoscale(True, ax2fmt, None)
        elif scale == 'log':
            self.axis.set_minor_locator(
                tck.LogLocator(subs=numpy.arange(2, 10)))
            logFmt = tck.LogFormatterSciNotation(base=10,
                                                 labelOnlyBase=False,
                                                 minor_thresholds=(4, 1))
            self.axis.set_minor_formatter(logFmt)
            self.ax.relim()
            self.ax.autoscale(True, ax2fmt, None)
        elif scale == 'symlog':
            axMin = min(abs(self.Data[self.Data != 0]))
            axMax = max(abs(self.Data))
            axRange = numpy.log10(axMax / axMin)
            if ax2fmt == 'x':
                self.set_scale(
                    'symlog', basex=10, subsx=numpy.arange(2, 10),
                    linthreshx=axMin * 10 ** (axRange / 2))
            elif ax2fmt == 'y':
                self.set_scale(
                    'symlog', basey=10, subsy=numpy.arange(2, 10),
                    linthreshy=axMin * 10 ** (axRange / 2))
            # Thomas Duvernay, 06/01/19
            # There seems to be a bug with the labelling of the 0 tick
            # when a 'symlog' is used as an axis scale. It looks like
            # it is considered as a minor tick.
            symLogLoc = tck.SymmetricalLogLocator(
                subs=numpy.arange(2, 10),
                linthresh=axMin * 10 ** (axRange / 2), base=10)
            self.axis.set_minor_locator(symLogLoc)
            logFmt = tck.LogFormatterSciNotation(
                base=10, labelOnlyBase=False, minor_thresholds=(4, 1),
                linthresh=axMin * 10 ** (axRange / 2))
            self.axis.set_minor_formatter(logFmt)
        self.ax.set_xlabel(self.xField, fontweight='bold', fontsize=20)
        self.ax.set_ylabel(self.yField, fontweight='bold', fontsize=20)
        self.ax.tick_params(which='major', length=7, labelsize=16, width=2)
        self.ax.tick_params(which='minor', length=4, labelsize=10, width=2,
                            colors='xkcd:scarlet', labelrotation=45)
        self.ax.xaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                            color='xkcd:black')
        self.ax.yaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                            color='xkcd:black')
        self.fig.set_tight_layout(True)

    def InitCompletion(self, comboBox):
        completion = Gtk.EntryCompletion.new()
        completion.set_text_column(0)
        completion.set_inline_completion(True)
        completion.set_inline_selection(True)
        completion.set_model(comboBox.get_model())
        comboBox.get_child().set_completion(completion)

    def KeyPressed(self, widget, event):
        key = event.string
        if key == 'r':
            self.RefreshData(self.statfile)
            self.PlotData('update', self.PlotType, self.ax.get_xscale(),
                          self.ax.get_yscale())
        elif key == 'q':
            self.destroy()
        elif key == 'x' or key == 'y':
            if key == 'x':
                self.get_scale = self.ax.get_xscale
                self.Data = self.xData
            elif key == 'y':
                self.get_scale = self.ax.get_yscale
                self.Data = self.yData
            if self.get_scale() == 'linear' \
                    and self.Data[self.Data != 0].size == 0:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the ' + key + ' axis '
                              + 'is null.', stacklevel=2)
                scale = 'linear'
            elif self.get_scale() == 'linear' \
                    and max(abs(self.Data)) / min(abs(self.Data)) < 100:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the ' + key + ' axis '
                              + 'has a range of variation smaller than two '
                              + 'orders of magnitude.', stacklevel=2)
                scale = 'linear'
            elif self.get_scale() == 'linear':
                axMin = min(abs(self.Data[self.Data != 0]))
                axMax = max(abs(self.Data))
                if axMin != axMax and min(self.Data) < 0:
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

    def PlotData(self, action, type, xscale, yscale):
        self.xField = self.xCombo.get_active_text()
        self.yField = self.yCombo.get_active_text()
        self.xData = self.entries[self.xField]
        self.yData = self.entries[self.yField]
        if numpy.unique(self.xData).size == 1 \
                and numpy.unique(self.yData).size == 1:
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

    def PopulateCombo(self, action, prevIterX=None, prevIterY=None):
        for entry in sorted(self.entries.Paths()):
            self.xCombo.append_text(entry)
            self.yCombo.append_text(entry)
        if action == 'load':
            if 'ElapsedTime' in self.entries.Paths():
                iterX = self.xCombo.get_model(). \
                    get_iter(sorted(self.entries.Paths()).index('ElapsedTime'))
                self.xCombo.set_active_iter(iterX)
                self.yCombo.set_active(0)
            else:
                self.xCombo.set_active(0)
                self.yCombo.set_active(1)
        elif action == 'reload':
            self.xCombo.set_active_iter(prevIterX)
            self.yCombo.set_active_iter(prevIterY)

    def ReadData(self, statfile):
        stats = []
        for i, filename in enumerate(statfile):
            failcount = 0
            while failcount < 4:
                try:
                    stats.append(fluidity_tools.Stat(filename))
                    break
                except (TypeError, ValueError):
                    time.sleep(0.2)
                    failcount += 1
            if failcount == 4:
                stats.append(fluidity_tools.Stat(filename))
        if len(stats) == 1:
            return stats[0]
        else:
            return fluidity_tools.JoinStat(*stats)

    def RefreshData(self, statfile):
        with self.xCombo.handler_block(self.xCon), \
                self.yCombo.handler_block(self.yCon):
            self.xCombo.remove_all()
            self.yCombo.remove_all()
            self.entries = self.ReadData(statfile)
            self.PopulateCombo('reload',
                               prevIterX=self.xCombo.get_active_iter(),
                               prevIterY=self.yCombo.get_active_iter())

    def ReleaseFocus(self, widget, event):
        if event.keyval == Gdk.KEY_Escape:
            self.canvas.grab_focus()


window = StatplotWindow(args.statfile)
window.connect('delete-event', Gtk.main_quit)
window.connect('destroy', Gtk.main_quit)
window.show_all()
Gtk.main()
