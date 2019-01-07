#!/usr/bin/python
#
# James Maddison
# Thomas Duvernay

import argparse
import numpy
import sys
import time
import warnings

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from matplotlib.backends.backend_gtk3cairo import FigureCanvasGTK3Cairo
from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

import fluidity.diagnostics.fluiditytools as fluidity_tools

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='GUI for Fluidity .stat outputs\n\nWhen the window is '
    + 'displayed, press the following keys for additional options:\n'
    + '- x, y -> Alternate the scale of the X, Y axis between linear and '
    + 'logarithmic.\n- l -> Alternate the representation of the data between '
    + 'a solid line and markers.\n- r -> Update the .stat file. Only relevant '
    + 'for simulations that are still running.\n- q -> Exit the GUI.')
parser.add_argument('statfile', nargs=1, help='path to the .stat file')
args = parser.parse_args(sys.argv[1:])


class StatplotWindow(Gtk.Window):
    def __init__(self, statfile):
        Gtk.Window.__init__(self, title=statfile[-1])
        self._file = statfile
        self.connect("key-press-event", self._KeyPressed)
        self._ReadData(statfile)
        self._vBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(self._vBox)
        self._hBox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL)
        self._vBox.pack_end(self._hBox, False, False, 0)
        self._xCombo = Gtk.ComboBoxText()
        self._yCombo = Gtk.ComboBoxText()
        paths = sorted(self._stat.Paths())
        for path in paths:
            self._xCombo.append_text(path)
            self._yCombo.append_text(path)
        if "ElapsedTime" in paths:
            ind = paths.index("ElapsedTime")
            iterX = self._xCombo.get_model().get_iter("%d" % ind)
            if ind == 0:
                iterY = self._yCombo.get_model().get_iter(1)
            else:
                iterY = self._yCombo.get_model().get_iter_first()
        else:
            iterX = self._xCombo.get_model().get_iter_first()
            iterY = self._yCombo.get_model().get_iter(1)
        self._xCombo.set_active_iter(iterX)
        self._yCombo.set_active_iter(iterY)
        self._hBox.pack_start(self._xCombo, True, True, 0)
        self._hBox.pack_end(self._yCombo, True, True, 0)
        self._fig, self._ax = plt.subplots(nrows=1, ncols=1, num=0)
        self._canvas = FigureCanvasGTK3Cairo(self._fig)
        self._toolbar = NavigationToolbar2GTK3(self._canvas, None)
        self._PlotType = 'line'
        self._PlotData(self._PlotType, 'linear', 'linear')
        self._fBox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self._fBox.pack_start(self._canvas, True, True, 0)
        self._fBox.pack_end(self._toolbar, False, False, 0)
        self._vBox.pack_start(self._fBox, True, True, 0)
        self._xCombo.connect("changed", self._ComboChangedX)
        self._yCombo.connect("changed", self._ComboChangedY)
        return

    def _ReadData(self, statfile):
        stats = []
        for i, filename in enumerate(statfile):
            failcount = 0
            while failcount < 5:
                try:
                    stats.append(fluidity_tools.Stat(filename))
                    break
                except (TypeError, ValueError):
                    time.sleep(0.2)
                    failcount += 1
            if failcount == 5:
                raise Exception("Could not open %s" % filename)
        if len(stats) == 1:
            self._stat = stats[0]
        else:
            self._stat = fluidity_tools.JoinStat(*stats)

    def _RefreshData(self, statfile, xData, yData):
        stats = []
        for i, filename in enumerate(statfile):
            failcount = 0
            while failcount < 5:
                try:
                    stats.append(fluidity_tools.Stat(filename))
                    break
                except (TypeError, ValueError):
                    time.sleep(0.2)
                    failcount += 1
            if failcount == 5:
                raise Exception("Could not open %s" % filename)
        if len(stats) == 1:
            self._stat = stats[0]
        else:
            self._stat = fluidity_tools.JoinStat(*stats)
        paths = sorted(self._stat.Paths())
        for path in paths:
            self._xCombo.append_text(path)
            self._yCombo.append_text(path)
        if xData in paths and yData in paths:
            indX = paths.index(xData)
            indY = paths.index(yData)
            iterX = self._xCombo.get_model().get_iter("%d" % indX)
            iterY = self._yCombo.get_model().get_iter("%d" % indY)
        else:
            iterX = self._xCombo.get_model().get_iter_first()
            iterY = self._yCombo.get_model().get_iter(1)
        self._xCombo.set_active_iter(iterX)
        self._yCombo.set_active_iter(iterY)

    def _ComboChangedX(self, widget):
        self._PlotData(self._PlotType, 'linear', self._ax.get_yscale())
        return

    def _ComboChangedY(self, widget):
        self._PlotData(self._PlotType, self._ax.get_xscale(), 'linear')
        return

    def _PlotData(self, type, xscale, yscale):
        self._xField = self._xCombo.get_active_text()
        self._yField = self._yCombo.get_active_text()
        self._xData = self._stat[self._xField]
        self._yData = self._stat[self._yField]
        self._ax.cla()
        if numpy.unique(self._xData).size == 1 \
                and numpy.unique(self._yData).size == 1:
            type = 'marker'
        if type == 'line':
            self._ax.plot(self._xData, self._yData, linewidth=2,
                          color='xkcd:ocean blue')
        elif type == 'marker':
            self._ax.plot(self._xData, self._yData, linestyle='none',
                          marker='D', markersize=7,
                          markerfacecolor='xkcd:ocean blue',
                          markeredgecolor='xkcd:black')
        majFor = tck.ScalarFormatter(useMathText=True)
        majFor.set_scientific(True)
        majFor.set_powerlimits((0, 0))
        if xscale == 'linear':
            self._ax.ticklabel_format(style='sci', axis='x',
                                      scilimits=(0, 0), useMathText=True)
            self._ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'Attempting')
                self._ax.set_xlim([min(self._xData) - 0.02 *
                                   (max(self._xData) - min(self._xData)),
                                   max(self._xData) + 0.02 *
                                   (max(self._xData) - min(self._xData))])
        elif xscale == 'log':
            self._ax.set_xscale('log')
            self._ax.xaxis.set_minor_locator(
                tck.LogLocator(subs=numpy.arange(2, 10)))
            logFmt = tck.LogFormatterSciNotation(base=10,
                                                 labelOnlyBase=False,
                                                 minor_thresholds=(4, 1))
            self._ax.xaxis.set_minor_formatter(logFmt)
            if min(self._xData) > 0:
                self._ax.set_xlim([10 ** (numpy.log10(min(self._xData)) -
                                          0.02 *
                                          numpy.log10(max(self._xData) /
                                                      min(self._xData))),
                                   10 ** (numpy.log10(max(self._xData)) +
                                          0.02 *
                                          numpy.log10(max(self._xData) /
                                                      min(self._xData)))])
            elif min(self._xData) == 0:
                self._ax.set_xlim([0, 10 ** (numpy.log10(max(self._xData)) +
                                             0.02 *
                                             numpy.log10(max(self._xData)))])
        elif xscale == 'symlog':
            xMin = min(abs(self._xData[self._xData != 0]))
            xMax = max(abs(self._xData))
            xRange = numpy.log10(xMax / xMin)
            self._ax.set_xscale('symlog', basex=10, subsx=numpy.arange(2, 10),
                                linthreshx=xMin * 10 ** (xRange / 2))
            symLogLoc = tck.SymmetricalLogLocator(
                subs=numpy.arange(2, 10), linthresh=xMin * 10 ** (xRange / 2),
                base=10)
            self._ax.xaxis.set_minor_locator(symLogLoc)
            logFmt = tck.LogFormatterSciNotation(
                base=10, labelOnlyBase=False, minor_thresholds=(4, 1),
                linthresh=xMin * 10 ** (xRange / 2))
            self._ax.xaxis.set_minor_formatter(logFmt)
        if yscale == 'linear':
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', 'Attempting')
                self._ax.set_ylim([min(self._yData) - 0.02 *
                                  (max(self._yData) - min(self._yData)),
                                  max(self._yData) + 0.02 *
                                  (max(self._yData) - min(self._yData))])
            self._ax.yaxis.set_major_formatter(majFor)
            self._ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
        elif yscale == 'log':
            self._ax.set_yscale('log')
            self._ax.yaxis. \
                set_minor_locator(tck.LogLocator(subs=numpy.arange(2, 10)))
            logFmt = tck.LogFormatterSciNotation(base=10,
                                                 labelOnlyBase=False,
                                                 minor_thresholds=(4, 1))
            self._ax.yaxis.set_minor_formatter(logFmt)
            if min(self._yData) > 0:
                self._ax.set_ylim([10 ** (numpy.log10(min(self._yData)) -
                                          0.02 *
                                          numpy.log10(max(self._yData) /
                                                      min(self._yData))),
                                   10 ** (numpy.log10(max(self._yData)) +
                                          0.02 *
                                          numpy.log10(max(self._yData) /
                                                      min(self._yData)))])
            elif min(self._yData) == 0:
                self._ax.set_ylim([0, 10 ** (numpy.log10(max(self._yData)) +
                                             0.02 *
                                             numpy.log10(max(self._yData)))])
        elif yscale == 'symlog':
            yMin = min(abs(self._yData[self._yData != 0]))
            yMax = max(abs(self._yData))
            yRange = numpy.log10(yMax / yMin)
            self._ax.set_yscale('symlog', basey=10, subsy=numpy.arange(2, 10),
                                linthreshy=yMin * 10 ** (yRange / 2))
            symLogLoc = tck.SymmetricalLogLocator(
                subs=numpy.arange(2, 10), linthresh=yMin * 10 ** (yRange / 2),
                base=10)
            self._ax.yaxis.set_minor_locator(symLogLoc)
            logFmt = tck.LogFormatterSciNotation(
                base=10, labelOnlyBase=False, minor_thresholds=(2, 1),
                linthresh=yMin * 10 ** (yRange / 2))
            self._ax.yaxis.set_minor_formatter(logFmt)
        self._ax.set_xlabel(self._xField, fontweight='bold', fontsize=20)
        self._ax.set_ylabel(self._yField, fontweight='bold', fontsize=20)
        self._ax.tick_params(which='major', length=7, labelsize=16, width=2)
        self._ax.tick_params(which='minor', length=4, labelsize=10, width=2,
                             colors='xkcd:scarlet', labelrotation=45)
        self._ax.xaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                             color='xkcd:black')
        self._ax.yaxis.get_offset_text().set(fontsize=13, fontweight='bold',
                                             color='xkcd:black')
        self._fig.set_tight_layout(True)
        self._fig.canvas.draw()
        return

    def _KeyPressed(self, widget, event):
        key = event.string
        if key == 'r':
            self._RefreshData(self._file, self._xCombo.get_active_text(),
                              self._yCombo.get_active_text())
            self._PlotData(self._PlotType, self._ax.get_xscale(),
                           self._ax.get_yscale())
        elif key == 'q':
            self.destroy()
        elif key == "x":
            if self._ax.get_xscale() == 'linear' \
                    and self._xData[self._xData != 0].size == 0:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the Y axis is a '
                              + 'constant.', stacklevel=2)
            elif self._ax.get_xscale() == 'linear' \
                    and max(abs(self._xData)) / min(abs(self._xData)) < 2:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the X axis has a range '
                              + 'of variation smaller than a factor two.',
                              stacklevel=2)
            elif self._ax.get_xscale() == 'linear':
                xMin = min(abs(self._xData[self._xData != 0]))
                xMax = max(abs(self._xData))
                xRange = numpy.log10(xMax / xMin)
                if xMin != xMax and min(self._xData) < 0:
                    self._ax.set_xscale('symlog', basex=10,
                                        linthreshx=xMin * 10 ** (xRange / 2),
                                        subsx=numpy.arange(2, 10))
                    symLogLoc = tck. \
                        SymmetricalLogLocator(subs=numpy.arange(2, 10),
                                              linthresh=xMin
                                              * 10 ** (xRange / 2), base=10)
                    self._ax.xaxis.set_minor_locator(symLogLoc)
                    logFmt = tck.LogFormatterSciNotation(
                        base=10, labelOnlyBase=False, minor_thresholds=(4, 1),
                        linthresh=xMin * 10 ** (xRange / 2))
                    self._ax.xaxis.set_minor_formatter(logFmt)
                    # Thomas Duvernay, 06/01/19
                    # There seems to be a bug with the labelling of the 0 tick
                    # when a 'symlog' is used as an axis scale. It looks like
                    # it is considered as a minor tick.
                elif xMin != xMax:
                    self._ax.set_xscale('log')
                    self._ax.xaxis.set_minor_locator(
                            tck.LogLocator(subs=numpy.arange(2, 10)))
                    logFmt = tck.LogFormatterSciNotation(
                        base=10, labelOnlyBase=False, minor_thresholds=(4, 1))
                    self._ax.xaxis.set_minor_formatter(logFmt)
                    if min(self._xData) > 0:
                        self._ax.set_xlim(
                            [10 ** (numpy.log10(min(self._xData)) - 0.02
                                    * numpy.log10(max(self._xData)
                                                  / min(self._xData))),
                             10 ** (numpy.log10(max(self._xData)) + 0.02
                                    * numpy.log10(max(self._xData)
                                                  / min(self._xData)))])
                    else:
                        self._ax.set_xlim(
                            [0, 10 ** (numpy.log10(max(self._xData)) + 0.02
                                       * numpy.log10(max(self._xData)))])
                else:
                    warnings.warn('Change to logarithmic scale denied: the '
                                  + 'selected variable for the X axis is a '
                                  + 'constant.', stacklevel=2)
            elif self._ax.get_xscale() in ['log', 'symlog']:
                self._ax.set_xscale('linear')
                self._ax.set_xlim([min(self._xData) - 0.02 *
                                   (max(self._xData) - min(self._xData)),
                                   max(self._xData) + 0.02 *
                                   (max(self._xData) - min(self._xData))])
                majFor = tck.ScalarFormatter(useMathText=True)
                majFor.set_scientific(True)
                majFor.set_powerlimits((0, 0))
                self._ax.xaxis.set_major_formatter(majFor)
                self._ax.xaxis.set_minor_locator(tck.AutoMinorLocator())
            self._fig.canvas.draw()
        elif key == "y":
            if self._ax.get_yscale() == 'linear' \
                    and self._yData[self._yData != 0].size == 0:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the Y axis is a '
                              + 'constant.', stacklevel=2)
            elif self._ax.get_yscale() == 'linear' \
                    and max(abs(self._yData)) / min(abs(self._yData)) < 2:
                warnings.warn('Change to logarithmic scale denied: the '
                              + 'selected variable for the Y axis has a range '
                              + 'of variation smaller than a factor two.',
                              stacklevel=2)
            elif self._ax.get_yscale() == 'linear':
                yMin = min(abs(self._yData[self._yData != 0]))
                yMax = max(abs(self._yData))
                yRange = numpy.log10(yMax / yMin)
                if yMin != yMax and min(self._yData) < 0:
                    self._ax.set_yscale('symlog', basey=10,
                                        linthreshy=yMin * 10 ** (yRange / 2),
                                        subsy=numpy.arange(2, 10))
                    symLogLoc = tck.SymmetricalLogLocator(
                        subs=numpy.arange(2, 10),
                        linthresh=yMin * 10 ** (yRange / 2), base=10)
                    self._ax.yaxis.set_minor_locator(symLogLoc)
                    logFmt = tck.LogFormatterSciNotation(
                        base=10, labelOnlyBase=False, minor_thresholds=(2, 1),
                        linthresh=yMin * 10 ** (yRange / 2))
                    self._ax.yaxis.set_minor_formatter(logFmt)
                    # Thomas Duvernay, 06/01/19
                    # There seems to be a bug with the labelling of the 0 tick
                    # when a 'symlog' is used as an axis scale. It looks like
                    # it is considered as a minor tick.
                elif yMin != yMax > 1:
                    self._ax.set_yscale('log')
                    self._ax.yaxis.set_minor_locator(
                            tck.LogLocator(subs=numpy.arange(2, 10)))
                    logFmt = tck.LogFormatterSciNotation(
                        base=10, labelOnlyBase=False, minor_thresholds=(4, 1))
                    self._ax.yaxis.set_minor_formatter(logFmt)
                    if min(self._yData) > 0:
                        self._ax.set_ylim(
                            [10 ** (numpy.log10(min(self._yData)) - 0.02
                                    * numpy.log10(max(self._yData)
                                                  / min(self._yData))),
                             10 ** (numpy.log10(max(self._yData)) + 0.02
                                    * numpy.log10(max(self._yData)
                                                  / min(self._yData)))])
                    else:
                        self._ax.set_ylim(
                            [0, 10 ** (numpy.log10(max(self._yData)) + 0.02
                                       * numpy.log10(max(self._yData)))])
                else:
                    warnings.warn('Change to logarithmic scale denied: the '
                                  + 'selected variable for the Y axis is a '
                                  + 'constant.', stacklevel=2)
            elif self._ax.get_yscale() in ['log', 'symlog']:
                self._ax.set_yscale('linear')
                self._ax.set_ylim([min(self._yData) - 0.02 *
                                   (max(self._yData) - min(self._yData)),
                                   max(self._yData) + 0.02 *
                                   (max(self._yData) - min(self._yData))])
                majFor = tck.ScalarFormatter(useMathText=True)
                majFor.set_scientific(True)
                majFor.set_powerlimits((0, 0))
                self._ax.yaxis.set_major_formatter(majFor)
                self._ax.yaxis.set_minor_locator(tck.AutoMinorLocator())
            self._fig.canvas.draw()
        elif key == 'l':
            self._PlotType = 'line' if self._PlotType == 'marker' else 'marker'
            self._PlotData(self._PlotType, self._ax.get_xscale(),
                           self._ax.get_yscale())
        return


window = StatplotWindow(args.statfile)
window.connect("delete-event", Gtk.main_quit)
window.set_default_size(1600, 900)
window.set_position(Gtk.WindowPosition.CENTER)
window.show_all()
Gtk.main()
