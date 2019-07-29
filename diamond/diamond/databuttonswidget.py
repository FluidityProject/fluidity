#!/usr/bin/env python

#    This file is part of Diamond.
#
#    Diamond is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Diamond is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Diamond.  If not, see <http://www.gnu.org/licenses/>.

from gi.repository import GObject as gobject
from gi.repository import Gtk as gtk

class DataButtonsWidget(gtk.HBox):

  __gsignals__ = { "revert" : (gobject.SignalFlags.RUN_LAST, gobject.TYPE_NONE, ()),
                   "store"  : (gobject.SignalFlags.RUN_LAST, gobject.TYPE_NONE, ())}

  def __init__(self):
    gtk.HBox.__init__(self)
    revertButton = gtk.Button()
    revertButton.set_label("Revert data")
    revertButton.connect("clicked", self._revert)

    storeButton = gtk.Button()
    storeButton.set_label("Store data")
    storeButton.connect("clicked", self._store)

    self.pack_start(revertButton, True, True, 0)
    self.pack_end(storeButton, True, True, 0)

    return

  def _revert(self, widget = None):
    self.emit("revert")

  def _store(self, widget = None):
    self.emit("store")

gobject.type_register(DataButtonsWidget)
