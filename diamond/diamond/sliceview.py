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

from . import attributewidget
from . import databuttonswidget
from . import datawidget
from . import mixedtree

class SliceView(gtk.Window):
  
  __gsignals__ = { "on-store" : (gobject.SignalFlags.RUN_LAST, gobject.TYPE_NONE, ()),
                   "update-name"  : (gobject.SignalFlags.RUN_LAST, gobject.TYPE_NONE, ())}

  def __init__(self, parent):
    gtk.Window.__init__(self)
    
    self.set_default_size(800, 600)
    self.set_title("Slice View")
    self.set_modal(True)
    self.set_transient_for(parent)

    mainvbox = gtk.VBox()
    self.vbox = gtk.VBox()

    scrolledWindow = gtk.ScrolledWindow()
    scrolledWindow.set_policy(gtk.PolicyType.NEVER, gtk.PolicyType.AUTOMATIC)
    scrolledWindow.add_with_viewport(self.vbox)

    self.databuttons = databuttonswidget.DataButtonsWidget()

    self.statusbar = gtk.Statusbar()

    mainvbox.pack_start(scrolledWindow, True, True, 0)
    mainvbox.pack_start(self.databuttons, False, True, 0)
    mainvbox.pack_start(self.statusbar, False, True, 0)

    self.add(mainvbox)
    self.show_all()

  def update(self, node, tree):
    nodes = self.get_nodes(node, tree)
    if not nodes:
      self.destroy()

    for n in nodes:
      self.vbox.pack_start(self.control(n), True, True, 0)

    maxwidth = 0
    for child in self.vbox.get_children():
      width, height = child.label.get_size_request()
      maxwidth = max(maxwidth, width)
    
    for child in self.vbox.get_children():
      child.label.set_size_request(maxwidth, -1)

    self.check_resize()

  def get_nodes(self, node, tree):
   nodes = []

   for child in tree.get_children():
     if child.active:
       if child.name == node.name and child.is_sliceable():
         nodes.append(child.get_mixed_data())
       nodes += self.get_nodes(node, child)
   
   return nodes

  def control(self, node):
    hbox = gtk.HBox()

    label = gtk.Label(node.get_name_path())
    hbox.label = label

    data = datawidget.DataWidget()
    data.geometry_dim_tree = self.geometry_dim_tree
    data.connect("on-store", self.on_store)
    data.set_buttons(self.databuttons)
    data.update(node)

    attributes = attributewidget.AttributeWidget()
    attributes.connect("on-store", self.on_store)
    attributes.connect("update-name", self.update_name)
    attributes.update(node)

    hbox.pack_start(label, True, True, 0)
    hbox.pack_start(data, True, True, 0)
    hbox.pack_start(attributes, True, True, 0)
    
    hbox.show_all()

    return hbox

  def on_store(self, widget = None):
    self.emit("on-store")

  def update_name(self, widget = None):
    self.emit("update-name") 

gobject.type_register(SliceView)
