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

import os
import os.path
import re
import sys
import tempfile
import webbrowser
import cStringIO as StringIO

import pango
import gobject
import gtk
import gtk.glade

import debug
import dialogs
import choice
import config
import plist
import schema
import scherror
import tree
import plugins
import StringIO
import TextBufferMarkup

from lxml import etree

try:
  gtk.Tooltip()
except:
  debug.deprint("Interface warning: Unable to use GTK tooltips")

"""
Here are some notes about the code:

Important fields:
  file_path: the directory containing the current file (the working directory or directory of last opened / saved file if no file is open)
  filename: output filename of the current file
  data_paths: paths (from the root node) to important Diamond data (e.g. geometry dimension)
  geometry_dim_tree: MixedTree, with parent equal to the geometry dimension tree and child equal to the geometry dimension data subtree
  gladefile: input Glade file
  gui: GUI GladeXML object
  logofile: the GUI logo file
  main_window: GUI toplevel window
  node_attrs: RHS attributes entry widget
  node_comment: RHS comment entry widget
  node_comment_interacted: used to determine if the comment widget has been interacted with without the comment being stored
  node_data: RHS data entry widget
  node_data_buttons_hbox: container for "Revert Data" and "Store Data" buttons
  node_data_interacted: used to determine if a node data widget has been interacted with without data being stored
  node_data_frame: frame containing data entry widgets
  node_desc: RHS description widget
  node_desc_link_bounds: a list of tuples corresponding to the start and end points of links in the current tree.Tree / MixedTree documentation
  options_tree_select_func_enabled: boolean, true if the options tree select function is enabled (used to overcome a nasty clash with the treeview clicked signal) - re-enabled on next options_tree_select_func call
  selected_node: a tree.Tree or MixedTree containing data to be displayed on the RHS
  selected_iter: last iter set by on_select_row
  s: current schema
  saved: boolean, false if the current file has been edited
  schemafile: the current RNG schema file
  schemafile_path: the directory containing the current schema file (the working directory if no schema is open)
  signals: dictionary containing signal handlers for signals set up in the Glade file
  statusbar: GUI status bar
  tree: LHS tree root
  treestore: the LHS tree model
  treeview: the LHS tree widget

Important routines:
  cellcombo_edited: called when a choice is selected on the left-hand pane
  init_options_frame: initialise the right-hand side
  init_treemodel: set up the treemodel and treeview
  on_treeview_clicked: when a row is clicked, process the consequences (e.g. activate inactive instance)
  set_treestore: stuff the treestore with a given tree.Tree
  on_find_find_button & search_treestore: the find functionality
  on_select_row: when a row is selected, update the options frame
  update_options_frame: paint the right-hand side

If there are bugs in reading in, see schema.read.
If there are bugs in writing out, see tree.write.
"""

class Diamond:
  def __init__(self, gladefile, schemafile = None, schematron_file = None, logofile = None, input_filename = None, 
      dim_path = "/geometry/dimension", suffix=None):
    self.gladefile = gladefile
    self.gui = gtk.glade.XML(self.gladefile)

    self.statusbar = DiamondStatusBar(self.gui.get_widget("statusBar"))
    self.find      = DiamondFindDialog(self, gladefile)
    self.popup = self.gui.get_widget("popupmenu")

    self.plugin_buttonbox = self.gui.get_widget("plugin_buttonbox")
    self.plugin_buttonbox.set_layout(gtk.BUTTONBOX_START)
    self.plugin_buttonbox.show()
    self.plugin_buttons = []

    self.scherror  = scherror.DiamondSchemaError(self, gladefile, schemafile, schematron_file)

    signals     =  {"on_new": self.on_new,
                    "on_quit": self.on_quit,
                    "on_open": self.on_open,
                    "on_open_schema": self.on_open_schema,
                    "on_save": self.on_save,
                    "on_save_as": self.on_save_as,
                    "on_validate": self.scherror.on_validate,
                    "on_validate_schematron": self.scherror.on_validate_schematron,
                    "on_expand_all": self.on_expand_all,
                    "on_collapse_all": self.on_collapse_all,
                    "on_find": self.find.on_find,
                    "on_go_to_node": self.on_go_to_node,
                    "on_console": self.on_console,
                    "on_display_properties_toggled": self.on_display_properties_toggled,
                    "on_about": self.on_about,
                    "on_copy_spud_path": self.on_copy_spud_path,
                    "on_copy": self.on_copy,
                    "on_paste": self.on_paste}
    self.gui.signal_autoconnect(signals)

    self.main_window = self.gui.get_widget("mainWindow")
    self.main_window.connect("delete_event", self.on_delete)

    self.logofile = logofile
    if self.logofile is not None:
      gtk.window_set_default_icon_from_file(self.logofile)

    self.init_treemodel()

    self.data_paths = {}
    self.data_paths["dim"] = dim_path

    self.suffix = suffix

    self.selected_node = None
    self.init_options_frame()

    self.file_path = os.getcwd()
    self.schemafile_path = os.getcwd()
    self.filename = None
    self.schemafile = None
    self.init_datatree()    
    self.set_saved(True)
    self.open_file(schemafile = None, filename = None)

    self.main_window.show()

    if not schemafile is None:
      self.open_file(schemafile = schemafile, filename = input_filename)

    # Hide specific menu items
    menu = self.gui.get_widget("menu")
    
    # Disable Find
    menu.get_children()[1].get_submenu().get_children()[0].set_property("visible", False)

    if schematron_file is None:
      # Disable Validate Schematron
      menu.get_children()[3].get_submenu().get_children()[1].set_property("sensitive", False)

    return

  def program_exists(self, name):
    ret = os.system("which %s > /dev/null" % name)
    return ret == 0

  ### MENU ###

  def update_title(self):
    """
    Update the Diamond title based on the save status of the currently open file.
    """

    title = "Diamond: "
    if not self.saved:
      title += "*"
    if self.filename is None:
      title += "(Unsaved)"
    else:
      title += os.path.basename(self.filename)
      if len(os.path.dirname(self.filename)) > 0:
        title += " (%s)" % os.path.dirname(self.filename)

    self.main_window.set_title(title)

    return

  def set_saved(self, saved, filename = ""):
    """
    Change the save status of the current file.
    """

    self.saved = saved
    if not filename == "":
      self.filename = filename
      if not filename is None:
        self.file_path = os.path.dirname(filename) + os.path.sep
    self.update_title()

    return

  def close_schema(self):
    if self.schemafile is None:
      return

    # clear the schema.
    self.s = None
    self.schemafile = None
    self.schemafile_path = None
    self.scherror.schema_file = None
       
    return

  def load_schema(self, schemafile):
    # so, if the schemafile has already been opened, then ..
    if schemafile == self.schemafile:
      self.statusbar.set_statusbar('Schema ' + schemafile + ' already loaded')
      return

    # if we aren't using a http schema, and we're passed a relative filename, we
    # need to absolut-ify it.
    if 'http' not in schemafile:
      schemafile = os.path.abspath(schemafile)

    self.statusbar.set_statusbar('Loading schema from ' + schemafile)

    # now, let's try and read the schema.
    try:
      s_read = schema.Schema(schemafile)
      self.s = s_read
      self.statusbar.set_statusbar('Loaded schema from ' + schemafile)
    except:
      dialogs.error_tb(self.main_window, "Unable to open schema file \"" + schemafile + "\"")
      self.statusbar.clear_statusbar()
      return

    self.schemafile = schemafile
    self.schemafile_path = os.path.dirname(schemafile) + os.path.sep
    self.scherror.schema_file = schemafile
    
    self.remove_children(None)
    self.init_datatree()
    
    return

  def close_file(self):
    if self.filename is None:
      return
  
    self.remove_children(None)
    self.init_datatree()

    self.filename = None
    
    return

  def load_file(self, filename):
    # if we have a relative path, make it absolute
    filename = os.path.abspath(filename)
  
    if filename == self.filename:
      return
  
    try:
      os.stat(filename)
    except OSError:
      self.filename = filename
      self.set_saved(False)
      
      self.remove_children(None)      
      self.init_datatree()
      
      return
      
    try:
      tree_read = self.s.read(filename)
      
      # Extract and display validation errors
      saved = True
      lost_eles, added_eles, lost_attrs, added_attrs = self.s.read_errors()
      if len(lost_eles) > 0 or len(added_eles) > 0 or len(lost_attrs) > 0 or len(added_attrs) > 0:
        saved = False
        msg = ""
        if len(lost_eles) > 0:
          msg += "Warning: lost xml elements:\n"
          for ele in lost_eles:
            msg += ele + "\n"
        if len(added_eles) > 0:
          msg += "Warning: added xml elements:\n"
          for ele in added_eles:
            msg += ele + "\n"
        if len(lost_attrs) > 0:
          msg += "Warning: lost xml attributes:\n"
          for ele in lost_attrs:
            msg += ele + "\n"
        if len(added_attrs) > 0:
          msg += "Warning: added xml attributes:\n"
          for ele in added_attrs:
            msg += ele + "\n"
      
        dialogs.long_message(self.main_window, msg)
        
      self.tree = tree_read
      self.filename = filename
    except:
      dialogs.error_tb(self.main_window, "Unable to open file \"" + filename + "\"")
      return

    self.set_saved(saved, filename)
      
    return

  def open_file(self, schemafile = "", filename = ""):
    """
    Handle opening or clearing of the current file and / or schema.
    """

    self.find.on_find_close_button()
    if schemafile is None:
      self.close_schema()
    elif not schemafile == "":
      self.load_schema(schemafile)
    if filename is None:
      self.close_file()
    elif not filename == "":
      self.load_file(filename)
      
    self.treeview.freeze_child_notify()
    self.treeview.set_model(None)
    self.set_treestore(None, [self.tree], True)
    self.treeview.set_model(self.treestore)
    self.treeview.thaw_child_notify()

    self.set_geometry_dim_tree()

    self.treeview.get_selection().unselect_all()

    self.selected_node = None
    self.update_options_frame()

    self.scherror.destroy_error_list()

    return

  def on_new(self, widget=None):
    """
    Called when new is clicked. Clear the treestore and reset the datatree.
    """

    if not self.saved:
      prompt_response = dialogs.prompt(self.main_window, "Unsaved data. Do you wish to continue?", gtk.MESSAGE_WARNING)
      if not prompt_response == gtk.RESPONSE_YES:
        return
    self.open_file(filename = None)
    self.filename = None

    return

  def on_open(self, widget=None):
    """
    Called when open is clicked. Open a user supplied file.
    """

    if not self.saved:
      prompt_response = dialogs.prompt(self.main_window, "Unsaved data. Do you want to save the current document before continuing?", gtk.MESSAGE_WARNING)
      if not prompt_response == gtk.RESPONSE_YES:
        return

    filter_names_and_patterns = {}
    if self.suffix is None:
      for xmlname in config.schemata:
        filter_names_and_patterns[config.schemata[xmlname][0]] = "*." + xmlname
    elif self.suffix in config.schemata.keys():
      filter_names_and_patterns[config.schemata[self.suffix][0]] = "*." + self.suffix
    else:
      filter_names_and_patterns[self.suffix] = "*." + self.suffix

    filename = dialogs.get_filename(title = "Open XML file", action = gtk.FILE_CHOOSER_ACTION_OPEN, filter_names_and_patterns = filter_names_and_patterns, folder_uri = self.file_path)

    if not filename is None:
      self.open_file(filename = filename)

    return

  def on_open_schema(self, widget=None):
    """
    Called when open schema is clicked. Clear the treestore and reset the schema.
    """

    if not self.saved:
      prompt_response = dialogs.prompt(self.main_window, "Unsaved data. Do you wish to continue?", gtk.MESSAGE_WARNING)
      if not prompt_response == gtk.RESPONSE_YES:
        return

    filename = dialogs.get_filename(title = "Open RELAX NG schema", action = gtk.FILE_CHOOSER_ACTION_OPEN, filter_names_and_patterns = {"RNG files":"*.rng"}, folder_uri = self.schemafile_path)
    if not filename is None:
      self.open_file(schemafile = filename)

    return

  def on_save(self, widget=None):
    """
    Write out to XML. If we don't already have a filename, open a dialog to get
    one.
    """

    self.node_data_store()

    if self.filename is None:
      return self.on_save_as(widget)
    else:
      self.statusbar.set_statusbar("Saving ...")
      self.main_window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.WATCH))
      try:
        self.tree.write(self.filename)
      except:
        dialogs.error_tb(self.main_window, "Saving to \"" + self.filename + "\" failed")
        self.statusbar.clear_statusbar()
        self.main_window.window.set_cursor(None)
        return False

      self.set_saved(True)

      self.statusbar.clear_statusbar()
      self.main_window.window.set_cursor(None)
      return True

    return False

  def on_save_as(self, widget=None):
    """
    Write out the XML to a file.
    """

    if self.schemafile is None:
      dialogs.error(self.main_window, "No schema file open")
      return False

    filter_names_and_patterns = {}
    if self.suffix is None:
      for xmlname in config.schemata:
        filter_names_and_patterns[config.schemata[xmlname][0]] = "*." + xmlname
    elif self.suffix in config.schemata.keys():
      filter_names_and_patterns[config.schemata[self.suffix][0]] = "*." + self.suffix
    else:
      filter_names_and_patterns[self.suffix] = "*." + self.suffix

    filename = dialogs.get_filename(title = "Save XML file", action = gtk.FILE_CHOOSER_ACTION_SAVE, filter_names_and_patterns = filter_names_and_patterns, folder_uri = self.file_path)

    if not filename is None:
      # Check that the selected file has a file extension. If not, add a .xml extension.
      if len(filename.split(".")) <= 1:
        filename += ".xml"

      # Save the file
      self.statusbar.set_statusbar("Saving ...")
      self.main_window.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.WATCH))
      self.tree.write(filename)
      self.set_saved(True, filename)
      self.statusbar.clear_statusbar()
      self.main_window.window.set_cursor(None)
      return True

    return False

  def on_delete(self, widget, event):
    """
    Called when the main window is deleted. Return "True" to prevent the deletion
    of the main window (deletion is handled by "on_quit").
    """

    self.on_quit(widget, event)

    return True

  def on_quit(self, widget, event = None):
    """
    Quit the program. Prompt the user to save data if the current file has been
    changed.
    """

    if self.saved:
      self.destroy()
    else:
      prompt_response = dialogs.prompt(self.main_window, "Unsaved data. Do you wish to save?", gtk.MESSAGE_WARNING, True)
      if prompt_response == gtk.RESPONSE_YES:
        if self.filename is None:
          retval = self.on_save_as()
        else:
          retval = self.on_save()
        if retval is True:
          self.destroy()
      elif prompt_response == gtk.RESPONSE_NO:
        self.destroy()

    return

  def destroy(self):
    """
    End the program.
    """

    try:
      gtk.main_quit()
    except:
      debug.dprint("Failed to quit - already quit?")

    return

  def on_display_properties_toggled(self, widget=None):
    self.options_frame.set_property("visible", not self.options_frame.get_property("visible"))
    return

  def set_cell_node(self, iter):
      node = self.get_painted_tree(iter)

      if isinstance(node, MixedTree):
        node = node.child

      if node.data is not None:
        self.data_renderer.set_property("foreground", "black")
        self.treestore.set_value(iter, 4, node.data)
      elif isinstance(node, tree.Tree) and not node.not_editable():
        self.data_renderer.set_property("foreground", "gray")

        datatype = ""

        if isinstance(node.datatype, plist.List):
          datatype = "(" + node.datatype.datatype.__name__ + ")"
        elif isinstance(node.datatype, tuple):
          datatype = str(node.datatype)
        else:
          datatype = node.datatype.__name__

        self.treestore.set_value(iter, 4, datatype)

  def expand_fill_data(self, model, path, iter, user_data):
      return self.set_cell_node(iter)

  def on_go_to_node(self, widget=None):
   """
   Go to a node, identified by an XPath
   """

   dialog = dialogs.GoToDialog(self)
   spudpath = dialog.run()

   return

  def on_expand_all(self, widget=None):
    """
    Show the whole tree.
    """

    self.treeview.expand_all()
    self.treestore.foreach(self.expand_fill_data, None)

    return

  def on_collapse_all(self, widget=None):
    """
    Collapse the whole tree.
    """

    self.treeview.collapse_all()

    return
    
  def on_console(self, widget = None):
    """
    Launch a python console
    """    
    
    # Construct the dictionary of locals that will be used by the interpreter
    locals = {}
    locals["interface"] = globals()
    locals["diamond_gui"] = self
  
    dialogs.console(self.main_window, locals)
    
    return

  def on_about(self, widget=None):
    """
    Tell the user how fecking great we are.
    """

    about = gtk.AboutDialog()
    about.set_name("Diamond")
    about.set_copyright("GPLv3")
    about.set_comments("A RELAX-NG-aware XML editor")
    about.set_authors(["Patrick E. Farrell", "James R. Maddison", "Matthew T. Whitworth", "Fraser J. Waters"])
    about.set_license("Diamond is free software: you can redistribute it and/or modify\n"+
                      "it under the terms of the GNU General Public License as published by\n"+
                      "the Free Software Foundation, either version 3 of the License, or\n"+
                      "(at your option) any later version.\n"+
                      "\n"+
                      "Diamond is distributed in the hope that it will be useful,\n"+
                      "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"+
                      "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"+
                      "GNU General Public License for more details.\n"+
                      "You should have received a copy of the GNU General Public License\n"+
                      "along with Diamond.  If not, see http://www.gnu.org/licenses/.")

    logo = gtk.gdk.pixbuf_new_from_file(self.logofile)

    try:
      image = about.get_children()[0].get_children()[0].get_children()[0]
      image.set_tooltip_text("Diamond: it's clearer than GEM")
    except:
      pass

    about.set_logo(logo)
    about.connect("destroy", dialogs.close_dialog)
    about.connect("response", dialogs.close_dialog)
    about.show()

    return


  def on_copy_spud_path(self, widget=None):
    path = self.get_selected_row(self.treeview.get_selection())
    if path is None:
      debug.deprint("No selection.")
      return
    iter = self.treestore.get_iter(path)
    active_tree = self.treestore.get_value(iter, 3)
    name = self.get_spudpath(active_tree)
    clipboard = gtk.clipboard_get()
    clipboard.set_text(name)
    clipboard.store()

  def on_copy(self, widget=None):
    if isinstance(self.selected_node, MixedTree):
      node = self.selected_node.parent
    else:
      node = self.selected_node    

    if node != None and node.active:
      ios = StringIO.StringIO()
      node.write(ios)
    
      clipboard = gtk.clipboard_get()
      clipboard.set_text(ios.getvalue())
      clipboard.store()

      ios.close()
    return

  def on_paste(self, widget=None):
    clipboard = gtk.clipboard_get()
    ios = StringIO.StringIO(clipboard.wait_for_text())
    
    if self.selected_iter is not None:    
      node = self.treestore.get_value(self.selected_iter, 3)

    if node != None:

      newnode = self.s.read(ios, node)

      if newnode is None:
        self.statusbar.set_statusbar("Trying to paste invalid XML.")
        return

      if not node.active:
        self.expand_tree(self.selected_iter)

      # Extract and display validation errors
      lost_eles, added_eles, lost_attrs, added_attrs = self.s.read_errors()
      if len(lost_eles) > 0 or len(added_eles) > 0 or len(lost_attrs) > 0 or len(added_attrs) > 0:
        saved = False
        msg = ""
        if len(lost_eles) > 0:
          msg += "Warning: lost xml elements:\n"
          for ele in lost_eles:
            msg += ele + "\n"
        if len(added_eles) > 0:
          msg += "Warning: added xml elements:\n"
          for ele in added_eles:
            msg += ele + "\n"
        if len(lost_attrs) > 0:
          msg += "Warning: lost xml attributes:\n"
          for ele in lost_attrs:
            msg += ele + "\n"
        if len(added_attrs) > 0:
          msg += "Warning: added xml attributes:\n"
          for ele in added_attrs:
            msg += ele + "\n"
      
        dialogs.long_message(self.main_window, msg)
 
      self.set_saved(False)     
   
      self.treeview.freeze_child_notify()
      iter = self.set_treestore(self.selected_iter, [newnode], True, True)
      self.treeview.thaw_child_notify()
      
      self.treeview.get_selection().select_iter(iter)

    return


  ## LHS ###

  def init_datatree(self):
    """
    Add the root node of the XML tree to the treestore, and its children.
    """

    if self.schemafile is None:
      self.set_treestore(None, [])
      self.tree = None
    else:
      l = self.s.valid_children(":start")

      self.tree = l[0]
      self.set_treestore(None, l)

    root_iter = self.treestore.get_iter_first()
    self.treeview.freeze_child_notify()
    self.treeview.set_model(None)
    self.expand_treestore(root_iter)
    self.treeview.set_model(self.treestore)
    self.treeview.thaw_child_notify()

    return

  def init_treemodel(self):
    """
    Set up the treestore and treeview.
    """

    self.treeview = optionsTree = self.gui.get_widget("optionsTree")
    self.treeview.connect("row-collapsed", self.on_treeview_row_collapsed)
    self.treeview.connect("key_press_event", self.on_treeview_key_press)
    self.treeview.connect("button_press_event", self.on_treeview_button_press)
    self.treeview.connect("popup_menu", self.on_treeview_popup)
    try:  # allow for possibility of no tooltips (like elsewhere)
      self.treeview.connect("query-tooltip", self.on_tooltip)
      self.treeview.set_property("has-tooltip", False)
    except:
      pass

    self.treeview.set_property("rules-hint", True)

    model = gtk.ListStore(str, str, gobject.TYPE_PYOBJECT)
    self.cellcombo = cellCombo = gtk.CellRendererCombo()
    cellCombo.set_property("model", model)
    cellCombo.set_property("text-column", 0)
    cellCombo.set_property("editable", True)
    cellCombo.set_property("has-entry", False)

    # Node column
    column = gtk.TreeViewColumn("Node", cellCombo, text=0)
    column.set_property("expand", True)
    column.set_resizable(True)
    column.set_cell_data_func(cellCombo, self.set_combobox_liststore)

    self.choicecell = choiceCell = gtk.CellRendererPixbuf()
    column.pack_end(choiceCell, expand=False)
    column.set_cell_data_func(choiceCell, self.set_cellpicture_choice)
    optionsTree.append_column(column)

    self.imgcell = cellPicture = gtk.CellRendererPixbuf()
    self.imgcolumn = imgcolumn = gtk.TreeViewColumn("", cellPicture)
    imgcolumn.set_property("expand", False)
    imgcolumn.set_property("fixed-width", 20)
    imgcolumn.set_property("sizing", gtk.TREE_VIEW_COLUMN_FIXED)
    imgcolumn.set_cell_data_func(cellPicture, self.set_cellpicture_cardinality)
    optionsTree.append_column(imgcolumn)

    # mtw07 - add column for quick preview of data.
    self.data_renderer = gtk.CellRendererText()
    self.data_renderer.set_property("editable", True)
    self.data_renderer.connect("edited", self.on_cell_edit)
    self.data_renderer.connect("editing-started", self.on_cell_edit_start)

#    self.data_col = data_col = gtk.TreeViewColumn("Data", self.data_renderer, text=4)
#    data_col.set_property("expand", True)
#    data_col.set_property("sizing", gtk.TREE_VIEW_COLUMN_AUTOSIZE)
#    optionsTree.append_column(data_col)

    # display name, gtk.ListStore containing the display names of possible choices, pointer to node in self.tree -- a choice or a tree, pointer to currently active tree and its data.
    self.treestore = gtk.TreeStore(str, gobject.TYPE_PYOBJECT, gobject.TYPE_PYOBJECT, gobject.TYPE_PYOBJECT, str)
    self.treeview.set_model(self.treestore)
#    self.treeview.set_grid_lines(gtk.TREE_VIEW_GRID_LINES_VERTICAL)

    optionsTree.get_selection().connect("changed", self.on_select_row)
    self.treeview.get_selection().set_select_function(self.options_tree_select_func)
    self.options_tree_select_func_enabled = True
    optionsTree.connect("button_press_event", self.on_treeview_clicked)
    optionsTree.connect("row-activated", self.on_activate_row)
    cellCombo.connect("edited", self.cellcombo_edited)

    self.treeview.set_enable_search(False)

    return

  def create_liststore(self, l):
    """
    Given a list of possible choices, create the liststore for the
    gtk.CellRendererCombo that contains the names of possible choices.
    """

    liststore = gtk.ListStore(str, gobject.TYPE_PYOBJECT)
    if l.__class__ is choice.Choice:
      l = l.l # I'm really sorry about this, blame dham

    if not isinstance(l, list):
      l = [l]

    # Ignoring the numerous l's involved in getting the choices,
    # l is now a list of possible names.
    for t in l:
      name = self.get_display_name(t)
      liststore.append([name, t])

    return liststore

  def set_treestore(self, iter=None, new_tree=[], recurse=False, replace=False):
    """
    Given a list of children of a node in a treestore, stuff them in the treestore.
    """

    if replace:
      replacediter = iter
      iter = self.treestore.iter_parent(replacediter)
    else:
      self.remove_children(iter)
    
    for t in new_tree:
      if t.__class__ is tree.Tree:
        if self.choice_or_tree_is_hidden(t):
          continue

        liststore = self.create_liststore(t)

        # Convert node data, if it exists, to a string
        data = ""
        node_data = t.data
        if node_data is not None:
#          if t.__class__ is str:
#            if len(t) > 4: # Trim the string if it's long
#              node_data = t[:4] + ".."
          
          data = str(node_data)
 
        if replace:
          child_iter = self.treestore.insert_before(iter, replacediter, [self.get_display_name(t), liststore, t, t, data])
        else:
          child_iter = self.treestore.append(iter, [self.get_display_name(t), liststore, t, t, data])
        
        if recurse and t.active: self.set_treestore(child_iter, t.children, recurse)
      elif t.__class__ is choice.Choice:
        liststore = self.create_liststore(t)
        ts_choice = t.get_current_tree()
        if self.choice_or_tree_is_hidden(ts_choice):
          continue
        if replace:
          child_iter = self.treestore.insert_before(iter, replacediter, [self.get_display_name(ts_choice), liststore, t, ts_choice, ""])
        else:
          child_iter = self.treestore.append(iter, [self.get_display_name(ts_choice), liststore, t, ts_choice, ""])
        if recurse and t.active: self.set_treestore(child_iter, ts_choice.children, recurse)

    if replace:
      self.treestore.remove(replacediter)
      return child_iter
   
    return

  def expand_choice_or_tree(self, choice_or_tree):
    """
    Query the schema for what valid children can live under this node, and add
    them to the choice or tree. This recurses.
    """

    if isinstance(choice_or_tree, choice.Choice):
      for opt in choice_or_tree.choices():
        self.expand_choice_or_tree(opt)
    else:
      l = self.s.valid_children(choice_or_tree.schemaname)
      l = choice_or_tree.find_or_add(l)
      for opt in l:
        self.expand_choice_or_tree(opt)

    return

  def expand_treestore(self, iter = None):
    """
    Query the schema for what valid children can live under this node, then set the
    treestore appropriately. This recurses.
    """

    if iter is None:
      iter = self.treestore.get_iter_first()
      if iter is None:
        self.set_treestore(iter, [])
        return

    choice_or_tree, active_tree = self.treestore.get(iter, 2, 3)
    if active_tree.active is False or choice_or_tree.active is False:
      return

    l = self.s.valid_children(active_tree.schemaname)
    l = active_tree.find_or_add(l)
    self.set_treestore(iter, l)

    child_iter = self.treestore.iter_children(iter)
    while child_iter is not None:
      # fix for recursive schemata!
      child_active_tree = self.treestore.get_value(child_iter, 3)
      if child_active_tree.schemaname == active_tree.schemaname:
        debug.deprint("Warning: recursive schema elements not supported: %s" % active_tree.name)
        child_iter = self.treestore.iter_next(child_iter)
        if child_iter is None: break

      self.expand_treestore(child_iter)
      child_iter = self.treestore.iter_next(child_iter)

    return

  def remove_children(self, iter):
    """
    Delete the children of iter in the treestore.
    """

    childiter = self.treestore.iter_children(iter)
    if childiter is None: return

    result = True

    while result is True:
      result = self.treestore.remove(childiter)

    return

  def set_combobox_liststore(self, column, cellCombo, treemodel, iter, user_data=None):
    """
    This hook function sets the properties of the gtk.CellRendererCombo for each
    row. It sets up the cellcombo to use the correct liststore for its choices,
    decides whether the cellcombo should be editable or not, and sets the
    foreground colour.
    """

    liststore = self.treestore.get_value(iter, 1)
    choice_or_tree = self.treestore.get_value(iter, 2)
    active_tree = self.treestore.get_value(iter, 3)
    data = self.treestore.get_value(iter, 4)

    # set the model for the cellcombo, where it gets the possible choices for the name
    cellCombo.set_property("model", liststore)

    # set the properties: colour, etc.
    if choice_or_tree.__class__ is tree.Tree:
      cellCombo.set_property("editable", False)
    elif choice_or_tree.__class__ is choice.Choice:
      cellCombo.set_property("editable", True)

    if self.treestore_iter_is_active(iter):
      if active_tree.valid is True:
        cellCombo.set_property("foreground", "black")
      else:
        cellCombo.set_property("foreground", "blue")
    else:
        cellCombo.set_property("foreground", "gray")

    return

  def set_cellpicture_choice(self, column, cell, model, iter):
    """
    This hook function sets up the other gtk.CellRendererPixbuf, the one that gives
    the clue to the user whether this is a choice or not.
    """

    choice_or_tree = self.treestore.get_value(iter, 2)
    if choice_or_tree.__class__ is tree.Tree:
      cell.set_property("stock-id", None)
    elif choice_or_tree.__class__ is choice.Choice:
      cell.set_property("stock-id", gtk.STOCK_GO_DOWN)

    return

  def set_cellpicture_cardinality(self, column, cell, model, iter):
    """
    This hook function sets up the gtk.CellRendererPixbuf on the extreme right-hand
    side for each row; this paints a plus or minus or nothing depending on whether
    something can be added or removed or has to be there.
    """

    choice_or_tree = self.treestore.get_value(iter, 2)
    if choice_or_tree.cardinality == "":
      cell.set_property("stock-id", None)
    elif choice_or_tree.cardinality == "?" or choice_or_tree.cardinality == "*":
      if choice_or_tree.active:
        cell.set_property("stock-id", gtk.STOCK_REMOVE)
      else:
        cell.set_property("stock-id", gtk.STOCK_ADD)
    elif choice_or_tree.cardinality == "+":
      parent_tree = choice_or_tree.parent
      count = parent_tree.count_children_by_schemaname(choice_or_tree.schemaname)
      
      if choice_or_tree.active and count == 2: # one active, one inactive
        cell.set_property("stock-id", None)
      elif choice_or_tree.active:
        cell.set_property("stock-id", gtk.STOCK_REMOVE)
      else:
        cell.set_property("stock-id", gtk.STOCK_ADD)

    return

  def on_treeview_row_collapsed(self, treeview, iter, path):
    """
    Called when a row in the LHS treeview is collapsed.
    """

    #self.treeview.get_column(0).queue_resize()
    #self.treeview.get_column(1).queue_resize()

    return

  def on_tooltip(self, widget, x, y, keyboard_mode, tooltip):
    y-=25 # It's hardcoded. Gtk doesn't offer a way to get at the column height without a lot more code.
    (tx, ty) = self.treeview.convert_bin_window_to_tree_coords(x, y)
    pathinfo = self.treeview.get_path_at_pos(x, y)

    if pathinfo is None:
      return False

    path = pathinfo[0]
    column = pathinfo[1]

    if path is None:
      return

    ctitle = column.get_title()

    iter = self.treestore.get_iter(path)

    # Get the tree or choice pointed to by the iterator.
    tree = self.treestore.get_value(iter, 2)

    if tree.__class__ is choice.Choice:
      tree = tree.get_current_tree()

    if ctitle == "Node":
      if tree.doc is None:
        text = "(No documentation)"
      else:
        text = self.render_whitespace(tree.doc)
    elif ctitle == "Data":
      comment_tree = self.get_comment(tree)

      if comment_tree is None:
        text = "(No comment)"
      else:
        text = comment_tree.data
    else:
      return False

    tooltip.set_text(text)    
    return True

  def on_cell_edit(self, cellr, path, new_text):
    iter = self.treestore.get_iter(path)
    node = self.selected_node

    if new_text == "":
      return

    # mtw07 - Check whether the data is valid or not. If so, update the treeview
    # and internal model. GTK automatically clears the text field if we don't
    # add the new data to the treestore.
    (invalid, data) = node.valid_data(node.datatype, new_text)

    if invalid:
      try:
        name = node.datatype.__name__
        dialogs.error(self.main_window, "Invalid data type. A %s is required." % name)
      except:
        pass
      return

    self.treestore.set_value(iter, 4, new_text)
    node.data = new_text

    # Update the node data box.
    self.node_data.get_buffer().set_text(new_text)

    # Update the validation errors list.
    if self.scherror.errlist_is_open():
      if self.scherror.errlist_type == 0:
         self.scherror.on_validate_schematron()
      else:
         self.scherror.on_validate()

    return

  def on_cell_edit_start(self, cellrenderer, editable, path):
    iter = self.treestore.get_iter(path)
    node = self.selected_node

    # If the cell has (<type>) in it, clear it for the user to enter data.
    if node.datatype is not None and node.data is None:
      editable.set_text("")

    return True

  def on_treeview_clicked(self, treeview, event):
    """
    This routine is called every time the mouse is clicked on the treeview on the
    left-hand side. It processes the "buttons" gtk.STOCK_ADD and gtk.STOCK_REMOVE
    in the right-hand column, activating, adding and removing tree nodes as
    necessary.
    """

    if event.button != 1:
      return

    pathinfo = treeview.get_path_at_pos(int(event.x), int(event.y))

    if pathinfo is None:
      return

    path = pathinfo[0]
    col = pathinfo[1]

    iter = self.treestore.get_iter(path)
    self.update_data_column(self.treestore, iter)

    if col is not self.imgcolumn:
      return

    self.toggle_tree(iter)

    self.on_select_row(self.treeview.get_selection())

    self.treeview.queue_draw()
    self.treeview.get_column(0).queue_resize()
    self.treeview.get_column(1).queue_resize()

    return

  def toggle_tree(self, iter):
    """
    Toggles the state of part of the tree.
    """

    choice_or_tree = self.treestore.get_value(iter, 2)

    if choice_or_tree.active:
      self.collapse_tree(iter)
    else:
      self.expand_tree(iter)

    return
  
  def collapse_tree(self, iter):
    """
    Collapses part of the tree.
    """

    choice_or_tree, = self.treestore.get(iter, 2)
    parent_iter = self.treestore.iter_parent(iter)

    if parent_iter == None:
      parent_tree = None
    else:
      parent_tree = self.treestore.get_value(parent_iter, 3)

    if not choice_or_tree.active:
      return

    if choice_or_tree.cardinality == "":
      return

    if choice_or_tree.cardinality == "?":
      choice_or_tree.active = False
      self.set_saved(False)
      self.remove_children(iter)

    elif choice_or_tree.cardinality == "*":
      # If this is the only one, just make it inactive.
      # Otherwise, just delete it.
      count = parent_tree.count_children_by_schemaname(choice_or_tree.schemaname)
      if count == 1:
        choice_or_tree.active = False
        self.set_saved(False)
        self.remove_children(iter)
      else:
        self.delete_tree(iter)

    elif choice_or_tree.cardinality == "+":
      count = parent_tree.count_children_by_schemaname(choice_or_tree.schemaname)
      if count == 2: # one active, one inactive
        # do nothing
        return
      else: # count > 2
        self.delete_tree(iter)
  
    parent_tree.recompute_validity()
    self.treeview.queue_draw()
    return

  def delete_tree(self, iter):
    choice_or_tree, = self.treestore.get(iter, 2)
    parent_iter = self.treestore.iter_parent(iter)
    isSelected = self.treeview.get_selection().iter_is_selected(iter)
    sibling = self.treestore.iter_next(iter)

    if parent_iter == None:
      parent_tree = None
    else:
      parent_tree = self.treestore.get_value(parent_iter, 3)

    confirm = dialogs.prompt(self.main_window, "Are you sure you want to delete this node?")
    if confirm == gtk.RESPONSE_YES:
      parent_tree.delete_child_by_ref(choice_or_tree)
      self.remove_children(iter)
      self.treestore.remove(iter)
      self.set_saved(False)
      
      if isSelected and sibling:
        self.treeview.get_selection().select_iter(sibling)
    return

  def expand_tree(self, iter):
    """
    Expands part of the tree.
    """

    choice_or_tree, active_tree = self.treestore.get(iter, 2, 3)
    parent_iter = self.treestore.iter_parent(iter)

    if parent_iter == None:
      parent_tree = None
    else:
      parent_tree = self.treestore.get_value(parent_iter, 3)

    if choice_or_tree.active:
      return

    if choice_or_tree.cardinality == "":
      return

    elif choice_or_tree.cardinality == "?":
      choice_or_tree.active = True
      self.set_saved(False)
      self.expand_treestore(iter)

    elif choice_or_tree.cardinality == "*" or choice_or_tree.cardinality == "+":
      # Make this active, and add a new inactive instance
      choice_or_tree.active = True
      new_tree = parent_tree.add_inactive_instance(choice_or_tree)
      liststore = self.create_liststore(new_tree)
      self.expand_treestore(iter)
      iter = self.treestore.insert_after(
        parent=parent_iter, sibling=iter, 
        row=[self.get_display_name(new_tree), liststore, new_tree, new_tree.get_current_tree(), ""])
      self.set_saved(False)

    parent_tree.recompute_validity()
    return

  def options_tree_select_func(self, info = None):
    """
    Called when the user selected a new item in the treeview. Prevents changing of
    node and attempts to save data if appropriate.
    """

    if not self.options_tree_select_func_enabled:
      self.options_tree_select_func_enabled = True
      return False

    if not self.node_data_store():
      return False

    if isinstance(self.selected_node, MixedTree) and not self.geometry_dim_tree is None and self.selected_node.parent is self.geometry_dim_tree.parent and not self.selected_node.data is None:
      self.geometry_dim_tree.set_data(self.selected_node.data)

    return True

  def on_treeview_key_press(self, treeview, event):
    """
    Called when treeview intercepts a key press. Collapse and expand rows.
    """

    if event.keyval == gtk.keysyms.Right:
      self.treeview.expand_row(self.get_selected_row(), open_all = False)

    if event.keyval == gtk.keysyms.Left:
      self.treeview.collapse_row(self.get_selected_row())

    if event.keyval == gtk.keysyms.Delete:
       self.collapse_tree(self.treestore.get_iter(self.get_selected_row()))
 
    return

  def on_treeview_button_press(self, treeview, event):
    if event.button == 3:
      x = int(event.x)
      y = int(event.y)
      path = treeview.get_path_at_pos(x, y)[0]
      if path is not None:
        treeview.get_selection().select_path(path)
        self.show_popup(None, event.button, event.time)
        return True
    return False

  def popup_location(self, widget, user_data):
    column = self.treeview.get_column(0)
    path = self.treeview.get_selection().get_selected_rows()[0]
    area = self.treeview.get_cell_area(path, column)
    tx, ty = area.x, area.y
    x, y = self.treeview.tree_to_widget_coords(tx, ty)
    return (x, y, True)
    
  def on_treeview_popup(self, treeview):
    self.show_popup(None, self.popup_location, gtk.get_current_event_time())
    return

  def show_popup(self, func, button, time):
    self.popup.popup( None, None, func, button, time)  
    return

  def on_select_row(self, selection=None):
    """
    Called when a row is selected. Update the options frame.
    """
    
    path = self.get_selected_row(self.treeview.get_selection())
    if path is None:
      return
    self.selected_iter = iter = self.treestore.get_iter(path)
    choice_or_tree = self.treestore.get_value(iter, 2)

    active_tree = self.treestore.get_value(iter, 3)
    debug.dprint(active_tree)

    self.selected_node = self.get_painted_tree(iter)
    self.update_options_frame()

    node = self.selected_node

    if isinstance(node, MixedTree):
      node = node.child

    # TODO: Handle tuple datatypes. Is it possible in a Gtk treeview column?
    if not isinstance(node.datatype, tuple):
      self.data_renderer.set_property("editable", not (node.not_editable()))

    name = self.get_spudpath(active_tree)
    self.statusbar.set_statusbar(name)
    self.current_spudpath = name
    self.current_xpath = self.get_xpath(active_tree)

    self.clear_plugin_buttons()

    for plugin in plugins.plugins:
      if plugin.matches(name):
        self.add_plugin_button(plugin)

    return

  def get_spudpath(self, active_tree):
    # get the name to paint on the statusbar
    name_tree = active_tree
    name = ""
    while name_tree is not None:
      if "name" in name_tree.attrs and name_tree.attrs["name"][1] is not None:
        used_name = name_tree.name + '::%s' % name_tree.attrs["name"][1]
      elif name_tree.parent is not None and name_tree.parent.count_children_by_schemaname(name_tree.schemaname) > 1:
        siblings = [x for x in name_tree.parent.children if x.schemaname == name_tree.schemaname]
        i = 0
        for sibling in siblings:
          if sibling is name_tree:
            break
          else:
            i = i + 1
        used_name = name_tree.name + "[%s]" % i
      else:
        used_name = name_tree.name

      name = "/" + used_name + name
      name_tree = name_tree.parent

    # and split off the root name:
    name = '/' + '/'.join(name.split('/')[2:])
    return name

  def get_xpath(self, active_tree):
    # get the name to paint on the statusbar
    name_tree = active_tree
    name = ""
    while name_tree is not None:
      if "name" in name_tree.attrs and name_tree.attrs["name"][1] is not None:
        used_name = name_tree.name + '[@name="%s"]' % name_tree.attrs["name"][1]
      elif name_tree.parent is not None and name_tree.parent.count_children_by_schemaname(name_tree.schemaname) > 1:
        siblings = [x for x in name_tree.parent.children if x.schemaname == name_tree.schemaname]
        i = 0
        for sibling in siblings:
          if sibling is name_tree:
            break
          else:
            i = i + 1
        used_name = name_tree.name + "[%s]" % i
      else:
        used_name = name_tree.name

      name = "/" + used_name + name
      name_tree = name_tree.parent

    return name

  def clear_plugin_buttons(self):
    for button in self.plugin_buttons:
      self.plugin_buttonbox.remove(button)
    
    self.plugin_buttons = []

  def add_plugin_button(self, plugin):
    button = gtk.Button(label=plugin.name)
    button.connect('clicked', self.plugin_handler, plugin)
    button.show()

    self.plugin_buttons.append(button)
    self.plugin_buttonbox.add(button)

  def plugin_handler(self, widget, plugin):
    f = StringIO.StringIO()
    self.tree.write(f)
    xml = f.getvalue()
    plugin.execute(xml, self.current_xpath)

  def get_selected_row(self, selection=None):
    """
    Get the iter to the selected row.
    """

    if (selection == None):
        selection = self.gui.get_widget("optionsTree").get_selection()

    (model, paths) = selection.get_selected_rows()
    if ((len(paths) != 1) or (paths[0] == None)):
      return None
    else:
      return paths[0]

  def update_data_column(self, model, itParent):
    """
    Update the data column in the treeview. Used when a user expands a row.
    """

    iter = model.iter_children(itParent)

    while iter is not None:
      self.set_cell_node(iter)
      iter = model.iter_next(iter)
   

  def on_activate_row(self, treeview, path, view_column):
    """
    Called when you double click or press Enter on a row.
    """

    iter = self.treestore.get_iter(path)
    
    self.expand_tree(iter)

    if path is None: 
      return
    
    if treeview.row_expanded(path):
      treeview.collapse_row(path)
    else:
      treeview.expand_row(path, False)

      # Update the data column for the newly expanded row.
      (model, itParent) = self.treeview.get_selection().get_selected()
      self.update_data_column(model, itParent)

    return

  def cellcombo_edited(self, cellrenderertext, path, new_text):
    """
    This is called when a cellcombo on the left-hand treeview is edited,
    i.e. the user chooses between more than one possible choice.
    """

    iter = self.treestore.get_iter(path)
    self.treestore.set(iter, 0, new_text)
    choice = self.treestore.get_value(iter, 2)

    # get the ref to the new active choice
    liststore = self.treestore.get_value(iter, 1)
    list_iter = liststore.get_iter_first()
    ref = None
    while list_iter is not None:
      list_text = liststore.get_value(list_iter, 0)
      if list_text == new_text:
        ref = liststore.get_value(list_iter, 1)
        break
      list_iter = liststore.iter_next(list_iter)

    # record the choice in the datatree
    choice.set_active_choice_by_ref(ref)
    new_active_tree = choice.get_current_tree()

    name = self.get_spudpath(new_active_tree)
    self.statusbar.set_statusbar(name)
    self.treestore.set(iter, 3, new_active_tree)
    self.current_spudpath = name
    xpath = self.get_xpath(new_active_tree)
    self.current_xpath = xpath

    self.clear_plugin_buttons()

    for plugin in plugins.plugins:
      if plugin.matches(xpath):
        self.add_plugin_button(plugin)

    self.remove_children(iter)
    self.expand_treestore(iter)
    self.treeview.expand_row(path, False)

    self.set_saved(False)
    self.selected_node = self.get_painted_tree(iter)
    self.update_options_frame()

    return

  def paint_validity(self):
    """
    Walk up the parental line, repainting the colour for their validity
    appropriately. This is called when a validity-changing event occurs.
    Repaint the whole tree if there is no selection.
    """

    def paint_iter_validity(iter):
      while iter is not None:
        active_tree = self.treestore.get_value(iter, 3)
        if active_tree.valid:
          self.cellcombo.set_property("foreground", "black")
        else:
          self.cellcombo.set_property("foreground", "blue")
        iter = self.treestore.iter_next(iter)

      return

    selection = self.treeview.get_selection()
    path = self.get_selected_row(selection)
    if path is None:
      paint_iter_validity(self.treestore.get_iter_first())
    else:
      paint_iter_validity(self.treestore.get_iter(path))

    self.treeview.queue_draw()

    return

  def get_display_name(self, active_tree):
    """
    This is a fluidity hack, allowing the name displayed in the treeview on the
    left to be different to the element name. If it has an attribute name="xxx",
    element_tag (xxx) is displayed.
    """

    if active_tree.__class__ is tree.Tree:
      displayname = active_tree.name
      if "name" in active_tree.attrs.keys():
        attrname = active_tree.attrs["name"][1]
        if attrname is not None:
          displayname = displayname + " (" + attrname + ")"
    elif active_tree.__class__ is choice.Choice:
      displayname = self.get_display_name(active_tree.get_current_tree())

    return displayname

  def get_treeview_iter(self, selection):
    """
    Get a treeview iterator object, given a selection.
    """

    path = self.get_selected_row(selection)
    if path is None:
      return self.get_treestore_iter_from_xmlpath(self.current_xpath)

    return self.treestore.get_iter(path)

  def update_painted_name(self):
    """
    This updates the treestore (and the liststore for the gtk.CellRendererCombo)
    with a new name, when the name="xxx" attribute is changed.
    """

    iter = self.get_treeview_iter(self.treeview.get_selection())
    if iter is None:
      return

    liststore = self.treestore.get_value(iter, 1)
    active_tree = self.treestore.get_value(iter, 3)
    new_name = self.get_display_name(active_tree)
    self.treestore.set_value(iter, 0, new_name)

    # find the liststore iter corresponding to the painted choice
    list_iter = liststore.get_iter_first()
    while list_iter is not None:
      liststore_tree = liststore.get_value(list_iter, 1)
      if liststore_tree is active_tree:
        liststore.set_value(list_iter, 0, new_name)
      list_iter = liststore.iter_next(list_iter)

    self.treeview.get_column(0).queue_resize()

    return

  def get_painted_tree(self, iter_or_tree, lock_geometry_dim = True):
    """
    Check if the given tree, or the active tree at the given iter in the treestore,
    have any children of the form *_value. If so, we need to make the node painted
    by the options tree a mix of the two: the documentation and attributes come from
    the parent, and the data from the child.

    Also check if it is the geometry node, validity of any tuple data, and, if an
    iter is supplied, check that the node is active.
    """

    if isinstance(iter_or_tree, tree.Tree):
      active_tree = iter_or_tree
    else:
      active_tree = self.treestore.get_value(iter_or_tree, 3)

    integers = [child for child in active_tree.children if child.name == "integer_value"]
    reals    = [child for child in active_tree.children if child.name == "real_value"]
    logicals = [child for child in active_tree.children if child.name == "logical_value"]
    strings  = [child for child in active_tree.children if child.name == "string_value"]

    child = None
    if len(integers) > 0:
      child = integers[0]
    if len(reals) > 0:
      child = reals[0]
    if len(logicals) > 0:
      child = logicals[0]
    if len(strings) > 0:
      child = strings[0]

    if child is None:
      painted_tree = active_tree
    else:
      painted_tree = MixedTree(active_tree, child)

    if not isinstance(iter_or_tree, tree.Tree) and not self.treestore_iter_is_active(iter_or_tree):
      painted_tree = tree.Tree(painted_tree.name, painted_tree.schemaname, painted_tree.attrs, doc = painted_tree.doc)
      painted_tree.active = False
    elif lock_geometry_dim and not self.geometry_dim_tree is None and not self.geometry_dim_tree.data is None:
      if active_tree is self.geometry_dim_tree:
        data_tree = tree.Tree(painted_tree.name, painted_tree.schemaname, datatype = "fixed")
        data_tree.data = painted_tree.data
        painted_tree = MixedTree(painted_tree, data_tree)
      elif isinstance(self.geometry_dim_tree, MixedTree) and active_tree is self.geometry_dim_tree.parent:
        data_tree = tree.Tree(painted_tree.child.name, painted_tree.child.schemaname, datatype = "fixed")
        data_tree.data = painted_tree.data
        painted_tree = MixedTree(painted_tree, data_tree)

    return painted_tree

  def get_treestore_iter_from_xmlpath(self, xmlpath):
    """
    Convert the given XML path to an iter into the treestore. For children of a
    single parent with the same names, only the first child is considered.
    """

    names = xmlpath.split("/")

    iter = self.treestore.get_iter_first()
    if iter is None:
      return None
    for name in names[1:len(names) - 1]:
      while not self.treestore.get_value(iter, 0) == name:
        iter = self.treestore.iter_next(iter)
        if iter is None:
          return None
      iter = self.treestore.iter_children(iter)
      if iter is None:
        return None

    return iter
    
  def set_geometry_dim_tree(self):
    """
    Find the iter into the treestore corresponding to the geometry dimension, and
    perform checks to test that the geometry dimension node is valid.
    """

    # The tree must exist
    if self.tree is None:
      self.geometry_dim_tree = None
      return

    # A geometry dimension element must exist
    iter = self.get_treestore_iter_from_xmlpath("/" + self.tree.name + self.data_paths["dim"])
    if iter is None:
      self.geometry_dim_tree = None
      return

    painted_tree = self.get_painted_tree(iter, False)
    if isinstance(painted_tree, MixedTree):
       # If the geometry dimension element has a hidden data element, it must
       # have datatype tuple or fixed
       if not isinstance(painted_tree.datatype, tuple) and not painted_tree.datatype == "fixed":
         self.geometry_dim_tree = None
         return
    elif not painted_tree.datatype == "fixed":
      # Otherwise, only fixed datatype is permitted
      self.geometry_dim_tree = None
      return

    # All parents of the geometry dimension element must have cardinality ""
    # (i.e. not ?, * or +).
    parent = painted_tree.parent
    while not parent is None:
      if not parent.cardinality == "":
        self.geometry_dim_tree = None
        return
      parent = parent.parent

    # All possible geometry dimensions must be positive integers
    if isinstance(painted_tree.datatype, tuple):
      possible_dims = painted_tree.datatype
    elif painted_tree.datatype == "fixed":
      possible_dims = [painted_tree.data]
    else:
      self.geometry_dim_tree = None
      return
    for opt in possible_dims:
      try:
        test = int(opt)
        assert test > 0
      except:
        self.geometry_dim_tree = None
        return
      
    # A valid geometry dimension element has been located
    self.geometry_dim_tree = painted_tree
    
    return

  def treestore_iter_is_active(self, iter):
    """
    Test whether the node at the given iter in the LHS treestore is active.
    """

    while not iter is None:
      choice_or_tree = self.treestore.get_value(iter, 2)
      active_tree = self.treestore.get_value(iter, 3)
      if not choice_or_tree.active or not active_tree.active:
        return False
      iter = self.treestore.iter_parent(iter)

    return True

  def choice_or_tree_is_hidden(self, choice_or_tree):
    """
    Tests whether the supplied choice or tree should be hidden from the LHS.
    """

    return self.choice_or_tree_is_comment(choice_or_tree) or choice_or_tree.name in ["integer_value", "real_value", "string_value", "logical_value"]

  def choice_or_tree_is_comment(self, choice_or_tree):
    """
    Test whether the given node is a comment node.
    """

    if not isinstance(choice_or_tree, tree.Tree):
      return False

    if not choice_or_tree.name == "comment":
      return False

    if not choice_or_tree.attrs == {}:
      return False

    if not choice_or_tree.children == []:
      return False

    if not choice_or_tree.datatype is str:
      return False

    if not choice_or_tree.cardinality == "?":
      return False

    return True

  def get_comment(self, choice_or_tree):
    """
    Return the first comment found as a child of the supplied node, or None if
    none found.
    """

    if choice_or_tree is None or isinstance(choice_or_tree, choice.Choice):
      return None

    for child in choice_or_tree.children:
      if self.choice_or_tree_is_comment(child):
        return child

    return None

    return

  def choice_or_tree_matches(self, text, choice_or_tree, recurse, search_active_subtrees = False):
    """
    See if the supplied node matches a given piece of text. If recurse is True,
    the node is deemed to match if any of its children match or, if the node is a
    choice and search_active_subtrees is True, any of the available trees in the
    choice match.
    """

    if self.choice_or_tree_is_hidden(choice_or_tree):
      return False
    elif isinstance(choice_or_tree, choice.Choice):
      if self.choice_or_tree_matches(text, choice_or_tree.get_current_tree(), False):
        return True
      elif recurse and self.find.search_gui.get_widget("searchInactiveChoiceSubtreesCheckButton").get_active():
        for opt in choice_or_tree.choices():
          if not search_active_subtrees and opt is choice_or_tree.get_current_tree():
            continue
          if opt.children == []:
            self.expand_choice_or_tree(opt)
          if self.choice_or_tree_matches(text, opt, recurse, True):
            return True
    else:
      if self.get_painted_tree(choice_or_tree).matches(text, self.find.search_gui.get_widget("caseSensitiveCheckButton").get_active()):
        return True
      else:
        if self.find.search_gui.get_widget("caseSensitiveCheckButton").get_active():
          text_re = re.compile(text)
        else:
          text_re = re.compile(text, re.IGNORECASE)
        comment = self.get_comment(choice_or_tree)
        if not comment is None and not comment.data is None and not text_re.search(comment.data) is None:
          return True
        elif recurse:
          for opt in choice_or_tree.children:
            if self.choice_or_tree_matches(text, opt, recurse, True):
              return True

      return False

  def search_treestore(self, text, iter = None):
    """
    Recursively search the tree for a node that matches a given piece of text.
    MixedTree.matches and choice_or_tree_matches decide what is a match (using
    tree.Matches).

    This uses lazy evaluation to only search as far as necessary; I love
    Python generators. If you don't know what a Python generator is and need to
    understand this, see PEP 255.
    """

    if iter is None:
      iter = self.treestore.get_iter_first()
    if iter is None:
      yield None
    choice_or_tree = self.treestore.get_value(iter, 2)

    if self.choice_or_tree_matches(text, choice_or_tree, isinstance(choice_or_tree, choice.Choice)):
      yield iter

    child_iter = self.treestore.iter_children(iter)
    while child_iter is not None:
      for iter in self.search_treestore(text, child_iter):
        yield iter
      child_iter = self.treestore.iter_next(child_iter)

    return

  ### RHS ###
  
  def render_whitespace(self, desc):
    ''' Render the line wrapping in desc as follows:
    
    * Newlines followed by 0-1 spaces are ignored.
    * Blank lines start new paragraphs.
    * Newlines followed by more than 1 space are honoured.
    '''

    prev_line_literal=False
    prev_line_new_para=False
    newdesc=""

    for line in desc.split("\n"):

        if (line[:1]==" "):
            # Literal line with leading blanks.
            newdesc=newdesc+"\n"+line
            prev_line_literal=True
            prev_line_new_para=False
            continue

        if (line.strip()==""):
            # New paragraph.

            # Collapse multiple new paragraphs into one, except when
            # following a literal line.
            if (prev_line_new_para and not prev_line_literal):
                continue

            newdesc=newdesc+"\n"
            prev_line_new_para=True                
            continue

        if prev_line_literal:
            newdesc=newdesc+"\n"
            prev_line_literal=False
            prev_line_new_para=False
            
        if prev_line_new_para:
            newdesc=newdesc+"   "
            prev_line_new_para=False
        
        # Default case
        newdesc=newdesc+line+" "
    
    return newdesc

  def link_bounds(self, text):
    """
    Return a list of tuples corresponding to the start and end points of links in
    the supplied string.
    """

    bounds = []

    text_split = text.lower().split("http://")
    if len(text_split) > 1:
      lbound = -7
      for i in range(len(text_split))[1:]:
        lbound += len(text_split[i - 1]) + 7
        ubound = lbound + len(text_split[i].split(" ")[0].split("\n")[0]) + 7
        while text[ubound - 1:ubound] in [".", ",", ":", ";", "\"", "'", ")", "]", "}"]:
          ubound -= 1
        bounds.append((lbound, ubound))

    return bounds

  def type_name(self, datatype):
    """
    Return a human readable version of datatype.
    """

    def familiar_type(type_as_printable):
      """
      Convert some type names to more familiar equivalents.
      """

      if type_as_printable == "decim":
        return "float"
      elif type_as_printable == "int":
        return "integer"
      elif type_as_printable == "str":
        return "string"
      else:
        return type_as_printable

    datatype_string = str(datatype)

    if datatype_string[:7] == "<type '" and datatype_string[len(datatype_string) - 2:] == "'>":
      value_type_split = datatype_string.split("'")
      return familiar_type(value_type_split[1])

    value_type_split1 = datatype_string.split(".")
    value_type_split2 = value_type_split1[len(value_type_split1) - 1].split(" ")
    if len(value_type_split2) == 1:
      return familiar_type(value_type_split2[0][0:len(value_type_split2[0]) - 6])
    else:
      return familiar_type(value_type_split2[0])

  def printable_type(self, datatype, bracket = True):
    """
    Create a string to be displayed in place of empty data / attributes.
    """

    if isinstance(datatype, plist.List):
      if (isinstance(datatype.cardinality, int) and datatype.cardinality == 1) or datatype.cardinality == "":
        type_as_printable = self.type_name(datatype.datatype).lower()
      else:
        type_as_printable = self.type_name(datatype).lower() + " of "
        list_type_as_printable = self.type_name(datatype.datatype).lower()
        if isinstance(datatype.cardinality, int):
          type_as_printable += str(datatype.cardinality) + " " + list_type_as_printable + "s"
        else:
          type_as_printable += list_type_as_printable + "s"
    else:
      type_as_printable = self.type_name(datatype).lower()

    if bracket:
      type_as_printable = "(" + type_as_printable + ")"

    return type_as_printable

  def init_options_frame(self):
    """
    Initialise the RHS.
    """

    # mtw07 - TODO: set the visibility of the options frame, depending on the value in the schema file.
    self.options_frame = self.gui.get_widget("optionsFrame")

    # Display the right hand side by default.
    self.options_frame.set_property("visible", True)

    self.node_desc = self.gui.get_widget("nodeDescription")
    self.node_desc.set_buffer(TextBufferMarkup.PangoBuffer())
    self.node_desc.connect("button-release-event", self.node_desc_mouse_button_release)
    self.node_desc.connect("motion-notify-event", self.node_desc_mouse_over)

    self.node_attrs = self.gui.get_widget("nodeAttributes")
    attrs_model = gtk.ListStore(gobject.TYPE_STRING, gobject.TYPE_STRING, gobject.TYPE_PYOBJECT)
    self.node_attrs.set_model(attrs_model)
    self.node_attrs.connect("motion-notify-event", self.node_attrs_mouse_over)
    key_renderer = gtk.CellRendererText()
    key_renderer.set_property("editable", False)
    attrs_col1 = gtk.TreeViewColumn("Name", key_renderer, text = 0)
    attrs_col1.set_cell_data_func(key_renderer, self.node_attrs_key_data_func)
    attrs_col1.set_property("min-width", 75)
    attrs_val_entry_renderer = gtk.CellRendererText()
    attrs_val_entry_renderer.connect("edited", self.node_attrs_edited)
    attrs_val_entry_renderer.connect("editing-started", self.node_attrs_entry_edit_start)
    attrs_val_combo_renderer = gtk.CellRendererCombo()
    attrs_val_combo_renderer.set_property("text-column", 0)
    attrs_val_combo_renderer.connect("edited", self.node_attrs_selected)
    attrs_val_combo_renderer.connect("editing-started", self.node_attrs_combo_edit_start)
    attrs_col2 = gtk.TreeViewColumn("Value", attrs_val_entry_renderer, text = 1)
    attrs_col2.pack_start(attrs_val_combo_renderer)
    attrs_col2.set_attributes(attrs_val_combo_renderer, text = 1)
    attrs_col2.set_cell_data_func(attrs_val_entry_renderer, self.node_attrs_entry_data_func)
    attrs_col2.set_cell_data_func(attrs_val_combo_renderer, self.node_attrs_combo_data_func)
    attrs_col2.set_property("expand", True)
    attrs_col2.set_property("min-width", 75)
    attrs_icon_renderer = gtk.CellRendererPixbuf()
    attrs_col3 = gtk.TreeViewColumn("", attrs_icon_renderer)
    attrs_col3.set_cell_data_func(attrs_icon_renderer, self.node_attrs_icon_data_func)
    self.node_attrs.append_column(attrs_col1)
    self.node_attrs.append_column(attrs_col2)
    self.node_attrs.append_column(attrs_col3)

    self.node_data_frame = self.gui.get_widget("dataFrame")

    self.node_data_buttons_hbox = self.gui.get_widget("dataButtonsHBox")

    data_revert_button = self.gui.get_widget("dataRevertButton")
    data_revert_button.connect("clicked", self.node_data_revert)

    data_store_button = self.gui.get_widget("dataStoreButton")
    data_store_button.connect("clicked", self.node_data_store)

    self.node_comment = self.gui.get_widget("nodeComment")
    self.node_comment.get_buffer().create_tag("comment_buffer_tag")
    self.node_comment.connect("focus-in-event", self.node_comment_focus_in)
    self.node_comment.connect("expose-event", self.node_comment_expose)

    return

  def update_options_frame(self):
    """
    Update the RHS.
    """

    if self.selected_node is None:
      self.set_node_desc("<span foreground=\"grey\">No node selected</span>")
    elif self.selected_node.doc is None:
      self.set_node_desc("<span foreground=\"red\">No documentation</span>")
    else:
      self.set_node_desc(self.selected_node.doc)

    self.update_node_attrs()

    if self.selected_node is None or not self.selected_node.active:
      self.set_node_data_entry()
    elif self.node_data_is_tensor() and not self.geometry_dim_tree.data is None:
      self.set_node_data_tensor()
    elif isinstance(self.selected_node.datatype, tuple):
      self.set_node_data_combo()
    else:
      self.set_node_data_entry()

    self.update_node_comment()

    self.gui.get_widget("optionsFrame").queue_resize()
    
    return

  def node_desc_mouse_over(self, widget, event):
    """
    Called when the mouse moves over the node description widget. Sets the cursor
    to a hand if the mouse hovers over a link.

    Based on code from HyperTextDemo class in hypertext.py from PyGTK 2.12 demos
    """

    if self.selected_node is None or self.selected_node.doc is None:
      return

    buffer_pos = self.node_desc.window_to_buffer_coords(gtk.TEXT_WINDOW_TEXT, int(event.x), int(event.y))
    char_offset = self.node_desc.get_iter_at_location(buffer_pos[0], buffer_pos[1]).get_offset()

    for bounds in self.node_desc_link_bounds:
      if char_offset >= bounds[0] and char_offset <= bounds[1]:
        self.node_desc.get_window(gtk.TEXT_WINDOW_TEXT).set_cursor(gtk.gdk.Cursor(gtk.gdk.HAND2))
        return

    self.node_desc.get_window(gtk.TEXT_WINDOW_TEXT).set_cursor(gtk.gdk.Cursor(gtk.gdk.XTERM))

    return

  def node_desc_mouse_button_release(self, widget, event):
    """
    Called when a mouse button is released over the node description widget.
    Launches a browser if the mouse release was over a link, the left mouse button
    was released and no text was selected.

    Based on code from HyperTextDemo class in hypertext.py from PyGTK 2.12 demos
    """

    if self.selected_node is None or self.selected_node.doc is None or not event.button == 1:
      return

    selection_bounds = self.node_desc.get_buffer().get_selection_bounds()
    if not len(selection_bounds) == 0:
      return

    buffer_pos = self.node_desc.window_to_buffer_coords(gtk.TEXT_WINDOW_TEXT, int(event.x), int(event.y))
    char_offset = self.node_desc.get_iter_at_location(buffer_pos[0], buffer_pos[1]).get_offset()

    for bounds in self.node_desc_link_bounds:
      if char_offset >= bounds[0] and char_offset <= bounds[1]:
        webbrowser.open(self.selected_node.doc[bounds[0]:bounds[1]])
        return

    return

  def set_node_desc(self, desc):
    """
    Set the node description.
    """

    self.node_desc_link_bounds = self.link_bounds(desc)

    desc = self.render_whitespace(desc)

    if not len(self.node_desc_link_bounds) == 0:
      new_desc = ""
      for i in range(len(self.node_desc_link_bounds)):
        if i == 0:
          new_desc += desc[:self.node_desc_link_bounds[i][0]]
        else:
          new_desc += desc[self.node_desc_link_bounds[i - 1][1]:self.node_desc_link_bounds[i][0]]
        new_desc += "<span foreground=\"blue\" underline=\"single\">" + desc[self.node_desc_link_bounds[i][0]:self.node_desc_link_bounds[i][1]] + "</span>"
      new_desc += desc[self.node_desc_link_bounds[len(self.node_desc_link_bounds) - 1][1]:]
      if self.node_desc_link_bounds[len(self.node_desc_link_bounds) - 1][1] == len(desc):
        new_desc += " "
      desc = new_desc

    self.node_desc.get_buffer().set_text(desc)

    return

  def update_node_attrs(self):
    """
    Update the RHS attributes widget.
    """

    self.node_attrs.get_model().clear()

    if self.selected_node is None:
      self.node_attrs.get_column(2).set_property("visible", False)
      self.node_attrs.get_column(0).queue_resize()
      self.node_attrs.get_column(1).queue_resize()
    elif len(self.selected_node.attrs.keys()) == 0:
      self.gui.get_widget("attributeFrame").set_property("visible", False)
    else:
      self.gui.get_widget("attributeFrame").set_property("visible", True)
      for key in self.selected_node.attrs.keys():
        iter = self.node_attrs.get_model().append()
        self.node_attrs.get_model().set_value(iter, 0, key)
        cell_model = gtk.ListStore(gobject.TYPE_STRING)
        self.node_attrs.get_model().set_value(iter, 2, cell_model)
        if isinstance(self.selected_node.attrs[key][0], tuple):
          if self.selected_node.attrs[key][1] is None:
            if isinstance(self.selected_node.attrs[key][0][0], tuple):
              self.node_attrs.get_model().set_value(iter, 1, "Select " + self.printable_type(self.selected_node.attrs[key][0][1]) + "...")
            else:
              self.node_attrs.get_model().set_value(iter, 1, "Select...")
          else:
            self.node_attrs.get_model().set_value(iter, 1, self.selected_node.attrs[key][1])
          if isinstance(self.selected_node.attrs[key][0][0], tuple):
            opts = self.selected_node.attrs[key][0][0]
          else:
            opts = self.selected_node.attrs[key][0]
          for opt in opts:
            cell_iter = cell_model.append()
            cell_model.set_value(cell_iter, 0, opt)
          self.node_attrs.get_column(2).set_property("visible", True)
        elif self.selected_node.attrs[key][0] is None:
          self.node_attrs.get_model().set_value(iter, 1, "No data")
        elif self.selected_node.attrs[key][1] is None:
          self.node_attrs.get_model().set_value(iter, 1, self.printable_type(self.selected_node.attrs[key][0]))
        else:
          self.node_attrs.get_model().set_value(iter, 1, self.selected_node.attrs[key][1])
      self.node_attrs.get_column(0).queue_resize()
      self.node_attrs.get_column(1).queue_resize()
      self.node_attrs.get_column(2).queue_resize()

    return

  def node_attrs_mouse_over(self, widget, event):
    """
    Called when the mouse moves over the node attributes widget. Sets the
    appropriate attribute widget tooltip.
    """

    path_info = self.node_attrs.get_path_at_pos(int(event.x), int(event.y))
    if path_info is None:
      try:
        self.node_attrs.set_tooltip_text("")
        self.node_attrs.set_property("has-tooltip", False)
      except:
        pass
      return

    path = path_info[0]
    col = path_info[1]
    if col is not self.node_attrs.get_column(1):
      try:
        self.node_attrs.set_tooltip_text("")
        self.node_attrs.set_property("has-tooltip", False)
      except:
        pass
      return

    iter = self.node_attrs.get_model().get_iter(path)
    iter_key = self.node_attrs.get_model().get_value(iter, 0)

    return

  def node_attrs_key_data_func(self, col, cell_renderer, model, iter):
    """
    Attribute name data function. Sets the cell renderer text colours.
    """

    iter_key = model.get_value(iter, 0)

    if not self.selected_node.active or self.selected_node.attrs[iter_key][0] is None or self.selected_node.attrs[iter_key][0] == "fixed":
      cell_renderer.set_property("foreground", "grey")
    elif self.selected_node.attrs[iter_key][1] is None:
      cell_renderer.set_property("foreground", "blue")
    else:
      cell_renderer.set_property("foreground", "black")

    return

  def node_attrs_entry_data_func(self, col, cell_renderer, model, iter):
    """
    Attribute text data function. Hides the renderer if a combo box is required,
    and sets colours and editability otherwise.
    """

    iter_key = model.get_value(iter, 0)

    if not self.selected_node.active or self.selected_node.attrs[iter_key][0] is None or self.selected_node.attrs[iter_key][0] == "fixed":
      cell_renderer.set_property("editable", False)
      cell_renderer.set_property("foreground", "grey")
      cell_renderer.set_property("visible", True)
    elif not isinstance(self.selected_node.attrs[iter_key][0], tuple):
      cell_renderer.set_property("editable", True)
      cell_renderer.set_property("visible", True)
      if self.selected_node.attrs[iter_key][1] is None:
        cell_renderer.set_property("foreground", "blue")
      else:
        cell_renderer.set_property("foreground", "black")
    else:
      cell_renderer.set_property("editable", False)
      cell_renderer.set_property("visible", False)

    return

  def node_attrs_combo_data_func(self, col, cell_renderer, model, iter):
    """
    Attribute combo box data function. Hides the renderer if a combo box is not
    required, and sets the combo box options otherwise. Adds an entry if required.
    """

    iter_key = model.get_value(iter, 0)

    if self.selected_node.active and isinstance(self.selected_node.attrs[iter_key][0], tuple):
      cell_renderer.set_property("editable", True)
      cell_renderer.set_property("visible", True)
      if isinstance(self.selected_node.attrs[iter_key][0][0], tuple):
        cell_renderer.set_property("has-entry", True)
      else:
        cell_renderer.set_property("has-entry", False)
      if self.selected_node.attrs[iter_key][1] is None:
        cell_renderer.set_property("foreground", "blue")
      else:
        cell_renderer.set_property("foreground", "black")
    else:
      cell_renderer.set_property("visible", False)
      cell_renderer.set_property("editable", False)
    cell_renderer.set_property("model", model.get_value(iter, 2))

    return

  def node_attrs_icon_data_func(self, col, cell_renderer, model, iter):
    """
    Attribute icon data function. Used to add downward pointing arrows for combo
    attributes, for consistency with the LHS.
    """

    iter_key = model.get_value(iter, 0)

    if self.selected_node.active and isinstance(self.selected_node.attrs[iter_key][0], tuple):
      cell_renderer.set_property("stock-id", gtk.STOCK_GO_DOWN)
    else:
      cell_renderer.set_property("stock-id", None)

    return

  def node_attrs_entry_edit_start(self, cell_renderer, editable, path):
    """
    Called when editing is started on an attribute text cell. Used to delete the
    printable_type placeholder.
    """

    iter = self.node_attrs.get_model().get_iter(path)
    iter_key = self.node_attrs.get_model().get_value(iter, 0)

    if self.selected_node.attrs[iter_key][1] is None:
      editable.set_text("")

    return

  def node_attrs_combo_edit_start(self, cell_renderer, editable, path):
    """
    Called when editing is started on an attribute combo cell. Used to delete the
    select placeholder for mixed entry / combo attributes.
    """

    iter = self.node_attrs.get_model().get_iter(path)
    iter_key = self.node_attrs.get_model().get_value(iter, 0)

    if isinstance(self.selected_node.attrs[iter_key][0][0], tuple) and self.selected_node.attrs[iter_key][1] is None:
      editable.child.set_text("")

    return

  def node_attrs_edited(self, cell_renderer, path, new_text):
    """
    Called when editing is finished on an attribute text cell. Updates data in the
    treestore.
    """

    iter = self.node_attrs.get_model().get_iter(path)
    iter_key = self.node_attrs.get_model().get_value(iter, 0)

    if self.selected_node.get_attr(iter_key) is None and new_text == "":
      return

    value_check = self.validity_check(new_text, self.selected_node.attrs[iter_key][0])

    if not value_check is None and not value_check == self.selected_node.attrs[iter_key][1]:
      if iter_key == "name" and not self.name_check(value_check):
        return

      self.node_attrs.get_model().set_value(iter, 1, value_check)
      self.selected_node.set_attr(iter_key, value_check)
      self.paint_validity()
      if iter_key == "name":
        self.update_painted_name()
        self.on_select_row()
      self.set_saved(False)

    return

  def node_attrs_selected(self, cell_renderer, path, new_text):
    """
    Called when an attribute combo box element is selected, or combo box entry
    element entry is edited. Updates data in the treestore.
    """

    iter = self.node_attrs.get_model().get_iter(path)
    iter_key = self.node_attrs.get_model().get_value(iter, 0)

    if new_text is None:
      return

    if isinstance(self.selected_node.attrs[iter_key][0][0], tuple) and not new_text in self.selected_node.attrs[iter_key][0][0]:
      if self.selected_node.get_attr(iter_key) is None and new_text == "":
        return
        
      new_text = self.validity_check(new_text, self.selected_node.attrs[iter_key][0][1])
      if iter_key == "name" and not self.name_check(new_text):
        return False
    if not new_text == self.selected_node.attrs[iter_key][1]:
      self.node_attrs.get_model().set_value(iter, 1, new_text)
      self.selected_node.set_attr(iter_key, new_text)
      self.paint_validity()
      if iter_key == "name":
        self.update_painted_name()
      self.set_saved(False)

    return

  def set_node_data_empty(self):
    """
    Empty the node data frame.
    """

    if len(self.node_data_frame.get_children()) > 1:
      if isinstance(self.node_data, gtk.TextView):
        self.node_data.handler_block_by_func(self.node_data_entry_focus_in)
      elif isinstance(self.node_data, gtk.ComboBox):
        self.node_data.handler_block_by_func(self.node_data_combo_focus_child)
      self.node_data_frame.remove(self.node_data_frame.child)

    self.node_data_interacted = False

    return

  def set_node_data_entry(self):
    """
    Create a text view for data entry in the node data frame.
    """

    self.set_node_data_empty()

    data_scrolled_window = gtk.ScrolledWindow()
    self.node_data_frame.add(data_scrolled_window)
    data_scrolled_window.show()

    try:
      import gtksourceview2
      buf = gtksourceview2.Buffer()
      lang_manager = gtksourceview2.LanguageManager()
      buf.set_highlight_matching_brackets(True)
      if self.node_data_is_python_code():
        python = lang_manager.get_language("python")
        buf.set_language(python)
        buf.set_highlight_syntax(True)
      self.node_data = gtksourceview2.View(buffer=buf)
      self.node_data.set_auto_indent(True)
      #self.node_data.set_highlight_current_line(True)
      self.node_data.set_insert_spaces_instead_of_tabs(True)
      self.node_data.set_tab_width(2)
      if self.node_data_is_python_code():
        self.node_data.set_show_line_numbers(True)
        font_desc = pango.FontDescription("monospace")
        if font_desc:
          self.node_data.modify_font(font_desc)
    except ImportError:
      self.node_data = gtk.TextView()

    data_scrolled_window.add(self.node_data)
    self.node_data.show()

    data_scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

    self.node_data.set_pixels_above_lines(2)
    self.node_data.set_pixels_below_lines(2)
    self.node_data.set_wrap_mode(gtk.WRAP_WORD)

    self.node_data.connect("focus-in-event", self.node_data_entry_focus_in)

    data_frame_packing = self.node_data_frame.get_property("parent").query_child_packing(self.node_data_frame)
    self.node_data_frame.get_property("parent").set_child_packing(self.node_data_frame, True, data_frame_packing[1], data_frame_packing[2], data_frame_packing[3])

    self.node_data.get_buffer().create_tag("node_data_buffer_tag")
    text_tag = self.node_data.get_buffer().get_tag_table().lookup("node_data_buffer_tag")
    if self.selected_node is None:
      self.node_data.set_cursor_visible(False)
      self.node_data.set_editable(False)
      self.node_data_buttons_hbox.hide()
      self.node_data.get_buffer().set_text("")
      text_tag.set_property("foreground", "grey")
    elif not self.selected_node.active:
      self.node_data.set_cursor_visible(False)
      self.node_data.set_editable(False)
      self.node_data_buttons_hbox.hide()
      self.node_data.get_buffer().set_text("Inactive node")
      text_tag.set_property("foreground", "grey")
    elif self.selected_node.datatype is None:
      self.node_data.set_cursor_visible(False)
      self.node_data.set_editable(False)
      self.node_data_buttons_hbox.hide()
      self.node_data.get_buffer().set_text("No data")
      text_tag.set_property("foreground", "grey")
    elif self.node_data_is_tensor():
      self.node_data.set_cursor_visible(False)
      self.node_data.set_editable(False)
      self.node_data_buttons_hbox.hide()
      self.node_data.get_buffer().set_text("Dimension not set")
      text_tag.set_property("foreground", "grey")
    elif self.selected_node.data is None:
      self.node_data.set_cursor_visible(True)
      self.node_data.set_editable(True)
      self.node_data_buttons_hbox.show()
      self.node_data.get_buffer().set_text(self.printable_type(self.selected_node.datatype))
      text_tag.set_property("foreground", "blue")
    else:
      self.node_data.get_buffer().set_text(self.selected_node.data)
      if self.selected_node.datatype == "fixed":
        self.node_data.set_cursor_visible(False)
        self.node_data.set_editable(False)
        self.node_data_buttons_hbox.hide()
        text_tag.set_property("foreground", "grey")
      else:
        self.node_data.set_cursor_visible(True)
        self.node_data.set_editable(True)
        self.node_data_buttons_hbox.show()
        #text_tag.set_property("foreground", "black")
    buffer_bounds = self.node_data.get_buffer().get_bounds()
    self.node_data.get_buffer().apply_tag(text_tag, buffer_bounds[0], buffer_bounds[1])

    return

  def set_node_data_tensor(self):
    """
    Create a table container packed with appropriate widgets for tensor data entry
    in the node data frame.
    """

    self.set_node_data_empty()

    data_scrolled_window = gtk.ScrolledWindow()
    self.node_data_frame.add(data_scrolled_window)
    data_scrolled_window.show()

    data_scrolled_window.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)

    dim1, dim2 = self.node_data_tensor_shape()
    self.node_data = gtk.Table(dim1, dim2)
    data_scrolled_window.add_with_viewport(self.node_data)
    self.node_data.show()

    data_scrolled_window.child.set_property("shadow-type", gtk.SHADOW_NONE)

    data_frame_packing = self.node_data_frame.get_property("parent").query_child_packing(self.node_data_frame)
    self.node_data_frame.get_property("parent").set_child_packing(self.node_data_frame, True, data_frame_packing[1], data_frame_packing[2], data_frame_packing[3])

    self.node_data_buttons_hbox.show()

    is_symmetric = self.node_data_is_symmetric_tensor()
    for i in range(dim1):
      for j in range(dim2):
        entry = gtk.Entry()
        self.node_data.attach(entry, dim2 - j - 1, dim2 - j, dim1 - i - 1, dim1 - i)
        if not is_symmetric or i >= j:
          entry.show()

          entry.connect("focus-in-event", self.node_data_tensor_element_focus_in, dim2 - j - 1, dim1 - i - 1)

          if self.selected_node.data is None:
            entry.set_text(self.printable_type(self.selected_node.datatype.datatype))
            entry.modify_text(gtk.STATE_NORMAL, gtk.gdk.color_parse("blue"))
          else:
            entry.set_text(self.selected_node.data.split(" ")[(dim2 - j - 1) + (dim1 - i - 1) * dim2])

    self.node_data_interacted = [False for i in range(dim1 * dim2)]

    return

  def set_node_data_combo(self):
    """
    Create a combo box for node data selection in the node data frame. Add an
    entry if required.
    """

    self.set_node_data_empty()

    if isinstance(self.selected_node.datatype[0], tuple):
      self.node_data = gtk.combo_box_entry_new_text()
    else:
      self.node_data = gtk.combo_box_new_text()
    self.node_data_frame.add(self.node_data)
    self.node_data.show()

    self.node_data.connect("set-focus-child", self.node_data_combo_focus_child)
    self.node_data.connect("scroll-event", self.node_data_combo_scroll)

    data_frame_packing = self.node_data_frame.get_property("parent").query_child_packing(self.node_data_frame)
    self.node_data_frame.get_property("parent").set_child_packing(self.node_data_frame, False, data_frame_packing[1], data_frame_packing[2], data_frame_packing[3])

    if isinstance(self.selected_node.datatype[0], tuple):
      self.node_data_buttons_hbox.show()
    else:
      self.node_data_buttons_hbox.hide()

    if self.selected_node.data is None:
      if isinstance(self.selected_node.datatype[0], tuple):
        self.node_data.child.set_text("Select " + self.printable_type(self.selected_node.datatype[1]) + "...")
      else:
        self.node_data.append_text("Select...")
        self.node_data.set_active(0)
      self.node_data.child.modify_text(gtk.STATE_NORMAL, gtk.gdk.color_parse("blue"))
      self.node_data.child.modify_text(gtk.STATE_PRELIGHT, gtk.gdk.color_parse("blue"))

    if isinstance(self.selected_node.datatype[0], tuple):
      options = self.selected_node.datatype[0]
    else:
      options = self.selected_node.datatype
    for i in range(len(options)):
      opt = options[i]
      self.node_data.append_text(opt)
      if self.selected_node.data == opt:
        self.node_data.set_active(i)

    if isinstance(self.selected_node.datatype[0], tuple) and not self.selected_node.data is None and not self.selected_node.data in self.selected_node.datatype[0]:
      self.node_data.child.set_text(self.selected_node.data)

    self.node_data.connect("changed", self.node_data_combo_changed)

    return

  def node_data_entry_focus_in(self, widget, event):
    """
    Called when a text view data entry widget gains focus. Used to delete the
    printable_type placeholder.
    """

    if not self.selected_node is None and not self.selected_node.datatype is None and not self.node_data_is_tensor() and self.selected_node.data is None and not self.node_data_interacted:
      self.node_data.get_buffer().set_text("")

    self.node_data_interacted = True

    return

  def node_data_tensor_element_focus_in(self, widget, event, row, col):
    """
    Called when a tensor data entry widget gains focus. Used to delete the
    printable_type placeholder.
    """

    dim1, dim2 = self.node_data_tensor_shape()
    if not self.node_data_interacted[col + row * dim2]:
      self.node_data_interacted[col + row * dim2] = True
      if self.node_data_is_symmetric_tensor():
        self.node_data_interacted[row + col * dim1] = True
      if self.selected_node.data is None:
        widget.set_text("")
        widget.modify_text(gtk.STATE_NORMAL, gtk.gdk.color_parse("black"))

    return

  def node_data_combo_focus_child(self, container, widget):
    """
    Called when a data selection widget gains focus. Used to delete the select
    placeholder.
    """

    if not self.node_data_interacted:
      self.node_data_interacted = True
      if self.selected_node.data is None:
        self.node_data_interacted = True
        if isinstance(self.selected_node.datatype[0], tuple):
          self.node_data.handler_block_by_func(self.node_data_combo_changed)
          self.node_data.child.set_text("")
          self.node_data.handler_unblock_by_func(self.node_data_combo_changed)
        else:
          self.node_data.set_active(1)
          self.node_data.remove_text(0)
        self.node_data.child.modify_text(gtk.STATE_NORMAL, gtk.gdk.color_parse("black"))
        self.node_data.child.modify_text(gtk.STATE_PRELIGHT, gtk.gdk.color_parse("black"))

    return

  def node_data_combo_changed(self, combo_box):
    """
    Called when a data combo box element is selected. Updates data in the
    treestore.
    """

    if not isinstance(self.selected_node.datatype[0], tuple) or not self.node_data.child.get_property("has-focus"):
      self.selected_node.set_data(self.node_data.get_active_text())
      self.paint_validity()
      self.set_saved(False)
      self.node_data_interacted = False
      iter = self.get_treeview_iter(self.treeview.get_selection())
      self.treestore.set_value(iter, 4, self.node_data.get_active_text())

    return

  def node_data_combo_scroll(self, widget, event):
    """
    Called when the data combo box is scrolled with the mouse wheel. Removes the
    select placeholder and updates data in the treestore.
    """

    self.node_data_combo_focus_child(self.node_data_frame, self.node_data)
    if not isinstance(self.selected_node.datatype[0], tuple) or not self.selected_node.data is None:
      self.node_data_combo_changed(self.node_data)

    return

  def node_data_revert(self, button = None):
    """
    "Revert Data" button click signal handler. Reverts data in the node data frame.
    """

    if self.node_data_is_tensor() and not self.geometry_dim_tree.data is None:
      self.set_node_data_tensor()
    elif not self.selected_node is None and isinstance(self.selected_node.datatype, tuple):
      self.set_node_data_combo()
    else:
      self.set_node_data_entry()

    return

  def node_data_store(self, button = None):
    """
    "Store Data" button click signal handler. Stores data from the node data frame
    in the treestore.
    """

    if self.node_data_is_tensor() and not self.geometry_dim_tree.data is None:
      store_success = self.node_data_tensor_store()
    elif not self.selected_node is None and isinstance(self.selected_node.datatype, tuple):
      store_success = self.node_data_combo_store()
    else:
      store_success = self.node_data_entry_store()

    if self.scherror.errlist_is_open():
      if self.scherror.errlist_type == 0:
         self.scherror.on_validate_schematron()
      else:
         self.scherror.on_validate()

    return store_success

  def node_data_entry_store(self):
    """
    Attempt to store data read from a textview packed in the node data frame.
    """

    if self.selected_node is None or self.selected_node.datatype in ["fixed", None] or self.node_data_is_tensor():
      return True

    data_buffer_bounds = self.node_data.get_buffer().get_bounds()
    new_data = self.node_data.get_buffer().get_text(data_buffer_bounds[0], data_buffer_bounds[1])

    if new_data == "":
      return True

    if self.selected_node.data is None and not self.node_data_interacted:
      return True
    else:
      value_check = self.validity_check(new_data, self.selected_node.datatype)
      if value_check is None:
        return False
      elif not value_check == self.selected_node.data:
        self.selected_node.set_data(value_check)
        if isinstance(self.selected_node, MixedTree) and "shape" in self.selected_node.child.attrs.keys() and self.selected_node.child.attrs["shape"][0] is int and isinstance(self.selected_node.datatype, plist.List) and self.selected_node.datatype.cardinality == "+":
          self.selected_node.child.set_attr("shape", str(len(value_check.split(" "))))
        self.paint_validity()
        
        iter = self.selected_iter
        self.treestore.set_value(iter, 4, new_data)
        self.set_saved(False)
        self.node_data_interacted = False

    return True

  def node_data_combo_store(self):
    """
    Attempt to store data read from a combo box entry packed in the node data
    frame.
    """

    if not isinstance(self.selected_node.datatype[0], tuple):
      return True

    new_data = self.node_data.child.get_text()

    if self.selected_node.data is None and not self.node_data_interacted:
      return True
    elif not new_data in self.selected_node.datatype[0]:
      new_data = self.validity_check(new_data, self.selected_node.datatype[1])
      if new_data is None:
        return False

    if not new_data == self.selected_node.data:
      self.selected_node.set_data(new_data)
      self.paint_validity()
      self.set_saved(False)
      self.node_data_interacted = False

    return True

  def node_data_tensor_store(self):
    """
    Attempt to store data read from tensor data entry widgets packed in the node
    data frame.
    """

    dim1, dim2 = self.node_data_tensor_shape()
    is_symmetric = self.node_data_is_symmetric_tensor()

    if not True in self.node_data_interacted:
      return True

    entry_values = []
    for i in range(dim1):
      for j in range(dim2):
        if is_symmetric and i > j:
          entry_values.append(self.node_data.get_children()[i + j * dim1].get_text())
        else:
          entry_values.append(self.node_data.get_children()[j + i * dim2].get_text())

    changed = False
    for i in range(dim1):
      for j in range(dim2):
        if self.node_data_interacted[j + i * dim2] and not entry_values[j + i * dim2] == "" and (self.selected_node.data is None or not self.selected_node.data.split(" ")[j + i * dim2] == entry_values[j + i * dim2]):
          changed = True
    if not changed:
      return True
    elif (self.selected_node.data is None and False in self.node_data_interacted) or "" in entry_values:
      dialogs.error(self.main_window, "Invalid value entered")
      return False

    new_data = ""
    for i in range(dim1):
      for j in range(dim2):
        new_data += " " + entry_values[j + i * dim2]

    value_check = self.validity_check(new_data, self.selected_node.datatype)
    if value_check is None:
      return False
    elif not value_check == self.selected_node.data:
      self.selected_node.set_data(value_check)

      dim1, dim2 = self.node_data_tensor_shape()
      if int(self.selected_node.child.attrs["rank"][1]) == 1:
        self.selected_node.child.set_attr("shape", str(dim1))
      else:
        self.selected_node.child.set_attr("shape", str(dim1) + " " + str(dim2))

      self.paint_validity()
      self.set_saved(False)
      self.node_data_interacted = [False for i in range(dim1 * dim2)]

    return True

  def update_node_comment(self):
    """
    Update the comment widget.
    """

    if self.selected_node is None or not self.selected_node.active:
      self.node_comment.get_buffer().set_text("")
      self.node_comment.set_cursor_visible(False)
      self.node_comment.set_editable(False)
      try:
        self.node_comment.set_tooltip_text("")
        self.node_comment.set_property("has-tooltip", False)
      except:
        pass

      return

    comment_tree = self.get_comment(self.selected_node)
    text_tag = self.node_comment.get_buffer().get_tag_table().lookup("comment_buffer_tag")
    if comment_tree is None:
      self.node_comment.get_buffer().set_text("No comment")
      self.node_comment.set_cursor_visible(False)
      self.node_comment.set_editable(False)
      text_tag.set_property("foreground", "grey")
      try:
        self.node_comment.set_tooltip_text("")
        self.node_comment.set_property("has-tooltip", False)
      except:
        pass
    else:
      if comment_tree.data is None:
        self.node_comment.get_buffer().set_text("(string)")
      else:
        self.node_comment.get_buffer().set_text(comment_tree.data)
      if self.selected_node.active:
        self.node_comment.set_cursor_visible(True)
        self.node_comment.set_editable(True)
        text_tag.set_property("foreground", "black")
      else:
        self.node_comment.set_cursor_visible(False)
        self.node_comment.set_editable(False)
        text_tag.set_property("foreground", "grey")

    buffer_bounds = self.node_comment.get_buffer().get_bounds()
    self.node_comment.get_buffer().apply_tag(text_tag, buffer_bounds[0], buffer_bounds[1])

    self.node_comment_interacted = False

    return

  def node_comment_focus_in(self, widget, event):
    """
    Called when the comment widget gains focus. Removes the printable_type
    placeholder.
    """

    comment_tree = self.get_comment(self.selected_node)
    if not comment_tree is None and not self.node_comment_interacted:
      self.node_comment_interacted = True
      if comment_tree.data is None:
        self.node_comment.get_buffer().set_text("")

    return

  def node_comment_expose(self, widget, event):
    """
    Called when the comment widget is repainted. Stores the comment if required.
    """

    self.node_comment_store()

    return

  def node_comment_store(self):
    """
    Store data in the node comment.
    """

    comment_tree = self.get_comment(self.selected_node)
    if comment_tree is None or not self.node_comment_interacted:
      return

    data_buffer_bounds = self.node_comment.get_buffer().get_bounds()
    new_comment = self.node_comment.get_buffer().get_text(data_buffer_bounds[0], data_buffer_bounds[1])

    if not new_comment == comment_tree.data:
      if new_comment == "":
        comment_tree.data = None
        comment_tree.active = False
      else:
        comment_tree.set_data(new_comment)
        comment_tree.active = True
        self.set_saved(False)

    return

  def validity_check(self, val, val_type):
    """
    Check to see if the supplied data with supplied type can be stored in a
    tree.Tree.
    """

    (invalid, data) = self.selected_node.valid_data(val_type, val)
    if not invalid and isinstance(data, str) and not data == "":
      if not data == val and self.validity_check(data, val_type) is None:
        return None
      else:
        return data
    else:
      dialogs.error(self.main_window, "Invalid value entered")
      return None

  def name_check(self, val):
    """
    Check to see if the supplied data is a valid tree name.
    """

    valid_chars = "_:[]1234567890qwertyuioplkjhgfdsazxcvbnmMNBVCXZASDFGHJKLPOIUYTREWQ"
    for char in val:
      if not char in valid_chars:
        dialogs.error(self.main_window, "Invalid value entered")
        return False

    return True

  def node_data_is_python_code(self):
    """
    Perform a series of tests on the current tree.Tree / MixedTree, to determine if
    it is intended to be used to store python code data.
    """

    try:
       lang = self.selected_node.get_attr("language")
       if lang == "python":
         return True
    except:
      pass
    
    if not isinstance(self.selected_node, MixedTree):
      return False
      
    if not self.selected_node.datatype is str:
      return False
  
    if "type" in self.selected_node.child.attrs.keys():
      return self.selected_node.child.attrs["type"][1] == "python"
    else:
      return False

  def node_data_is_tensor(self):
    """
    Perform a series of tests on the current tree.Tree / MixedTree, to determine if
    it is intended to be used to store tensor or vector data.
    """

    # Check that a geometry is defined
    if self.geometry_dim_tree is None:
      return False

    # Check that this element has calculable and positive dimensions
    if isinstance(self.geometry_dim_tree.datatype, tuple):
      possible_dims = self.geometry_dim_tree.datatype
    else:
      possible_dims = [self.geometry_dim_tree.data]
    for opt in possible_dims:
      try:
        dim1, dim2 = self.node_data_tensor_shape(int(opt))
        assert dim1 > 0
        assert dim2 > 0
      except:
        return False

    # All tensor elements must be of MixedTree type
    if not isinstance(self.selected_node, MixedTree):
      return False

    # The element must have dim1, rank and shape attributes
    if not "dim1" in self.selected_node.child.attrs.keys() or not "rank" in self.selected_node.child.attrs.keys() or not "shape" in self.selected_node.child.attrs.keys():
      return False
    # The dim1 and rank attributes must be of fixed type
    if not self.selected_node.child.attrs["dim1"][0] == "fixed" or not self.selected_node.child.attrs["rank"][0] == "fixed":
      return False

    if "dim2" in self.selected_node.child.attrs.keys():
      # If a dim2 attribute is specified, it must be of fixed type and the rank must be 2
      # Also, the shape attribute must be a list of integers with cardinality equal to the rank
      if not self.selected_node.child.attrs["dim2"][0] == "fixed" or not self.selected_node.child.attrs["rank"][1] == "2" or not isinstance(self.selected_node.child.attrs["shape"][0], plist.List) or not self.selected_node.child.attrs["shape"][0].datatype is int or not str(self.selected_node.child.attrs["shape"][0].cardinality) == self.selected_node.child.attrs["rank"][1]:
        return False
    # Otherwise, the rank must be one and the shape an integer
    elif not self.selected_node.child.attrs["rank"][1] == "1" or not self.selected_node.child.attrs["shape"][0] is int: 
      return False

    # The data for the element must be a list of one or more
    if not isinstance(self.selected_node.datatype, plist.List) or not self.selected_node.datatype.cardinality == "+":
      return False

    # If the shape has been set, check that it has a valid value
    if not self.selected_node.child.attrs["shape"][1] == None:
      if self.geometry_dim_tree.data is None:
        return False

      dim1, dim2 = self.node_data_tensor_shape()
      if "dim2" in self.selected_node.child.attrs.keys():
        if not self.selected_node.child.attrs["shape"][1] == str(dim1) + " " + str(dim2):
          return False
      elif not self.selected_node.child.attrs["shape"][1] == str(dim1):
        return False

    return True

  def node_data_tensor_shape(self, geometry_dim = None):
    """
    Read the tensor shape for tensor or vector data in the current MixedTree.
    """

    if geometry_dim == None:
      geometry_dim = int(self.geometry_dim_tree.data)

    dim1 = 1
    dim2 = 1
    if "dim1" in self.selected_node.child.attrs.keys():
      dim1 = int(eval(self.selected_node.child.attrs["dim1"][1], {"dim":geometry_dim}))
      if "dim2" in self.selected_node.child.attrs.keys():
        dim2 = int(eval(self.selected_node.child.attrs["dim2"][1], {"dim":geometry_dim}))

    return (dim1, dim2)

  def node_data_is_symmetric_tensor(self):
    """
    Read if the tensor data in the current MixedTree is symmetric.
    """

    dim1, dim2 = self.node_data_tensor_shape()
    if not dim1 == dim2:
      return False

    if "symmetric" in self.selected_node.child.attrs.keys() and self.selected_node.child.attrs["symmetric"][1] == "true":
      return True
    else:
      return False

class MixedTree:
  def __init__(self, parent, child):
    """
    The .doc and .attrs comes from parent, and the .data comes from child. This
    is used to hide integer_value etc. from the left hand side, but for its data
    entry to show up on the right.
    """

    self.parent = parent
    self.child = child

    self.name = parent.name
    self.schemaname = parent.schemaname
    self.attrs = self.parent.attrs
    self.children = parent.children
    self.datatype = child.datatype
    self.data = child.data
    self.doc = parent.doc
    self.active = parent.active

    return

  def set_attr(self, attr, val):
    self.parent.set_attr(attr, val)

    return

  def get_attr(self, attr):
    return self.parent.get_attr(attr)

  def set_data(self, data):
    self.child.set_data(data)
    self.datatype = self.child.datatype
    self.data = self.child.data

    return

  def valid_data(self, datatype, data):
    return self.parent.valid_data(datatype, data)

  def matches(self, text, case_sensitive = False):
    old_parent_data = self.parent.data
    self.parent.data = None
    parent_matches = self.parent.matches(text, case_sensitive)
    self.parent.data = old_parent_data

    if parent_matches:
      return True

    if case_sensitive:
      text_re = re.compile(text)
    else:
      text_re = re.compile(text, re.IGNORECASE)

    if not self.child.data is None and not text_re.search(self.child.data) is None:
      return True
    else:
      return False

class DiamondFindDialog:
  def __init__(self, parent, gladefile):
    self.parent = parent
    self.gladefile = gladefile
    self.search_dialog = None

    return

  def on_find(self, widget=None):
    """
    Open up the find dialog. It has to be created each time from the glade file.
    """

    if not self.search_dialog is None:
      return

    signals =      {"on_find_dialog_close": self.on_find_close_button,
                    "on_close_clicked": self.on_find_close_button,
                    "on_find_clicked": self.on_find_find_button}

    self.search_gui = gtk.glade.XML(self.gladefile, root="find_dialog")
    self.search_dialog = self.search_gui.get_widget("find_dialog")
    self.search_gui.signal_autoconnect(signals)
    search_entry = self.search_gui.get_widget("search_entry")
    search_entry.connect("activate", self.on_find_find_button)

    # reset the search parameters
    self.search_generator = None
    self.search_text = ""
    self.search_count = 0
    self.search_dialog.show()
    self.parent.statusbar.set_statusbar("")

    return

  def on_find_find_button(self, button):
    """
    Search. Each time "Find" is clicked, we compare the stored search text to the
    text in the entry box. If it's the same, we find next; if it's different, we
    start a new search. self.search_treestore does the heavy lifting.
    """

    search_entry = self.search_gui.get_widget("search_entry")

    self.parent.statusbar.clear_statusbar()

    text = search_entry.get_text()
    if text == "":
      self.parent.statusbar.set_statusbar("No text")
      return

    # check if we've started a new search
    if text != self.search_text:
      # started a new search
      self.search_generator = None
      self.search_generator = self.parent.search_treestore(text)
      self.search_text = text
      self.search_count = 0

    try:
      # get the iter of the next tree that matches
      iter = self.search_generator.next()
      path = self.parent.treestore.get_path(iter)
      # scroll down to it, expand it, and select it
      self.parent.treeview.expand_to_path(path)
      self.parent.treeview.get_selection().select_iter(iter)
      self.parent.treeview.scroll_to_cell(path, use_align=True, col_align=0.5)
      # count how many hits we've had
      self.search_count = self.search_count + 1
    except StopIteration:
      # reset the search and cycle
      self.search_text = ""
      # if something was found, go through again
      if self.search_count > 0:
        self.on_find_find_button(button)
      else:
        self.parent.statusbar.set_statusbar("No results")

    return

  def on_find_close_button(self, button = None):
    """
    Close the search widget.
    """

    if not self.search_dialog is None:
      self.search_dialog.hide()
      self.search_dialog = None
    self.parent.statusbar.clear_statusbar()

    return

class DiamondStatusBar:
  def __init__(self, statusbar):
    self.statusbar = statusbar
    self.context_id = statusbar.get_context_id("Messages")

    return

  def set_statusbar(self, msg):
    """
    Set the status bar message.
    """

    self.statusbar.push(self.context_id, msg)

    return

  def clear_statusbar(self):
    """
    Clear the status bar.
    """

    self.statusbar.push(self.context_id, "")

    return
