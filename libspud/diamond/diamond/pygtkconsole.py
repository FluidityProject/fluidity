#  GTK Interactive Console
#  (C) 2003, Jon Anderson
#  (C) 2007 Imperial College London and others.
# Please see the AUTHORS file for a full list of contributors.

#  See www.python.org/2.2/license.htm for
#  license details.
#

import code
import sys

from gi.repository import Gtk as gtk
from gi.repository import Gdk as gdk
from gi.repository import Pango as pango

import builtins
import __main__

from . import debug

banner = """GTK Interactive Python Console
%s
""" % sys.version

class Completer:
  """
  Taken from rlcompleter, with readline references stripped, and a local dictionary to use.
  """
  def __init__(self, locals):
    self.locals = locals

  def complete(self, text, state):
    """Return the next possible completion for 'text'.
    This is called successively with state == 0, 1, 2, ... until it
    returns None.  The completion should begin with 'text'.

    """
    if state == 0:
      if "." in text:
        self.matches = self.attr_matches(text)
      else:
        self.matches = self.global_matches(text)
    try:
      return self.matches[state]
    except IndexError:
      return None

  def global_matches(self, text):
    """Compute matches when text is a simple name.

    Return a list of all keywords, built-in functions and names
    currently defines in __main__ that match.

    """
    import keyword
    matches = []
    n = len(text)
    for list in [keyword.kwlist,list(builtins.__dict__.keys()),list(__main__.__dict__.keys()), list(self.locals.keys())]:
      for word in list:
        if word[:n] == text and word != "__builtins__":
          matches.append(word)
    return matches

  def attr_matches(self, text):
    """Compute matches when text contains a dot.

    Assuming the text is of the form NAME.NAME....[NAME], and is
    evaluatable in the globals of __main__, it will be evaluated
    and its attributes (as revealed by dir()) are used as possible
    completions.  (For class instances, class members are are also
    considered.)

    WARNING: this can still invoke arbitrary C code, if an object
    with a __getattr__ hook is evaluated.

    """
    import re
    m = re.match(r"(\w+(\.\w+)*)\.(\w*)", text)
    if not m:
      return
    expr, attr = m.group(1, 3)
    object = eval(expr, __main__.__dict__, self.locals)
    words = dir(object)
    if hasattr(object,'__class__'):
      words.append('__class__')
      words = words + get_class_members(object.__class__)
    matches = []
    n = len(attr)
    for word in words:
      if word[:n] == attr and word != "__builtins__":
        matches.append("%s.%s" % (expr, word))
    return matches

def get_class_members(klass):
  ret = dir(klass)
  if hasattr(klass,'__bases__'):
     for base in klass.__bases__:
       ret = ret + get_class_members(base)
  return ret



class OutputStream:
  """
  A Multiplexing output stream.
  It can replace another stream, and tee output to the original stream and too
  a GTK textview.
  """
  def __init__(self,view,old_out,style):
    self.view = view
    self.buffer = view.get_buffer()
    self.mark = self.buffer.create_mark("End",self.buffer.get_end_iter(), False )
    self.out = old_out
    self.style = style
    self.tee = 1

  def write(self,text):
    if self.tee:
      debug.dwrite(self.out, text, newline = False, flush = False)

    end = self.buffer.get_end_iter()

    if not self.view  == None:
      self.view.scroll_to_mark(self.mark, 0, True, 1, 1)

    self.buffer.insert_with_tags(end,text,self.style)
    
  def flush(self):
    self.view.queue_draw()
  
    return

class GTKInterpreterConsole(gtk.ScrolledWindow):
  """
  An InteractiveConsole for GTK. It's an actual widget,
  so it can be dropped in just about anywhere.
  """
  def __init__(self, locals = None):
    gtk.ScrolledWindow.__init__(self)
    self.set_hexpand(True)
    self.set_vexpand(True)
    self.set_policy (gtk.PolicyType.AUTOMATIC,gtk.PolicyType.AUTOMATIC)

    self.text = gtk.TextView()
    self.text.set_wrap_mode(True)

    self.interpreter = code.InteractiveInterpreter(locals)

    self.completer = Completer(self.interpreter.locals)
    self.buffer = []
    self.history = []
    self.banner = banner
    self.ps1 = ">>> "
    self.ps2 = "... "

    self.text.add_events( gdk.EventMask.KEY_PRESS_MASK )
    self.text.connect( "key_press_event", self.key_pressed )

    self.current_history = -1

    self.mark = self.text.get_buffer().create_mark("End",self.text.get_buffer().get_end_iter(), False )

            #setup colors
    self.style_banner = gtk.TextTag.new("banner")
    self.style_banner.set_property( "foreground", "saddle brown" )

    self.style_ps1 = gtk.TextTag.new("ps1")
    self.style_ps1.set_property( "foreground", "DarkOrchid4" )
    self.style_ps1.set_property( "editable", False )
    self.style_ps1.set_property("font", "courier" )

    self.style_ps2 = gtk.TextTag.new("ps2")
    self.style_ps2.set_property( "foreground", "DarkOliveGreen" )
    self.style_ps2.set_property( "editable", False  )
    self.style_ps2.set_property("font", "courier" )

    self.style_out = gtk.TextTag.new("stdout")
    self.style_out.set_property( "foreground", "midnight blue" )
    self.style_err = gtk.TextTag.new("stderr")
    self.style_err.set_property( "style", pango.Style.ITALIC )
    self.style_err.set_property( "foreground", "red" )

    self.text.get_buffer().get_tag_table().add(self.style_banner)
    self.text.get_buffer().get_tag_table().add(self.style_ps1)
    self.text.get_buffer().get_tag_table().add(self.style_ps2)
    self.text.get_buffer().get_tag_table().add(self.style_out)
    self.text.get_buffer().get_tag_table().add(self.style_err)

    self.stdout = OutputStream(self.text,sys.stdout,self.style_out)
    self.stderr = OutputStream(self.text,sys.stderr,self.style_err)

    sys.stderr = self.stderr
    sys.stdout = self.stdout

    self.current_prompt = None

    self.write_line(self.banner, self.style_banner)
    self.prompt_ps1()

    self.add(self.text)
    self.text.show()



  def reset_history(self):
    self.history = []

  def reset_buffer(self):
    self.buffer = []

  def prompt_ps1(self):
    self.current_prompt = self.prompt_ps1
    self.write_line(self.ps1,self.style_ps1)

  def prompt_ps2(self):
    self.current_prompt = self.prompt_ps2
    self.write_line(self.ps2,self.style_ps2)

  def write_line(self,text,style=None):
    start,end = self.text.get_buffer().get_bounds()
    if style==None:
      self.text.get_buffer().insert(end,text)
    else:
      self.text.get_buffer().insert_with_tags(end,text,style)

    self.text.scroll_to_mark(self.mark, 0, True, 1, 1)

  def push(self, line):

    self.buffer.append(line)
    if len(line) > 0:
      self.history.append(line)

    source = "\n".join(self.buffer)

    more = self.interpreter.runsource(source, "<<console>>")

    if not more:
      self.reset_buffer()

    return more

  def key_pressed(self,widget,event):
    if event.keyval == gdk.keyval_from_name('Return'):
      return self.execute_line()

    if event.keyval == gdk.keyval_from_name('Up'):
      self.current_history = self.current_history - 1
      if self.current_history < - len(self.history):
        self.current_history = - len(self.history)
      return self.show_history()
    elif event.keyval == gdk.keyval_from_name('Down'):
      self.current_history = self.current_history + 1
      if self.current_history > 0:
        self.current_history = 0
      return self.show_history()
    elif event.keyval == gdk.keyval_from_name( 'Home'):
      l = self.text.get_buffer().get_line_count() - 1
      start = self.text.get_buffer().get_iter_at_line_offset(l,4)
      self.text.get_buffer().place_cursor(start)
      return True
    elif event.keyval == gdk.keyval_from_name( 'space') and event.state & gdk.ModifierType.CONTROL_MASK:
      return self.complete_line()
    return False

  def show_history(self):
    if self.current_history == 0:
      return True
    else:
      self.replace_line( self.history[self.current_history] )
      return True

  def current_line(self):
    start,end = self.current_line_bounds()
    return self.text.get_buffer().get_text(start,end, True)

  def current_line_bounds(self):
    txt_buffer = self.text.get_buffer()
    l = txt_buffer.get_line_count() - 1

    start = txt_buffer.get_iter_at_line(l)
    if start.get_chars_in_line() >= 4:
      start.forward_chars(4)
    end = txt_buffer.get_end_iter()
    return start,end

  def replace_line(self,txt):
    start,end = self.current_line_bounds()
    self.text.get_buffer().delete(start,end)
    self.write_line(txt)

  def execute_line(self):
    line = self.current_line()

    self.write_line("\n")

    more = self.push(line)

    self.text.get_buffer().place_cursor(self.text.get_buffer().get_end_iter())

    if more:
        self.prompt_ps2()
    else:
        self.prompt_ps1()


    self.current_history = 0

    return True

  def complete_line(self):
    line = self.current_line()
    tokens = line.split()
    token = tokens[-1]

    completions = []
    p = self.completer.complete(token,len(completions))
    while p != None:
      completions.append(p)
      p = self.completer.complete(token, len(completions))

    if len(completions) != 1:
      self.write_line("\n")
      self.write_line("\n".join(completions), self.style_ps1)
      self.write_line("\n")
      self.current_prompt()
      self.write_line(line)
    else:
      i = line.rfind(token)
      line = line[0:i] + completions[0]
      self.replace_line(line)

    return True
