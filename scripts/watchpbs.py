#!/usr/bin/env python

import gtk
import gobject
import dbus
import sys
import copy
import os
import optparse

o = old_state = {}

class Job(object):
  def __init__(self, line):
    self.id = line[0].split(".")[0]
    self.user = line[1]
    self.queue = line[2]
    self.jobname = line[3]
    self.maxtime = line[8]
    self.state = line[9]
    self.runtime = line[10]

def watchpbs(widget, filter, host):
  global o
  try:
    (c,r) = fetch_pbs_info(host, filter)
    set_tooltip(widget, c)

    message = process_keys(c, o)

    if message:
      notify_message(widget, message)

    o = c
    widget.set_visible(True)
  except:
    print "Polling failed."
    widget.set_visible(False)
  return True

def set_tooltip(widget, c):
  tooltip = ""
  ck = copy.copy(c.keys()); ck.sort()
  for key in ck:
    j = c[key]
    tooltip += "%s %s %s %-15s\t %s %s %s\n" % (j.id, j.user, j.queue, j.jobname, j.maxtime, j.state,
      j.runtime)

  widget.set_tooltip(tooltip)

def notify_message(widget, message):
  r = widget.get_geometry()[1]
  session_bus = dbus.SessionBus()

  # Use notifications object
  notifications_object = session_bus.get_object('org.freedesktop.Notifications',
      '/org/freedesktop/Notifications')
  n = dbus.Interface(notifications_object,
      'org.freedesktop.Notifications')

  # Sample empty notification
  notification_id = n.Notify(sys.argv[0], 0, '', 'PBS status',
      message, dbus.Array([], signature='s'), {'x': r.x - 8, 'y': r.y + 15}, -1)

def process_keys(c, o):

  ok = copy.copy(o.keys()); ok.sort()
  ck = copy.copy(c.keys()); ck.sort()
  message = ""

  # Case 0. Job finished.
  for key in ok:
    if key not in ck:
      message += "Job terminated: (%s, %s)\n" % (o[key].id, o[key].jobname)

 # Case 1. Change of state.
  for key in ok:
    if key in ck:
      if o[key].state != c[key].state:
        message += "Job state change: (%s, %s): %s -> %s\n" % (o[key].id,
           o[key].jobname, o[key].state, c[key].state)

  # Case 2. New job.
  for key in ck:
    if key not in ok:
      if c[key].state == 'R':
        message += "Job running: (%s, %s)\n" % (c[key].id, c[key].jobname)
      else:
        message += "Job queued: (%s, %s)\n" % (c[key].id, c[key].jobname)

  return message

def fetch_pbs_info(host, filter):
  if filter:
    raw = os.popen("ssh %s '/usr/pbs/bin/qstat -an | grep ^[0-9] | grep %s'" % (host, filter))
  else:
    raw = os.popen("ssh %s '/usr/pbs/bin/qstat -an | grep ^[0-9]'" % host)

  c = {}
  r = raw.read()
  if not r: raise Exeption
  print r
  for line in r.split('\n'):
    if not line: continue
    job = Job(line.split())
    c[job.id] = job

  return (c, r)

def quit_cb(widget, data = None):
  if data:
    data.set_visible(False)
  gtk.main_quit()

def popup_menu_cb(widget, button, time, data = None):
  if button == 3:
    if data:
      data.show_all()
      data.popup(None, None, None, 3, time)
  pass

def initialise_notify():
   s = gtk.StatusIcon()
   s.set_from_stock(gtk.STOCK_INFO)
   s.set_tooltip("No polling done yet ..")
   s.set_visible(True)

   return s


if __name__ == "__main__":
  usage = "usage: %prog [--host=HOST] [--filter=FILTER] [--time=SECONDS]"
  parser = optparse.OptionParser(usage)
  parser.add_option("--host", dest="host", default="cx1",
                    help="Server to ssh to")
  parser.add_option("--filter", dest="filter", default="",
                    help="Regular expression to filter qstat output")
  parser.add_option("--time", dest="time", default=60, type="int",
                    help="Number of seconds to wait between polls")

  (options, args) = parser.parse_args()

  s = initialise_notify()
  m = gtk.Menu()
  menuItem = gtk.ImageMenuItem(gtk.STOCK_QUIT)
  menuItem.connect('activate', quit_cb, s)
  m.append(menuItem)
  s.connect('popup-menu', popup_menu_cb, m)
  watchpbs(s, options.filter, options.host)
  gobject.timeout_add(options.time * 1000, watchpbs, s, options.filter, options.host)
  gtk.main()
