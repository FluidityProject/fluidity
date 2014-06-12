#!/usr/bin/env python


indent = '   '
verbose_state = False

import sys
import os

def default_fluidity_path():
    """Ensures the location of the Fluidity tree is known via
    FLUIDITYPATH, a system variable.  FLUIDITYPATH is either set by the
    user beforehand or defaulted here based on the present path.
    Updates the Python search paths accordingly."""
    try:
        flpath = os.environ['FLUIDITYPATH']
    except KeyError:
        flpath = "../../"
        os.environ['FLUIDITYPATH'] = flpath
    pypath = flpath + "/python/"
    os.environ['PYTHONPATH'] = pypath
    sys.path.append(pypath)
    
# setter/getter of verbose_state
def verbose(new_verbose_state=None):
    global verbose_state
    if new_verbose_state is not None:
        verbose_state = new_verbose_state
    else:
        return verbose_state

    
class Command:
    """Class to encapsulate a family of commands corresponding to a tree
    of Handlers (see below).  The commands are coded up for different
    levels of the tree in the 'execute' method.
    """
    
    def execute(self, level_name, value, indent):
        """The level_name argument tells the Command where it is in the
        Handler tree.  The Command can then choose what to do
        accordingly.  It may use additional information supplied in the
        value argument.
        """
        pass
        

class CommandList(Command):
    """Lumps multiple commands together."""
    
    def __init__(self, commands):
        self.commands = commands
        
    def execute(self, level_name, value, indent):
        for cmd in self.commands:
            cmd.execute(level_name, value, indent)


class Handler:
    """A composite to enable recursion over arbitrary combinations of
    simulation options (e.g. domain dimensions, grid resolutions).  The
    'handle' method combines the recursion with executing of a Command.
    """
    
    def __init__(self, level_name, value):
        self.level_name = level_name
        self.value = value
        self.level_index = 0
        self.beg_char = ''
        self.end_char = '\n'
        
    def set_level_index(self, level_index):
        pass

    def handle(self, command):
        pass

        
class CompositeHandler(Handler):
    """A tree of handlers."""

    def __init__(self, level_name, value, children):
        Handler.__init__(self, level_name, value)
        try:
            # if children is a HandlerLevel (see below), convert it back
            self.children = children.get_handlers()
        except AttributeError:
            # already a Python list of Handlers
            self.children = children
        for child in self.children:
            child.set_level_index(self.level_index + 1)
        
    def set_level_index(self, level_index):
        self.level_index = level_index
        if level_index > 0:
            self.beg_char = '\n'
            self.end_char = ''
        # recurse
        for child in self.children:
            child.set_level_index(self.level_index + 1)
        
    def handle(self, command):
        if verbose(): sys.stdout.write(self.beg_char+self.level_index*indent + \
           self.level_name + ': ' + str(self.value))
        # try executing the command at this level
        command.execute(self.level_name, self.value,
                        self.level_index*indent)
        # recurse
        for child in self.children:
            child.handle(command)
        if verbose(): sys.stdout.write(self.end_char)
        
        
class LeafHandler(Handler):
    """A single handler."""

    def __init__(self, level_name, value):
        Handler.__init__(self, level_name, value)
        
    def set_level_index(self, level_index):
        self.level_index = level_index
        if level_index > 0:
            self.beg_char = '\n'
            self.end_char = ''

    def handle(self, command):
        if verbose(): sys.stdout.write(self.beg_char+self.level_index*indent+\
            self.level_name+': '+str(self.value))
        # try executing the command at this level
        command.execute(self.level_name, self.value,
                        self.level_index*indent)
        if verbose(): sys.stdout.write(self.end_char)

        
class HandlerLevel:
    """A Python list of Handlers.  It is different to a CompositeHandler in
    that there is no parent-child relationship -- here there is simply
    some children.  It has extra methods for building trees and doesn't
    need to be polymorphic.

    """

    def __init__(self, level_name, values):
        """Using LeafHandler, initialises from a Python list of generic values.
        """        
        self.handlers = []
        for value in values:
            self.handlers.append( LeafHandler(level_name, value) )

    def add_sub(self, handlers_or_handler_level):
        """Using CompositeHandler, takes another HandlerLevel, or list of
        Handlers, and appends copies to each element in its own list of
        Handlers, effectively adding a sub-level to each Handler.
        I.e. does self.handlers[i].append(handlers).  Returns a fresh
        HandlerLevel; does not self-modify.

        """
        new_handler_level = HandlerLevel(None, [])
        new_handlers = []
        try:
            sub_handlers = handlers_or_handler_level.handlers
        except AttributeError:
            sub_handlers = handlers_or_handler_level
        for handler in self.handlers:
            new_handlers.append(
                CompositeHandler(handler.level_name, handler.value, 
                                 sub_handlers))
        new_handler_level.handlers = new_handlers
        return new_handler_level
        
    def splice(self, handler_levels):
        # DEPRECATED - not useful enough; will only cause confusion
        """Using CompositeHandler, takes a list of other HandlerLevels and
        appends the contents of each HandlerLevel to the corresponding
        element in its own list of Handlers, effectively doing a splice.
        I.e. does self.handlers[i].append(handler_levels[i]).  Returns a
        fresh HandlerLevel; does not self-modify.

        """
        new_handler_level = HandlerLevel(None, [])
        new_handlers = []
        for i, handler in enumerate(self.handlers):
            new_handlers.append(
                CompositeHandler(handler.level_name, handler.value, 
                                 handler_levels[i].handlers))
        new_handler_level.handlers = new_handlers
        return new_handler_level

    def get_handlers(self):
        return self.handlers

    def handle(self, request):
        """Simple delegation"""
        for handler in self.handlers:
            handler.handle(request)
    
class TestSuite:
    def __init__(self, head_handler):
        self.head_handler = head_handler
