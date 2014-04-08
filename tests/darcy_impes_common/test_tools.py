#!/usr/bin/env python

import sys
sys.path.append('../../python')

import os
os.environ['PYTHONPATH'] = "../../python/"

verbose = True
def set_verbose(new_verbose):
    global verbose
    verbose = new_verbose

class Command:
    """Class to encapsulate a family of commands corresponding to a tree of
    Handlers (see below).  The commands are coded up for different
    levels of the tree in the 'execute' method.
    """
    
    def execute(self, level_name, value):
        """The level_name argument tells the Command where it is in the Handler tree.
        The Command can then choose what to do accordingly.  It may use
        additional information supplied in the value argument.
        """
        pass


class Handler:
    """A composite to enable recursion over arbitrary combinations of
    simulation options (e.g. domain dimensions, grid resolutions).  The
    'handle' method combines the recursion with executing of a Command.
    """
    
    def __init__(self, level_name, value):
        self.level_name = level_name
        self.value = value
        self.level_index = 0
        
    def set_level_index(self, level_index):
        pass

    def handle(self, command):
        pass


class CompositeHandler(Handler):
    """A tree of handlers."""

    def __init__(self, level_name, value, children):
        Handler.__init__(self, level_name, value)
        try:
            # if children is a HandlerList (see below), convert it back
            self.children = children.get_handlers()
        except AttributeError:
            # already a Python list of Handlers
            self.children = children
        for child in self.children:
            child.set_level_index(self.level_index + 1)
        
    def set_level_index(self, level_index):
        self.level_index = level_index
        # recurse
        for child in self.children:
            child.set_level_index(self.level_index + 1)
        
    def handle(self, command):
        if verbose: print self.level_index*'  ' + \
           self.level_name + ': ' + str(self.value)
        # try executing the command at this level
        command.execute(self.level_name, self.value)
        # recurse
        for child in self.children:
            child.handle(command)

        
class LeafHandler(Handler):
    """A single handler."""

    def __init__(self, level_name, value):
        Handler.__init__(self, level_name, value)
        
    def set_level_index(self, level_index):
        self.level_index = level_index

    def handle(self, command):
        if verbose:
            self.level_name+': '+str(self.value)
        # try executing the command at this level
        command.execute(self.level_name, self.value)

        
class HandlerList:
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

    def expand(self, handler_list):
        """Using CompositeHandler, takes another HandlerList and appends copies
        of it to each element in its own list of Handlers, effectively
        multiplying the HandlerLists together.  I.e. does
        self.handlers[i].append(handler_list).  Returns a fresh
        HandlerList; does not self-modify.
        """
        new_handler_list = HandlerList(None, [])
        new_handlers = []
        for handler in self.handlers:
            new_handlers.append(
                CompositeHandler(handler.level_name, handler.value, 
                                 handler_list.handlers))
        new_handler_list.handlers = new_handlers
        return new_handler_list
    
    def splice(self, handler_lists):
        """Using CompositeHandler, takes a (Python) list of other HandlerLists
        and appends each HandlerList to the corresponding element in its
        own list of Handlers, effectively doing a splice.  I.e. does
        self.handlers[i].append(handler_lists[i]).  Returns a fresh
        HandlerList; does not self-modify.
        """
        new_handler_list = HandlerList(None, [])
        new_handlers = []
        for i, handler in enumerate(self.handlers):
            new_handlers.append(
                CompositeHandler(handler.level_name, handler.value, 
                                 handler_lists[i].handlers))
        new_handler_list.handlers = new_handlers
        return new_handler_list

    def get_handlers(self):
        return self.handlers

    def handle(self, request):
        """Simple delegation"""
        for handler in self.handlers:
            handler.handle(request)
    
