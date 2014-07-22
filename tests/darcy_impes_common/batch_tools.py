#!/usr/bin/env python

"""This module mainly implements a Command-Composite pattern for
executing arbitrary commands over arbitrary combinations of options.

More documentation to come.
"""

import sys
import os
from copy import deepcopy
from collections import OrderedDict
from re import sub


## HELPERS

# in this module, 'verbosity' is an integer; 'verbose' is True/False.
default_verbosity = 1
indent = '   '

class Handler:
    """A composite to enable recursion over arbitrary combinations of
    simulation options (e.g. domain dimensions, grid resolutions).  It
    varies from the conventional Composite pattern in that the
    parent-child functionality is shared between 'Level' and 'Node'
    components, and these components refer to instances of one another
    rather than the abstract class.

    """

    def set_verbose(self, verbose):
        """Sets verbosity."""
        pass

    def set_level_name(self, level_name):
        """Assigns a level name."""
        pass

    def set_level_index(self, level_index):
        """Recurses through the tree assigning level numbers."""
        pass
        
    def handle(self, command):
        """Recurses through the tree, executing command at each node
        (whether branch or leaf)."""
        pass

    def make_level(self):
        """Helper for add_sub.  Returns a new (Level)Handler and
        deep-copies any components for safety."""
        pass

    def add_sub(self, other):
        """Adds a further level of nodes to each leaf node in the present tree.
        """
        pass

    def __mul__(self, other):
        """Like add_sub, but provides a new Handler instead of modifying
        the current object.  This means the function can be inlined neatly.
        """
        result = deepcopy(self)
        result.add_sub(other)
        return result

    def get_level_names(self, level_names=None):
        """Recurses down the data structure recording level names as
        they appear in the hierarchy.  Returns an ordered list.
        Does not check consistency of ordering."""
        pass

    
class LevelHandler(Handler):
    """A lightweight class containing a list of nodes at a certain
    level.  Most of the methods are just delegations to the nodes."""

    def __init__(self, level_name, node_list):
        self.node_list = node_list
        # the supplied nodes may be strings, in which case they must be
        # converted into NodeHandlers first.
        for i, node in enumerate(node_list):
            if not isinstance(node, NodeHandler):
                node_list[i] = NodeHandler(node)

        self.set_level_name(level_name)
        # default to level 0.  This will change if the present handler
        # gets assigned a parent.
        self.set_level_index(0)

    def set_verbose(self, verbose):
        # iterate and delegate
        for node in self.node_list:
            node.set_verbose(verbose)
        
    def set_level_name(self, level_name):
        self.level_name = level_name
        # it is convenient to have the nodes store level_name
        for node in self.node_list:
            node.set_level_name(level_name)
                
    def set_level_index(self, level_index):
        self.level_index = level_index
        # iterate and delegate
        for node in self.node_list:
            node.set_level_index(level_index)
    
    def handle(self, command):
        # iterate and delegate
        for node in self.node_list:
            node.handle(command)

    def make_level(self):
        # already a level, so just clone self and return
        return deepcopy(self)
    
    def add_sub(self, other):
        # iterate and delegate
        for node in self.node_list:
            node.add_sub(other)

    def get_level_names(self, level_names=None):
        # iterate and delegate
        for node in self.node_list:
            level_names = node.get_level_names(level_names)
        return level_names
        
        
class NodeHandler(Handler):
    """Has zero or one LevelHandler children.  If zero, acts as a leaf
    component.  If one, has a sub-level that will be recursed over."""
    
    def __init__(self, node_name, sublevel=None, 
                 verbose=(default_verbosity>0)):
        self.level_name = None
        self.node_name = node_name
        # default to level 0.  This will change if the present handler
        # gets assigned a parent.
        self.sublevel = sublevel
        self.set_level_index(0)
        self.set_verbose(verbose)

    def set_verbose(self, verbose):
        self.verbose = verbose
        # go to next level if it exists
        try:
            self.sublevel.set_verbose(verbose)
        except:
            pass
    
    def set_level_name(self, level_name):
        self.level_name = level_name
    
    def set_level_index(self, level_index):
        self.level_index = level_index
        # set characters to help print diagnostics nicely
        if level_index == 0:
            self.beg_char = ''
            self.end_char = '\n'
        else:
            self.beg_char = '\n'
            self.end_char = ''
        # go to next level if it exists
        try:
            self.sublevel.set_level_index(level_index + 1)
        except:
            pass

    def handle(self, command):
        at_leaf = not isinstance(self.sublevel, Handler)
        if self.verbose and self.node_name:
            if self.level_name:
                lev = self.level_name + ': '
            else:
                lev = ''
            # start printing diagnostics
            sys.stdout.write(self.beg_char+self.level_index*indent + \
                                 lev + str(self.node_name))
            sys.stdout.flush()
        # prep the command and try executing it at this level
        command.inform(self.level_name, self.node_name, at_leaf, \
                       self.level_index*indent, self.verbose)
        command.execute()
        # go to next level if it exists
        if not at_leaf:
            self.sublevel.handle(command)
        if self.verbose and self.node_name:
            # finish printing diagnostics
            sys.stdout.write(self.end_char)
        sys.stdout.flush()
        
    def make_level(self):
        # return a LevelHandler containing a copy of this one node
        return LevelHandler(self.level_name, [deepcopy(self)])

    def add_sub(self, other):
        # assign a new sublevel - or, if one already exists, recurse
        # downwards
        if self.sublevel is None:
            self.sublevel = other.make_level()
            # reset the subtree's level indices
            self.sublevel.set_level_index(self.level_index + 1)
        else:
            self.sublevel.add_sub(other)

    def get_level_names(self, level_names=None):
        if level_names==None:
            level_names = []
        if self.level_name not in level_names:
            # needs to be in the list, so insert here
            level_names.insert(self.level_index, self.level_name)
        # recurse downwards
        if self.sublevel is not None:
            level_names = self.sublevel.get_level_names(level_names)
        return level_names
    

## FOR PUBLIC USE 
    
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


def make_string(string_list):
    """Formats and joins a bunch of strings """
    result = None
    for s in string_list:
        # any nil values are ignored
        if s is None:
            continue
        # convert any non-strings
        if s is not str:
            s = str(s)
        # chain the values together with underscores
        if result is None:
            result = s
        else:
            result = result + '_' + s
        # for safety (because keys will form filenames), replace
        # any points with 'p'
        result = sub('\.', 'p', result)
    return result

    
def new_handler(arg1, arg2=None, arg3=None, verbose=(default_verbosity>0)):
    """Creates a Handler based on a flexible set of arguments.  Either a
    LevelHandler or NodeHandler will be returned.  The possibilities
    are:
       Create a level:
          new_handler(level_name, node_list)
          new_handler(level_name, node_name_list)
       Create a node:
          new_handler(node_name)
          new_handler(node_name, sublevel)
          new_handler(level_name, node_name)
          new_handler(level_name, node_name, sublevel)
    """

    # are we creating a level or a node?
    is_level = False
    if arg2 is not None:
        if isinstance(arg2, list):
            # we are creating a level and arg2 is the list of nodes
            is_level = True

    if is_level:
        # the supplied nodes may be represented as strings, in which
        # case they must be converted into NodeHandlers first.
        for i, node in enumerate(arg2):
            if not isinstance(node, NodeHandler):
                arg2[i] = NodeHandler(node)
        new_handler = LevelHandler(arg1, arg2)
    else:
        # we are creating a node
        if isinstance(arg2, str):
            # arg1 is the level name, arg2 is the node name and there
            # may be a sublevel in arg3
            new_handler = NodeHandler(arg2, arg3)
            new_handler.set_level_name(arg1)
        else:
            # arg1 is the node name and there may be a sublevel in arg2
            new_handler = NodeHandler(arg1, arg2)
            
    new_handler.set_verbose(verbose)
    return new_handler

        
class Command:
    """Purely abstract base class for encapsulating a command
    corresponding to a tree of Handlers.  The body of the command is
    coded up for different levels of the tree in the 'execute' method.
    """
    
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        """Passes information to the Command about where it is in the
        handler tree.  The Command can then choose what to do upon
        execution. """
        pass

    def execute(self):
        """Triggers the command."""
        pass

    def handle(self, other=None, message=None):
        """Enables handling of other commands or handling of self (see
        below)."""
        pass


class Tracker:
    """Helper class for Command implementations.  Clients may want to
    store an instance of this class in each concrete Command.  It
    manages information as to where the command is in the handler tree.
    It also has a method for making a unique ID.
    """
    
    def __init__(self, key_formatting_dict={}):
        self.level_name = ''
        self.level_dict = OrderedDict()
        self.at_root = None
        self.at_leaf = None
        self.key_formatting_dict = key_formatting_dict

    def inform(self, level_name, node_name, at_leaf, indent):
        """Commands should delegate inform(...) to this method."""
        # record the level-node pair and reset all below it
        self.level_name = level_name
        self.level_dict[level_name] = node_name
        while self.level_dict.keys()[-1] != level_name:
            self.level_dict.popitem()
        self.at_leaf = at_leaf
        self.at_root = (indent=='') 

    def update(self, level_name, node_name):
        """Modifies a key-value entry in the embedded dictionary"""
        self.level_dict[level_name] = node_name

    def get(self, what):
        """If what=='level', returns current level name.  If what=='node',
        returns current node name.  Otherwise assumes the argument is a
        level name from which a node name is to be returned.  Return
        value defaults to None.

        """
        if what is 'level':
            return self.level_name
        elif what is 'node':
            return self.level_dict[self.level_name]
        else:
            try:
                return self.level_dict[what]
            except KeyError:
                return None

    def at(self, where):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        return \
            (where=='root' and self.at_root) or \
            (where=='leaf' and self.at_leaf) or \
            (where==self.level_name)
    
    def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.  The list defaults to all the keys
        in the embedded dictionary.  Optional list argument 'exclude'
        specifies levels to be skipped over.  Optional dictionary
        'substitute' specifies levels to be overridden with alternative
        node values.
        """
        if level_names is None:
            level_names = self.level_dict.keys()
        node_names = []
        for ln in level_names:
            # any nonpresent or excluded values are ignored
            if ln not in self.level_dict or ln in exclude:
                continue
            if ln in substitute.keys():
                nn = substitute[ln]
            else:
                nn = self.level_dict[ln]
            # annotate with extra formatting if specified
            if ln in self.key_formatting_dict.keys():
                nn = self.key_formatting_dict[ln].format(nn)
            node_names.append(nn)
        # pass node_names to the make_string free function
        return make_string(node_names)

            
class Validate(Command):
    """A supporting Command that ensures compatibility with a handler."""
    
    def __init__(self, client_name, requisite_level_names):
        """client_name refers to the TestCommand object that is sending
        the Validate object to a handler in order to test compatibility.
        requisite_level_names is a list of level names that must be
        encountered, and in the given order."""
        self.tracker = Tracker()
        self.client_name = client_name
        self.requisite_level_names = requisite_level_names
        self.checklist = [False for i in range(len(requisite_level_names))]
        self.root_node_count = 0

    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        # delegate
        self.tracker.inform(level_name, node_name, at_leaf, indent)
        self.verbose = verbose

    def execute(self):
        # there should only be one node at the root corresponding to
        # level 'case' or whatever is topmost.  This will be used by
        # TestCommand objects to form filenames.
        if self.tracker.at('root'):
            self.root_node_count += 1
        if self.root_node_count > 1:
            raise RuntimeError(
                "There should only be one top-level node in the handler for {0}".\
                    format(self.client_name))
        
        # check off this level if it appears in the list
        l = self.tracker.get('level')
        if l in self.requisite_level_names:
            n = self.requisite_level_names.index(l)
            self.checklist[n] = True
        else:
            n = -1

        if self.tracker.at('leaf'):
            # we have reached the end of a branch.  Are all the
            # requisite levels present?
            for i, c in enumerate(self.checklist):
                if not c:
                    l = self.requisite_level_names[i]
                    raise RuntimeError(
                        "'{0}' is required in the handler for {1}".format(
                            l, self.client_name))

        else:
            # we are at a branch node.  If the current level appeared in
            # the list, does it appear in the correct order?  Check the
            # section of requisite_level_names traversed so far
            if n >= 0:
                for i, c in enumerate(self.checklist[0:n+1]):
                    if not c:
                        raise RuntimeError(
                            ("'{0}' is required to appear before '{1}' "+\
                                 "in the handler for {2}").format(
                                    self.requisite_level_names[i],
                                    self.tracker.get('level'),
                                    self.client_name))

        if n > -1 and self.verbose:
            sys.stdout.write("   ---OK")

        # when done, clear all below this level from the list
        while n > -1 and n+1 < len(self.checklist):
            self.checklist[n+1] = False
            n += 1


class SelfHandlingCommand(Command):
    """This abstract class is associated with a Handler.  The
    Handler is validated so that the client doesn't have to worry so
    much about Handler-Command compatibility.  The client can
    initiate the handle-execute process by calling handle()."""

    def __init__(self, requisite_level_names, handler, message,
                 verbosity=default_verbosity):
        self.verbosity = verbosity
        self.requisite_level_names = requisite_level_names
        self.associate(handler)
        self.level_names = handler.get_level_names([])
        self.message = message
        self.indent = ''

    def handle(self, other=None, message=None):
        """Sends itself (or some other command) to its associated
        handler."""
        if message is None:
            message = self.message
        if message:
            self.write_message('\n'+message+'\n')
        if other is None:
            self.handler.handle(self)
        else:
            self.handler.handle(other)
            
    def associate(self, handler):
        """Associates a handler with the current command."""
        # validate
        self_name = self.__class__.__name__
        # if verbosity is high, print info about the validation
        handler.set_verbose(self.verbosity > 1)
        if self.verbosity > 1:
            self.write_message("\nValidating " + self_name + " handler\n")
        handler.handle( Validate(self_name, self.requisite_level_names) )
        # if we have got here, validation was successful
        self.handler = handler
        self.handler.set_verbose(self.verbosity > 0)

    def write_message(self, message, target=sys.stdout, with_newline=False):
        """Called by the subclasses for verbose printing purposes."""
        if self.verbosity > 0:
            if with_newline:
                target.write('\n')
                target.write(self.indent)
            else:
                target.write(' ')
            target.write(message)
            target.flush()
    

class DoNothing(SelfHandlingCommand):
    """A blank for use in method stubs etc."""
    def __init__(self, message=''):
        SelfHandlingCommand.__init__(self, [], new_handler(''), message)
