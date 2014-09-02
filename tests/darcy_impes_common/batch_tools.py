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

## since it is likely that the client will want to give all the classes
# in this module the same verbosity, verbosity may as well be a module
# variable.  The TestSuite class below can be used to set/get the
# variable.  It is an integer to allow different levels of verbosity.
_global_verbosity = 1
_indent_str = '   '

def set_global_verbosity(verbosity_level):
    global _global_verbosity
    _global_verbosity = verbosity_level
    
    
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


# # [see note 1]
# 
# class Handler:
#     """A composite to enable recursion over arbitrary combinations of
#     simulation options (e.g. domain dimensions, grid resolutions).  It
#     varies from the conventional Composite pattern in that the
#     parent-child functionality is shared between 'Level' and 'Node'
#     components, and these components refer to instances of one another
#     rather than the abstract class.
#     """
#     def set_verbose(self, verbose):
#         """Sets verbosity."""
#         pass
#     def set_level_name(self, level_name):
#         """Assigns a level name."""
#         pass
#     def set_level_index(self, level_index):
#         """Recurses through the tree assigning level numbers."""
#         pass
#     def handle(self, command):
#         """Recurses through the tree, executing command at each node
#         (whether branch or leaf)."""
#         pass
#     def make_level(self):
#         """Helper for add_sub.  Returns a new (Level)Handler and
#         deep-copies any components for safety."""
#         pass
#     def add_sub(self, other):
#         """Adds a further level of nodes to each leaf node in the present tree.
#         """
#         pass
#     def __mul__(self, other):
#         """Like add_sub, but provides a new Handler instead of modifying
#         the current object.  This means the function can be inlined neatly.
#         """
#         result = deepcopy(self)
#         result.add_sub(other)
#         return result
#     def get_level_names(self, level_names=None):
#         """Recurses down the data structure recording level names as
#         they appear in the hierarchy.  Returns an ordered list.
#         Does not check consistency of ordering."""
#         pass

    
class _Handler:
    """A composite to enable recursion over arbitrary combinations of
    simulation options (e.g. domain dimensions, grid resolutions).  It
    varies from the conventional Composite pattern in that the
    parent-child functionality is shared between 'Level' and 'Node'
    components, and these components refer to instances of one another
    rather than the abstract class.
    """
    def __mul__(self, other):
        """Like add_sub in the subclasses that follow, but provides a
        new Handler instead of modifying the current object.  This means
        the function can be inlined neatly.
        """
        result = deepcopy(self)
        result.add_sub(other)
        return result

    
class _LevelHandler(_Handler):
    """A lightweight class containing a list of nodes at a certain
    level.  Most of the methods are just delegations to the nodes."""

    def __init__(self, level_name, node_list):
        # copy node_list to prevent side effects
        self.__node_list = deepcopy(node_list)
        # the supplied nodes may be strings, in which case they must be
        # converted into NodeHandlers first.
        for i, node in enumerate(self.__node_list):
            if not isinstance(node, _NodeHandler):
                self.__node_list[i] = _NodeHandler(node)

        self.set_level_name(level_name)
        # default to level 0.  This will change if the present handler
        # gets assigned a parent.
        self.set_level_index(0)

    def set_level_name(self, level_name):
        self.__level_name = level_name
        # it is convenient to have the nodes store level_name
        for node in self.__node_list:
            node.set_level_name(level_name)
                
    def set_level_index(self, level_index):
        self.level_index = level_index
        # iterate and delegate
        for node in self.__node_list:
            node.set_level_index(level_index)
    
    def handle(self, command, verbose=(_global_verbosity>0)):
        # iterate and delegate
        for node in self.__node_list:
            node.handle(command, verbose)

    def make_level(self):
        # already a level, so just clone self and return
        return deepcopy(self)
    
    def add_sub(self, other):
        # iterate and delegate
        for node in self.__node_list:
            node.add_sub(other)

    def get_level_names(self, level_names=None):
        # iterate and delegate
        for node in self.__node_list:
            level_names = node.get_level_names(level_names)
        return level_names
        
        
class _NodeHandler(_Handler):
    """Has zero or one LevelHandler children.  If zero, acts as a leaf
    component.  If one, has a sub-level that will be recursed over."""
    
    def __init__(self, node_name, sublevel=None):
        self.__level_name = None
        self.__node_name = node_name
        # default to level 0.  This will change if the present handler
        # gets assigned a parent.
        self.__sublevel = sublevel
        self.set_level_index(0)
    
    def set_level_name(self, level_name):
        self.__level_name = level_name
    
    def set_level_index(self, level_index):
        self.level_index = level_index
        # set characters to help print diagnostics nicely
        if level_index == 0:
            self.__beg_char = ''
            self.__end_char = '\n'
        else:
            self.__beg_char = '\n'
            self.__end_char = ''
        # go to next level if it exists
        try:
            self.__sublevel.set_level_index(level_index + 1)
        except:
            pass

    def handle(self, command, verbose=(_global_verbosity>0)):
        at_leaf = not isinstance(self.__sublevel, _Handler)
        if verbose and self.__node_name:
            if self.__level_name:
                lev = self.__level_name + ': '
            else:
                lev = ''
            # start printing diagnostics
            sys.stdout.write(self.__beg_char+self.level_index*_indent_str + \
                                 lev + str(self.__node_name))
            sys.stdout.flush()
        # prep the command and try executing it at this level
        command.inform(self.__level_name, self.__node_name, at_leaf, \
                       self.level_index*_indent_str, verbose)
        command.execute()
        # go to next level if it exists
        if not at_leaf:
            self.__sublevel.handle(command, verbose)
        if verbose and self.__node_name:
            # finish printing diagnostics
            sys.stdout.write(self.__end_char)
        sys.stdout.flush()
        
    def make_level(self):
        # return a LevelHandler containing a copy of this one node
        return _LevelHandler(self.__level_name, [deepcopy(self)])

    def add_sub(self, other):
        # assign a new sublevel - or, if one already exists, recurse
        # downwards
        if self.__sublevel is None:
            self.__sublevel = other.make_level()
            # reset the subtree's level indices
            self.__sublevel.set_level_index(self.level_index + 1)
        else:
            self.__sublevel.add_sub(other)

    def get_level_names(self, level_names=None):
        if level_names==None:
            level_names = []
        if self.__level_name not in level_names:
            # needs to be in the list, so insert here
            level_names.insert(self.level_index, self.__level_name)
        # recurse downwards
        if self.__sublevel is not None:
            level_names = self.__sublevel.get_level_names(level_names)
        return level_names

def make_string(string_list):
    """Formats and joins a bunch of strings """
    result = ''
    for s in string_list:
        # any nil values are ignored
        if s is None:
            continue
        # convert any non-strings
        if s is not str:
            s = str(s)
        # chain the values together with underscores
        if result=='' or s=='':
            link = ''
        else:
            link = '_'
        result = link.join([result, s])
        # for safety (because keys will form filenames), replace
        # any points with 'p'
        result = sub('\.', 'p', result)
    return result

    
def new_handler(arg1, arg2=None, arg3=None):
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
        # case they must be converted into NodeHandlers first.  Copy
        # arg2 to prevent side effects.
        arg2 = deepcopy(arg2)
        for i, node in enumerate(arg2):
            if not isinstance(node, _NodeHandler):
                arg2[i] = _NodeHandler(node)
        new_handler = _LevelHandler(arg1, arg2)
    else:
        # we are creating a node
        if isinstance(arg2, str):
            # arg1 is the level name, arg2 is the node name and there
            # may be a sublevel in arg3
            new_handler = _NodeHandler(arg2, arg3)
            new_handler.set_level_name(arg1)
        else:
            # arg1 is the node name and there may be a sublevel in arg2
            new_handler = _NodeHandler(arg1, arg2)
            
    return new_handler


# # [see note 1]
# 
# class Command:
#     """Class for encapsulating a command corresponding to a tree of
#     Handlers.  The body of the command is coded up for different levels
#     of the tree in the 'execute' method.
#     """
#     def inform(self, level_name, node_name, at_leaf, indent, verbose):
#         """Passes information to the Command about where it is in the
#         handler tree.  The Command can then choose what to do upon
#         execution. """
#         pass
#     def execute(self):
#         """Triggers the command."""
#         pass
#     def handle(self, other=None, message=None):
#         """Enables handling of other commands or handling of self (see
#         below)."""
#         pass
    

class Tracker:
    """Helper class for command implementations.  Clients may want to
    store an instance of this class in each command.  It manages
    information as to where the command is in the handler tree.  It also
    has a method for making a unique ID.
    """
    
    def __init__(self, key_formatting_dict={}):
        self.__level_name = ''
        self.__level_dict = OrderedDict()
        self.__at_root = None
        self.__at_leaf = None
        self.__key_formatting_dict = key_formatting_dict

    def inform(self, level_name, node_name, at_leaf, indent):
        """Commands should delegate inform(...) to this method."""
        # record the level-node pair and reset all below it
        self.__level_name = level_name
        self.__level_dict[level_name] = node_name
        while self.__level_dict.keys()[-1] != level_name:
            self.__level_dict.popitem()
        self.__at_root = (level_name == self.__level_dict.keys()[0])
        self.__at_leaf = at_leaf

    def update(self, level_name, node_name):
        """Modifies a key-value entry in the embedded dictionary"""
        self.__level_dict[level_name] = node_name

    def get(self, what):
        """Based on the value of the argument, returns the following:
        'root'  - first level
        'level' - current level
        'node'  - current node
        anything else - assumes the argument is a level name from which
        a node name is to be returned.  Return value defaults to None.
        """
        if what is 'root':
            return self.__level_dict.values()[0]
        if what is 'level':
            return self.__level_name
        elif what is 'node':
            return self.__level_dict[self.__level_name]
        else:
            try:
                return self.__level_dict[what]
            except KeyError as e:
                e.args += (
                    "The handler may not have reached level \'{0}\' yet.  ".\
                        format(what) + \
                        "So far, the levels are {0}.  ".\
                        format(self.__level_dict.keys()) + \
                        "Go up the stack and modify the level query.",)
                raise
            
    def at(self, where):
        """Returns True if where=='root' and we are at the base of the
        handler tree, if where=='leaf' and we are at the end of a
        branch, or if where is the current level_name.
        """
        return \
            (where=='root' and self.__at_root) or \
            (where=='leaf' and self.__at_leaf) or \
            (where==self.__level_name)
    
    def make_key(self, level_names=None, exclude=[], substitute={}):
        """Forms a string made up of node names corresponding to the
        supplied list of level names.  The list defaults to all the keys
        in the embedded dictionary.  Optional list argument 'exclude'
        specifies levels to be skipped over.  Optional dictionary
        'substitute' specifies levels to be overridden with alternative
        node values.
        """
        # default level names
        if level_names is None:
            level_names = self.__level_dict.keys()
        node_names = []
        
        # the string list args can contain the special word 'root',
        # which is synonymous with the first level name.
        def correct(string_list):
            try:
                i = string_list.index('root')
                string_list[i] = self.__level_dict.keys()[0]
            except:
                pass
            return string_list
        level_names = correct(level_names)
        exclude_names = correct(level_names)

        # make a list of node names
        for ln in level_names:
            # any nonpresent or excluded values are ignored
            if ln not in self.__level_dict or ln in exclude:
                continue
            if ln in substitute.keys():
                nn = substitute[ln]
            else:
                nn = self.__level_dict[ln]
            # annotate with extra formatting if specified
            if ln in self.__key_formatting_dict.keys():
                nn = self.__key_formatting_dict[ln].format(nn)
            node_names.append(nn)

        # pass node_names to the make_string free function
        return make_string(node_names)

    
class Command:
    """Class for encapsulating a command corresponding to a tree of
    Handlers.  The body of the command is coded up for different levels
    of the tree in the 'execute' method.
    """
    def __init__(self, verbosity_threshold=1):
        self.__verbosity_threshold = verbosity_threshold
        self.__verbose = _global_verbosity >= verbosity_threshold
        self.__indent = ''

    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        """Passes information to the Command about where it is in the
        handler tree.  The Command can then choose what to do upon
        execution. """
        self.__verbose = verbose
        self.__indent = indent

    def is_verbose(self):
        return self.__verbose and \
            _global_verbosity >= self.__verbosity_threshold

    def write_message(self, message, target=sys.stdout, with_newline=False):
        """For verbose printing purposes."""
        if self.is_verbose():
            if with_newline:
                target.write('\n')
                target.write(self.__indent)
            else:
                target.write(' ')
            target.write(message)
            target.flush()

            
class Validate(Command):
    """A supporting command that ensures compatibility with a handler."""
    
    def __init__(self, requisite_level_names, client_name):
        """client_name refers to the object that is sending the Validate
        object to a handler in order to test compatibility.
        requisite_level_names is a list of level names that must be
        encountered, and in the given order.  If 'root' appears at the
        beginning of this list, it signifies that there should only be
        one node at the top level."""
        Command.__init__(self, verbosity_threshold=2)
        self.__client_name = client_name
        self.__tracker = Tracker()
        self.__requisite_level_names = requisite_level_names
        self.__checklist = [False for i in range(len(requisite_level_names))]
        self.__root_node_countdown = -1
        try:
            if self.__requisite_level_names[0]=='root':
                self.__requisite_level_names.pop(0)
                self.__root_node_countdown = 1
        except:
            pass
        
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        # delegate to superclass and tracker component
        Command.inform(self, level_name, node_name, at_leaf, indent, verbose)
        self.__tracker.inform(level_name, node_name, at_leaf, indent)

    def execute(self):
        # do single root node check 
        if self.__root_node_countdown == 0:
            raise RuntimeError(
                "There should only be one top-level node in the handler for {0}".\
                    format(self.__client_name))
        if self.__tracker.at('root'):
            self.__root_node_countdown -= 1
            
        # check off this level if it appears in the list
        l = self.__tracker.get('level')
        if l in self.__requisite_level_names:
            n = self.__requisite_level_names.index(l)
            self.__checklist[n] = True
        else:
            n = -1

        if self.__tracker.at('leaf'):
            # we have reached the end of a branch.  Are all the
            # requisite levels present?
            for i, c in enumerate(self.__checklist):
                if not c:
                    l = self.__requisite_level_names[i]
                    raise RuntimeError(
                        "'{0}' is required in the handler for {1}".format(
                            l, self.__client_name))

        else:
            # we are at a branch node.  If the current level appeared in
            # the list, does it appear in the correct order?  Check the
            # section of requisite_level_names traversed so far
            if n >= 0:
                for i, c in enumerate(self.__checklist[0:n+1]):
                    if not c:
                        raise RuntimeError(
                            ("'{0}' is required to appear before '{1}' "+\
                                 "in the handler for {2}").format(
                                    self.__requisite_level_names[i],
                                    self.__tracker.get('level'),
                                    self.__client_name))

        if n > -1 and self.is_verbose():
            sys.stdout.write(" ... OK")

        # when done, clear all below this level from the list
        while n > -1 and n+1 < len(self.__checklist):
            self.__checklist[n+1] = False
            n += 1


class SelfHandlingCommand(Command):
    """This command is associated with a handler.  The client can
    initiate the handle-execute process by calling handle()."""

    def __init__(self, requisite_level_names, handler,
                 message='', verbosity_threshold=1 ):
        Command.__init__(self, verbosity_threshold)
        self.__message = message
        self.__handler = handler
        self.validate(requisite_level_names)
            
    def validate(self, requisite_level_names):
        val = Validate(requisite_level_names, self.__class__.__name__)
        if val.is_verbose():
            msg = 'Validating {0}'.format(self.__class__.__name__)
        else:
            msg = ''
        self.handle(val, msg)

    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        Command.inform(self, level_name, node_name, at_leaf, indent, verbose)

    def handle(self, other=None, message=None):
        """Sends itself (or some other command) to its associated
        handler."""
        if message is None:
            message = self.__message
        if message:
            self.write_message('\n'+message+'\n')
        if other is None:
            self.__handler.handle(self, self.is_verbose())
        else:
            self.__handler.handle(other, other.is_verbose())

    
# class SelfHandlingCommandList:
#     """Chains multiple (self-handling) commands together.  Simply
#     iterates and delegates for each method."""

#     def __init__(self, commands):
#         self.commands = commands

#     def inform(self, level_name, node_name, at_leaf, indent, verbose):
#         for cmd in self.commands:
#             cmd.inform(level_name, node_name, at_leaf, indent, verbose)

#     def execute(self):
#         for cmd in self.commands:
#             cmd.execute()

#     def handle(self, other=None, message=None):
#         for cmd in self.commands:
#             cmd.handle()

#     def write_message(self, message, target=sys.stdout, with_newline=False):
#         for cmd in self.commands:
#             cmd.write_message(message, target, with_newline)
    

class DoNothing(Command):
    """A blank SelfHandlingCommand for use in stubs etc."""
    def __init__(self, message='', verbosity_threshold=1):
        Command.__init__(self, verbosity_threshold)
        self.__message = message
            
    def validate(self, requisite_level_names):
        pass
    
    def inform(self, level_name, node_name, at_leaf, indent, verbose):
        pass
    
    def handle(self, other=None, message=None):
        if message is None:
            message = self.__message
        if message:
            self.write_message('\n'+message+'\n')

    
# [1] Python's duck typing renders abstract classes and methods
#     redundant.  They will be phased out completely once I have
#     strengthened the documentation.
