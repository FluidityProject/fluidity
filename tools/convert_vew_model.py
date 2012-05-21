#!/usr/bin/env python
import sys
import re
from pyparsing import Literal, Word, Group, OneOrMore, ZeroOrMore, alphanums, nums, ParseException, Forward, replaceWith, MatchFirst, Or
from xml.dom.minidom import parse

### taken from: http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split("([0-9]+)", key) ] 
    return sorted(l, key = alphanum_key)


### Base grammar ###
ident = Word( alphanums + "_" + "$" )
number = Word( nums + "." + "-" + "e" + "E" )
comma = Literal( ",").suppress()
lpar = Literal( "{" ).suppress()
rpar = Literal( "}" ).suppress()

expr = Forward()
bool_expr = Forward()
stmt = Forward()

assign = Literal(r"\assign")
val = Literal(r"\val")

var = Literal(r"\var")
variable = Group( var + lpar + ident + rpar )
assignment = Group( assign + lpar + variable + comma + expr + rpar )

si_value = Literal(r"\sival") + lpar + number + comma + number + rpar
unit = Literal(r"\unit") + lpar + number + comma + Word(alphanums) + comma + number + rpar
value = Group( val + lpar + si_value + comma + unit + rpar )

### Math expressions
math_ops = {
  r"\add"   : lambda t: "(" + t[0] + "".join( [ " + " + s for s in t[1:]] ) + ")",
  r"\mul"   : lambda t: "(" + t[0] + "".join( [ " * " + s for s in t[1:]] ) + ")",
  r"\sub"   : lambda t: "(" + t[0] + " - " + t[1] + ")",
  r"\div"   : lambda t: "(" + t[0] + " / " + t[1] + ")",
  r"\pow"   : lambda t: "math.pow(" + t[0] +", " + t[1] + ")",
  r"\max"   : lambda t: "max(" + t[0] +", " + t[1] + ")",
  r"\min"   : lambda t: "min(" + t[0] +", " + t[1] + ")",
  r"\exp"   : lambda t: "math.exp(" + t[0] + ")", 
  r"\log10" : lambda t: "numpy.log10(" + t[0] + ")", 
  r"\minus" : lambda t: "-(" + t[0] + ")", 
  r"\abs"   : lambda t: "abs( " + t[0] + " )", 
  r"\rnd"   : lambda t: "TODO RND(" + t[0] + ")"
}
math_op = MatchFirst( map(Literal, math_ops.keys()) )
math_expr = Group( math_op + lpar + expr + ZeroOrMore( comma + expr ) + rpar )

### Booleans and conditional expressions
bool_ops = {
  r"\equal"        : lambda t: "(" + t[0] + " == " + t[1] + ")",
  r"\neq"          : lambda t: "(" + t[0] + " != " + t[1] + ")",
  r"\greaterequal" : lambda t: "(" + t[0] + " >= " + t[1] + ")",
  r"\greater"      : lambda t: "(" + t[0] + " > " + t[1] + ")",
  r"\lessequal"    : lambda t: "(" + t[0] + " <= " + t[1] + ")",
  r"\less"         : lambda t: "(" + t[0] + " < " + t[1] + ")",
  r"\and"          : lambda t: "(" + t[0] + " and " + t[1] + ")",
  r"\or"           : lambda t: "(" + t[0] + " or " + t[1] + ")"
}
bool_op = Or( map(Literal, bool_ops.keys()) )
bool_expr << Group( bool_op + lpar + expr + comma + expr + rpar ) 

cond = Literal( r"\conditional" )
ifthen = Literal( r"\ifthen" )
conditional = Group( cond + lpar + bool_expr + comma + expr + comma + expr + rpar )
cond_stmt = Group( ifthen + lpar + bool_expr + comma + stmt + rpar )

# Special system variables
system_variables = {
    "TimeStep" : "dt_in_hours",
    "PI" : "math.pi",
    "Temp" : "env['Temperature']",
    "Vis_Irrad" : "env['Irradiance']",
    "Density" : "env['Density']",
    "S_t" : "vars['S_t']",
    "z" : "vars['z']",
    "IngestedCells" : "vars['PIngestedCells']",
    "MLDepth" : "param['MLDepth']",
    "Max_MLD" : "param['Max_MLD']",
    "d_year" : "param['d_year']"
}

# Variable class for converting names 
# according to variable type 
class Variable:
  def __init__(self, tok):
    global fgroup
    if tok[0] != var:
      raise Exception("Not a variable: " + tok)
    self.token = tok
    self.name = str(tok[1])
    self.state = False
    self.variety = False

    # VEW system variables 
    if self.name in system_variables.keys():
      self.name = system_variables[self.name]

    self.rhs = self.name
    self.lhs = self.name

    if self.name.endswith("$Pool"):
      self.name = self.name.split("$")[0]
      self.rhs = "vars['" + self.name + "']"
      self.lhs = self.name + "_new"
      self.state = True

    if self.name in fgroup.state_vars:
      self.lhs = self.name + "_new"
      self.rhs = "vars['" + self.name + "']"
      self.state = True

    if self.name.endswith("$Ingested"):
      self.name = self.name.split("$")[0] + "Ingested"
      self.rhs = "vars['" + self.name + "']"

    if self.name.endswith("$Conc"):
      self.name = self.name.split("$")[0]
      self.rhs = "env['Dissolved" + self.name + "']"

    if self.name in fgroup.parameters:
      self.rhs = "param['" + self.name + "']"

    if self.name in fgroup.variety_local:
      self.rhs = self.name + "[variety]"
      self.lhs = self.name + "[variety]"
      self.variety = True

    if self.name in fgroup.variety_conc:
      self.name = "env['" + self.name + "']"
      self.rhs = self.name + "[variety]"
      self.variety = True

    if self.name in fgroup.variety_param:
      self.rhs = "param[variety]['" + self.name + "']"
      self.variety = True

  def __eq__(self, other):
    return self.name == other.name



### Planktonica-specifics ###
#############################

divide = Literal( r"\divide" )
uptake = Literal( r"\uptake" )
release = Literal( r"\release" )
change = Literal( r"\change" )
pchange = Literal( r"\pchange" )
stage = Literal( r"\stage" )
ingest = Literal( r"\ingest" )
set = Literal( r"\set" )
varietysum = Literal( r"\varietysum" )
varhist = Literal( r"\varhist" )
visIrradAt = Literal( r"\visIrradAt" )
integrate = Literal( r"\integrate" )
create = Literal( r"\create" )
cell_division = Group( divide + lpar + value + rpar )
chemical_uptake = Group( uptake + lpar + expr + comma + variable + rpar )
chemical_release = Group( release + lpar + expr + comma + variable + rpar )
change_stage = Group( change + lpar + stage + lpar + ident + rpar + rpar )
pchange_stage = Group( pchange + lpar + stage + lpar + ident + rpar + comma + expr + rpar )
variety_vector_sum = Group( varietysum + lpar + variable + rpar )
variable_history = Group( varhist + lpar + variable + comma + value + rpar )
sample_irradiance = Group( visIrradAt + lpar + value + rpar )
integration = Group( integrate + lpar + expr + rpar )
ingestion = Group( ingest + lpar + expr + comma + expr + comma + expr + rpar )
set_variable = Group( set + lpar + variable + comma + expr + rpar )
creation = Group( create + lpar + stage + lpar + ident + rpar + comma + expr + OneOrMore( comma + set_variable ) + rpar )

expr << ( variable | value | bool_expr | math_expr | conditional |
          variety_vector_sum | variable_history | sample_irradiance | integration )
stmt << ( assignment | cond_stmt | cell_division | chemical_uptake | chemical_release | 
          change_stage | pchange_stage | ingestion | creation | set_variable )
planktonica = ( stmt )

fgroup = None

### Static indentation counter
class Indent:
  ind = 2
  @staticmethod
  def line():
    return "".join( [" " for i in range(Indent.ind)] )

  @staticmethod
  def inc():
    Indent.ind = Indent.ind + 2

  @staticmethod
  def dec():
    Indent.ind = Indent.ind - 2


### The big eval function....
def eval_expr(t):
  global fgroup
  out = ""

  if t[0] == var:
    v = Variable(t)
    if v.name in fgroup.local_vars.keys():
      fgroup.used_vars.append(v.rhs)
    out = v.rhs

  elif t[0] == val:
    base = float(t[2])
    e_pow = int(t[3])
    out = str(base)
    if e_pow != 0:
      out = out + "e" + str(e_pow)

  # Math functions
  elif t[0] in math_ops.keys():
    evt = [ eval_expr(tok) for tok in t[1:] ]
    out = math_ops[t[0]](evt)

  # Conditionals
  elif t[0] in bool_ops.keys():
    evt = [ eval_expr(tok) for tok in t[1:] ]
    out = bool_ops[t[0]](evt)
  elif t[0] == cond:
    out = "((" + eval_expr(t[2]) + ") if " + eval_expr(t[1]) + " else (" + eval_expr(t[3]) + "))"

  elif t[0] == varietysum:
    v = Variable(t[1])
    out = "sum(" + v.name + ".values())"

  elif t[0] == varhist:
    hist_ind = int(float(eval_expr(t[2])))
    v = Variable(t[1])
    out = "vars['" + v.name + "_" + str(hist_ind) + "']"

  elif t[0] == visIrradAt:
    out = "param['surface_irradiance']"

  # integration should be done by the framework
  elif t[0] == integrate:
    out = eval_expr(t[1])

  else:
    raise Exception("Unkown expression: " + str(t))

  return out


def eval_stmt(t):
  global fgroup
  out = ""

  if t[0] == assign:
    # When assigning state or pool variables, 
    # we create a buffer var and write the assignment at the end
    v = Variable(t[1])
    fgroup.initialised_vars.append(v.name)  

    if v.state and not v in fgroup.pool_update_vars:
      fgroup.pool_update_vars.append(v)      

    code = ""
    # Open a scope for variety-based calculations
    if v.variety:      
      code = code + v.name + " = {}\n" + Indent.line()
      Indent.inc()
      code = code + "for variety in env['P'].keys():\n" + Indent.line()

    fgroup.used_vars = []
    term = eval_expr(t[2])

    # Initialise all unitialised variables on the rhs
    for used_local in fgroup.used_vars:                
      if not used_local in fgroup.initialised_vars:
        code = code + used_local + " = 0.0\n" + Indent.line()
        fgroup.initialised_vars.append(used_local)

    if v.variety:
      Indent.dec()

    out = code + v.lhs + " = " + term

  elif t[0] == ifthen:
    Indent.inc()
    out = "if " + eval_expr(t[1]) + ":\n" + Indent.line() + eval_stmt(t[2])
    Indent.dec()

  # Planktonica
  elif t[0] == divide:
    out = "vars['Size'] = vars['Size'] * " + eval_expr(t[1])
  elif t[0] == change:
    s = str(t[2])
    out = "vars['Stage'] = stage_id('" + fgroup.name + "', '" + s + "')"
  elif t[0] == uptake:
    v = Variable(t[2])
    out = "vars['" + v.name + "Uptake'] = " + eval_expr(t[1])
  elif t[0] == release:
    v = Variable(t[2])
    out = "vars['" + v.name + "Release'] = " + eval_expr(t[1])

  elif t[0] == ingest:
    species_conc = eval_expr(t[1])
    ing_threshold = eval_expr(t[2])
    ing_amount = eval_expr(t[3])
    
    #code = "vars['P'] = {}\n" + Indent.line()
    Indent.inc()
    code = "for variety in vars['PRequest'].keys():\n" + Indent.line()
    Indent.dec()
    code = code + "vars['PRequest'][variety] = (dt * " + ing_amount + ") if (" + species_conc + " > " + ing_threshold + ") else 0.0"
    out = code

  elif t[0] == create:
    s = str(t[2])
    c = "new_agent_vars = {}\n"
    c = c + Indent.line() + "new_agent_vars['Stage'] = stage_id('" + fgroup.name + "', '" + s + "')\n"
    c = c + Indent.line() + "new_agent_vars['Size'] = " + eval_expr(t[3])
    for i in range(4,len(t)):
      c = c + "\n" + Indent.line() + "new_agent_" + eval_expr(t[i][1]) + " = " + eval_expr(t[i][2]) 
    c = c + "\n" + Indent.line() + "add_agent('" + fgroup.name + "', new_agent_vars, [-vars['z']])"
    out = c

  elif t[0] == pchange:
    # TODO This does not support multiple pchanges per update yet!!!
    s = str(t[2])
    out = "new_agent_vars = {}\n" + Indent.line()
    out = out + "new_agent_vars.update(vars)\n" + Indent.line()
    out = out + "new_agent_vars['Stage'] = stage_id('" + fgroup.name + "', '" + s + "')\n" + Indent.line()
    out = out + "new_agent_vars['Size'] = vars['Size'] * " + eval_expr(t[3]) + "\n" + Indent.line()
    out = out + "vars['Size'] = vars['Size'] - new_agent_vars['Size']\n" + Indent.line()
    out = out + "add_agent('" + fgroup.name + "', new_agent_vars, [-vars['z']])\n" + Indent.line()

  else:
    raise Exception("Unkown statement: " + str(t))

  return out
### End of big eval function


### Metamodel definitions ###
#############################

class Stage:
  counter = 0
  def __init__(self, dom_element):
    global id_counter
    self.dom = dom_element
    self.name = self.dom.getElementsByTagName("name")[0].firstChild.data
    self.id = Stage.counter
    self.functions = []
    self.function_eqns = {}
    self.assigned_locals = []
    Stage.counter = Stage.counter + 1

  def add_function(self, fname, function):
    self.functions.append(fname)
    self.function_eqns[fname] = function.getElementsByTagName("equation")

class Species:
  def __init__(self, dom_element):
    self.dom = dom_element
    self.name = self.dom.getAttribute("name").replace(' ', '_')
    self.fg_name = self.dom.getAttribute("fg").replace(' ', '_')

    self.parameter = {}
    for p in self.dom.getElementsByTagName("param"):
      p_name = p.getAttribute("name")
      p_val = p.getAttribute("a")
      self.parameter[p_name] = float(p_val)

class FoodSet:
  def __init__(self, dom_element):
    self.dom = dom_element
    self.species_name = self.dom.getAttribute("name").split(" : ")[0].replace(' ', '_')
    self.name = self.species_name + "_" + self.dom.getAttribute("name").split(" : ")[1]    

    self.food = {}
    self.target_species = {}
    for f in self.dom.getElementsByTagName("food"):
      stage = f.getAttribute("stage")
      self.food[stage] = {}
      for p in f.getElementsByTagName("param"):
        p_name = p.getAttribute("name")
        p_val = p.getAttribute("a")
        self.food[stage][p_name] = float(p_val)      

class FGroup:
  def __init__(self, dom_element):
    self.dom = dom_element
    self.name = self.dom.getElementsByTagName("name")[0].firstChild.data.replace(" ", "_")
    print "FGroup: " + self.name

    self.species = []
    self.foodsets = []

    self.parameters = []
    for parameter in self.dom.getElementsByTagName("parameter"):
      param_name = parameter.getElementsByTagName("name")[0].firstChild.data
      self.parameters.append(param_name)

    self.state_vars = []
    for state_var in self.dom.getElementsByTagName("variable"):
      varname = state_var.getElementsByTagName("name")[0].firstChild.data
      self.state_vars.append(varname)

    self.local_vars = {}
    for local in self.dom.getElementsByTagName("local"):
      varname = local.getElementsByTagName("name")[0].firstChild.data
      valtag = local.getElementsByTagName("value")
      varval = 0.0
      if len(valtag) > 0:
        varval = float(valtag[0].firstChild.data)
      self.local_vars[varname] = varval

    self.variety_param = {}
    for vparam in self.dom.getElementsByTagName("varietyparameter"):
      varname = vparam.getElementsByTagName("name")[0].firstChild.data
      val_str = vparam.getElementsByTagName("value")[0].firstChild.data
      self.variety_param[varname] = float(val_str)

    self.variety_local = []
    for vlocal in self.dom.getElementsByTagName("varietylocal"):
      varname = vlocal.getElementsByTagName("name")[0].firstChild.data
      self.variety_local.append(varname)

    self.variety_conc = []
    for vlocal in self.dom.getElementsByTagName("varietyconcentration"):
      varname = vlocal.getElementsByTagName("name")[0].firstChild.data
      self.variety_conc.append(varname)

    self.stages = {}
    for s in self.dom.getElementsByTagName("stage"):
      sname = s.getElementsByTagName("name")[0].firstChild.data
      self.stages[sname] = Stage(s)

    for function in fg.getElementsByTagName("function"):
      fname = function.getElementsByTagName("name")[0].firstChild.data
      for stage in function.getElementsByTagName("calledin"):
        sname = stage.firstChild.data
        self.stages[sname].add_function(fname, function)

  def write_parameters(self, file):
    file.write("\n# Parameters for FGroup " + self.name)

    # Write species parameter
    for s in self.species:
      file.write("\n# Species: " + s.name )
      file.write("\nspecies_" + s.name + " = {\n")
      for p in sorted_nicely( s.parameter.keys() ):
        file.write("    '" + p + "' : " + str(s.parameter[p]) + ",\n")
      file.write("}\n")

    for fs in self.foodsets:
      file.write("# Foodset: " + fs.name)
      file.write("\nfoodset_" + fs.name + " = {\n")
      for f_stage, food in fs.food.iteritems():
        file.write("    '" + f_stage + "' : {")
        for p in sorted_nicely( food.keys() ):
          file.write("\n        '" + p + "' : " + str(food[p]) + ",")
        file.write("    },\n")
      file.write("}\n")

  def write_update_kernel(self, file, stage):

    print "  Writing Stage: " + stage.name
    file.write("\ndef update_" + stage.name + "_" + self.name + "(param, vars, env, dt):\n")
    file.write('  """ FGroup:  ' + self.name + '\n')
    file.write('      Stage:   ' + stage.name + '\n  """\n')
    file.write("  dt_in_hours = dt / 3600.0\n")

    self.pool_update_vars = []
    self.initialised_vars = []    
    for fname in stage.functions:
      #print "  Function: " + fname
      file.write("\n  ### " + fname + " ###\n")
      for equation in stage.function_eqns[fname]:
        eq_name = equation.getElementsByTagName("name")[0].firstChild.data

        #f.write("  # " + eq_name + "\n")
        eq_string = equation.getElementsByTagName("eq")[0].firstChild.data
        eq_code = ""
        try:
          tokens = planktonica.parseString( eq_string )

          for t in tokens:
            eq_code = eq_code + eval_stmt(t)
        except ParseException, err:
          print "    Eq: " + eq_name + " ... Fail"
          print "    Parse Failure: ", eq_string
          print err
          raise
        except Exception, err:
          print "    Eq: " + eq_name + " ... Fail"
          print "    Evaluation failure, tokens: " + str(tokens)
          print err
          raise
        finally:
          file.write(Indent.line() + eq_code + "\n")

    # Add the housekeeping for pool updates
    file.write("\n  ### Setting pool variables\n")
    for v in self.pool_update_vars:
      eq_code = v.rhs + " = " + v.lhs
      file.write("  " + eq_code + "\n")

### Main model parsing ###
filename = sys.argv[1]
out_filename = filename.split(".")[0].strip() + '.py'
f = open(out_filename, "w")
f.write("import math\n")
f.write("import numpy\n")
f.write("from lebiology import stage_id, add_agent\n")

dom = parse(filename)
fgroup_doms = dom.getElementsByTagName("functionalgroup")
foodsets = dom.getElementsByTagName("foodsets")[0]

# Create Species...
species = []
for sdom in dom.getElementsByTagName("species")[1:]:
  species.append( Species(sdom) )

# Create FoodSets...
foodsets = []
foodsets_dom = dom.getElementsByTagName("foodsets")[0]
for fsdom in foodsets_dom.getElementsByTagName("foodset"):
  foodsets.append( FoodSet(fsdom) )

# Now create Functional Groups
for fg in fgroup_doms:
  fgroup = FGroup(fg)

  # Add species
  for s in species:
    if s.fg_name == fgroup.name:
      print "Adding species: " + s.name
      fgroup.species.append(s)


  # Add foodset
  for s in fgroup.species:
    for fs in foodsets:
      if s.name == fs.species_name:
        print "Adding foodset: " + fs.name
        fgroup.foodsets.append(fs)

  fgroup.write_parameters(f)

  sorted_stages = sorted(fgroup.stages.iteritems(), key=lambda stage: stage[1].id)
  for (sname, stage) in sorted_stages:
    fgroup.write_update_kernel(f, stage)


