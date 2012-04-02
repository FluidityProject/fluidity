#!/usr/bin/env python
import sys
import re
from pyparsing import Literal, Word, Group, OneOrMore, ZeroOrMore, alphanums, nums, ParseException, Forward
from xml.dom.minidom import parse

### taken from: http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split("([0-9]+)", key) ] 
    return sorted(l, key = alphanum_key)


### VEW grammar ###
ident = Word( alphanums + "_" + "$" )
number = Word( nums + "." + "-" + "e" + "E" )
comma = Literal( ",").suppress()
lpar = Literal( "{" ).suppress()
rpar = Literal( "}" ).suppress()

assign = Literal( r"\assign" )
var = Literal( r"\var" )
val = Literal( r"\val" )
sival = Literal( r"\sival" )
unit = Literal( r"\unit" )
add = Literal( r"\add" )
sub = Literal( r"\sub" )
mul = Literal( r"\mul" )
div = Literal( r"\div" )
exp = Literal( r"\exp" )
pow = Literal( r"\pow" )
max = Literal( r"\max" )
min = Literal( r"\min" )
log10 = Literal( r"\log10" )
abs = Literal( r"\abs" )
minus = Literal( r"\minus" )
rnd = Literal( r"\rnd" )
cond = Literal( r"\conditional" )
ifthen = Literal( r"\ifthen" )
eq = Literal( r"\equal" )
neq = Literal( r"\neq" )
geq = Literal( r"\greaterequal" )
gt = Literal( r"\greater" )
leq = Literal( r"\lessequal" )
lt = Literal( r"\less" )
boolor = Literal( r"\or" )
booland = Literal( r"\and" )

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

expr = Forward()
bool_expr = Forward()
stmt = Forward()

variable = Group( var + lpar + ident + rpar )
si_value = sival + lpar + number + comma + number + rpar
unit_def = unit + lpar + number + comma + Word(alphanums) + comma + number + rpar
value = Group( val + lpar + si_value + comma + unit_def + rpar )
assignment = Group( assign + lpar + variable + comma + expr + rpar )

# Should be a OneOrMore, but there is a dummy addition in LERM-PS_1!!!
addition = Group( add + lpar + expr + ZeroOrMore( comma + expr ) + rpar )
subtraction = Group( sub + lpar + expr + comma + expr + rpar )
multiplication = Group( mul + lpar + expr + OneOrMore( comma + expr ) + rpar )
division = Group( div + lpar + expr + comma + expr + rpar )
exponential = Group( exp + lpar + expr + rpar )
power = Group( pow + lpar + expr + comma + expr + rpar )
maximum = Group( max + lpar + expr + comma + expr + rpar )
minimum = Group( min + lpar + expr + comma + expr + rpar )
logarithm = Group( log10 + lpar + expr + rpar )
negation = Group( minus + lpar + expr + rpar )
absolute = Group( abs + lpar + expr + rpar )
random = Group( rnd + lpar + expr + rpar )

bool_op = ( geq | leq | eq | gt | lt | neq )

bool_expr << ( Group( bool_op + lpar + expr + comma + expr + rpar ) 
             | Group( boolor + lpar + bool_expr + comma + bool_expr + rpar ) 
             | Group( booland + lpar + bool_expr + comma + bool_expr + rpar ) 
             )
conditional = Group( cond + lpar + bool_expr + comma + expr + comma + expr + rpar )
cond_stmt = Group( ifthen + lpar + bool_expr + comma + stmt + rpar )

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

expr << ( variable | value | addition | subtraction | multiplication | division | absolute | 
          exponential | conditional | power | maximum | minimum | random | logarithm | negation | 
          variety_vector_sum | variable_history | sample_irradiance | integration )
stmt << ( assignment | cond_stmt | cell_division | chemical_uptake | chemical_release | 
          change_stage | pchange_stage | ingestion | creation | set_variable )
planktonica = ( stmt )

indent = 2
id_counter = 0
stage_id = {}
pool_update_vars = {}
current_fgroup = None

def eval_token( t ):
  global indent, stage_id, pool_update_vars, current_fgroup

  if t[0] == var:
    v = str(t[1])
    if v == "TimeStep":
      v = "dt_in_hours"
    elif v == "PI":
      v = "math.pi"
    elif v == "Temp":
      v = "env['Temperature']"
    elif v == "Vis_Irrad":
      v = "env['Irradiance']"
    elif v == "Density":
      v = "env['Density']"
    elif v.endswith("$Pool"):
      v = "vars['" + v.split("$")[0] + "']"
    elif v.endswith("$Ingested"):
      v = "vars['" + v.split("$")[0] + "Ingested" + "']"
    elif v.endswith("$Conc"):
      v = "env['Dissolved" + v.split("$")[0] + "']"

    elif v in current_fgroup.parameters.keys():
      v = "param['" + v + "']"
    elif v in current_fgroup.state_vars:
      v = "vars['" + v + "']"
    elif v in current_fgroup.local_vars:
      v = v
    # TODO
    elif v in current_fgroup.variety_param:
      v = "VARIETYPARAM:" + v
    elif v == "I_gv":
      v = "#TODO SPECIAL:" + v

    # Dev exceptions...
    elif v in ['z', 'MLDepth', 'Max_MLD', 'P', 'I_gv', 'IngestedCells', 'd_year', 'S_t']:
      v = "SPECIAL:" + v
    else:
      raise Exception("Unkown variable: " + v)
    return v

  elif t[0] == assign:
    # When assigning state or pool variables, 
    # we create a buffer var and write the assignment at the end
    if t[1][0] == var:
      v = str(t[1][1])
      if v.endswith("$Pool"):
        v = v.split("$")[0]
        pool_update_vars[v + "_new"] = t[1]
        return v + "_new = " + eval_token(t[2])
      elif v in current_fgroup.state_vars:
        pool_update_vars[v + "_new"] = t[1]
        return v + "_new = " + eval_token(t[2])
      else:
        return eval_token(t[1]) + " = " + eval_token(t[2])
    else:
      return eval_token(t[1]) + " = " + eval_token(t[2])

  elif t[0] == val:
    return str(float(t[2]))
  elif t[0] == add:
    c = "(" + eval_token(t[1])
    for summand in t[2:]:
      c = c + " + " + eval_token(summand)
    return c + ")"
  elif t[0] == sub:
    return "(" + eval_token(t[1]) + " - " + eval_token(t[2]) + ")"
  elif t[0] == mul:
    c = "(" + eval_token(t[1])
    for factor in t[2:]:
      c = c + " * " + eval_token(factor)
    return c + ")"
  elif t[0] == div:
    return "(" + eval_token(t[1]) + " / " + eval_token(t[2]) + ")"
  elif t[0] == exp:
    return "math.exp(" + eval_token(t[1]) + ")"
  elif t[0] == log10:
    return "math.log10(" + eval_token(t[1]) + ")"
  elif t[0] == pow:
    return "math.pow(" + eval_token(t[1]) + ", " + eval_token(t[2]) + ")"
  elif t[0] == max:
    return "max(" + eval_token(t[1]) + ", " + eval_token(t[2]) + ")"
  elif t[0] == min:
    return "min(" + eval_token(t[1]) + ", " + eval_token(t[2]) + ")"
  elif t[0] == minus:
    return "-" + eval_token(t[1])
  elif t[0] == abs:
    return "abs( " + eval_token(t[1]) + " )"

  elif t[0] == cond:
    return "((" + eval_token(t[2]) + ") if (" + eval_token(t[1]) + ") else (" + eval_token(t[3]) + "))"
  elif t[0] == boolor:
    return "(" + eval_token(t[1]) + ") or (" + eval_token(t[2]) + ")"
  elif t[0] == booland:
    return "(" + eval_token(t[1]) + ") and (" + eval_token(t[2]) + ")"
  elif t[0] == eq:
    return "(" + eval_token(t[1]) + " == " + eval_token(t[2]) + ")"
  elif t[0] == neq:
    return "(" + eval_token(t[1]) + " != " + eval_token(t[2]) + ")"
  elif t[0] == gt:
    return "(" + eval_token(t[1]) + " > " + eval_token(t[2]) + ")"
  elif t[0] == geq:
    return "(" + eval_token(t[1]) + " >= " + eval_token(t[2]) + ")"
  elif t[0] == lt:
    return "(" + eval_token(t[1]) + " < " + eval_token(t[2]) + ")"
  elif t[0] == leq:
    return "(" + eval_token(t[1]) + " <= " + eval_token(t[2]) + ")"

  elif t[0] == ifthen:
    c = "if " + eval_token(t[1]) + ":\n"
    indent = indent + 2
    for i in range(indent):
      c = c + " "
    c = c + eval_token(t[2])
    indent = indent - 2
    return c

  elif t[0] == divide:
    return "vars['Size'] = vars['Size'] * " + eval_token(t[1])
  elif t[0] == change:
    s = str(t[2])
    return "vars['Stage'] = " + str(float(stage_id[s])) + "  # " + s
  elif t[0] == uptake:
    chem = str(t[2][1]).split("$")[0]
    return "vars['" + chem + "Uptake'] = " + eval_token(t[1])
  elif t[0] == release:
    chem = str(t[2][1]).split("$")[0]
    return "vars['" + chem + "Release'] = " + eval_token(t[1])

  # TODO
  elif t[0] == varietysum:
    return "TODO VARIETYSUM( " + eval_token(t[1]) + " )"
  elif t[0] == varhist:
    return "TODO VARHIST( " + eval_token(t[1]) + ", " + eval_token(t[2]) + " )"
  elif t[0] == visIrradAt:
    return "TODO VISUAL_IRRADIANCE_AT( " + eval_token(t[1]) + " )"
  elif t[0] == integrate:
    return "TODO INTEGRATE( " + eval_token(t[1]) + " )"
  elif t[0] == ingest:
    return "#TODO INGEST( " + eval_token(t[1]) + ", " + eval_token(t[2]) + ", " + eval_token(t[3]) + " )"
  elif t[0] == create:
    s = str(t[2])
    return "pass\n    #TODO CREATE( " + s + ", ...)!!! " 
  elif t[0] == pchange:
    s = str(t[2])
    return "#TODO PCHANGE( " + s + ", " + eval_token(t[3]) + " )" 

  elif t[0] == rnd:
    return "#TODO RND( " + eval_token(t[1]) + " )"

  else:
    raise Exception("Unkown token: " + str(t))

### Metamodel definitions
class Stage:
  def __init__(self, dom_element):
    global id_counter
    self.dom = dom_element
    self.name = self.dom.getElementsByTagName("name")[0].firstChild.data
    self.id = id_counter
    self.functions = []
    self.function_eqns = {}
    id_counter = id_counter + 1
    stage_id[self.name] = self.id # Global record for evaluating variables
    #print "  Stage: " + self.name + "    ID: " + str(self.id)

  def add_function(self, fname, function):
    self.functions.append(fname)
    self.function_eqns[fname] = function.getElementsByTagName("equation")

class FGroup:
  def __init__(self, dom_element):
    self.dom = dom_element
    self.name = self.dom.getElementsByTagName("name")[0].firstChild.data.replace(" ", "_")
    print "FGroup: " + self.name

    self.species = {}
    self.parameters = {}
    for parameter in self.dom.getElementsByTagName("parameter"):
      param_name = parameter.getElementsByTagName("name")[0].firstChild.data
      val_str = parameter.getElementsByTagName("value")[0].firstChild.data
      self.parameters[param_name] = float(val_str)

    self.state_vars = []
    for state_var in self.dom.getElementsByTagName("variable"):
      varname = state_var.getElementsByTagName("name")[0].firstChild.data
      self.state_vars.append(varname)

    self.local_vars = []
    for local in self.dom.getElementsByTagName("local"):
      varname = local.getElementsByTagName("name")[0].firstChild.data
      self.local_vars.append(varname)

    self.variety_param = []
    for vparam in self.dom.getElementsByTagName("varietyparameter"):
      varname = vparam.getElementsByTagName("name")[0].firstChild.data
      self.variety_param.append(varname)

    self.stages = {}
    for s in self.dom.getElementsByTagName("stage"):
      sname = s.getElementsByTagName("name")[0].firstChild.data
      self.stages[sname] = Stage(s)

    for function in fg.getElementsByTagName("function"):
      fname = function.getElementsByTagName("name")[0].firstChild.data
      for stage in function.getElementsByTagName("calledin"):
        sname = stage.firstChild.data
        self.stages[sname].add_function(fname, function)

  def add_species(self, name, species):
    self.species[name] = {}
    for parameter in species.getElementsByTagName("param"):
      param_name = parameter.getAttribute("name")
      val_str = parameter.getAttribute("a")
      self.species[name][param_name] = float(val_str)

  def write_parameters(self, file, species):
    file.write("\n# Parameters for FGroup " + self.name)
    file.write("\n# Species: " + species )
    file.write("\nparams_" + species + " = {\n")
    for p in sorted_nicely( self.species[species].keys() ):
      file.write("    '" + p + "' : " + str(self.species[species][p]) + ",\n")
    file.write("}\n")

  def write_update_kernels(self, file):
    global fg_motion_functions, pool_update_vars, current_fgroup
    sorted_stages = sorted(self.stages.iteritems(), key=lambda stage: stage[1].id)

    for (sname, stage) in sorted_stages:
      print "  Writing Stage: " + stage.name
      #file.write("\n## FG: " + fgroup.name + ";   Stage: " + stage.name + ";   ID: " + str(stage.id) )
      file.write("\ndef update_" + stage.name + "_" + fgroup.name + "(param, vars, env, dt):\n")
      file.write('  """ FGroup:  ' + self.name + '\n')
      file.write('      Stage:   ' + stage.name + '   ID: ' + str(stage.id) + '\n  """\n')
      file.write("  dt_in_hours = dt / 3600.0\n")

      pool_update_vars = {}
      for fname in stage.functions:
        if fname in fg_motion_functions[self.name]:
          continue

        #print "  Function: " + fname
        file.write("\n  ### " + fname + " ###\n")
        for equation in stage.function_eqns[fname]:
          eq_name = equation.getElementsByTagName("name")[0].firstChild.data

          #f.write("  # " + eq_name + "\n")
          eq_string = equation.getElementsByTagName("eq")[0].firstChild.data
          eq_code = ""
          try:
            tokens = planktonica.parseString( eq_string )
            current_fgroup = self

            for t in tokens:
              eq_code = eq_code + eval_token(t)
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
            file.write("  " + eq_code + "\n")

      # Add the housekeeping for pool updates
      file.write("\n  ### Setting pool variables\n")
      for temp_var, token in pool_update_vars.iteritems():
        eq_code = eval_token(token) + " = " + temp_var
        file.write("  " + eq_code + "\n")

### Main model parsing ###

fg_motion_functions = {
  "Diatom" : [ "Motion" ],
  "Copepod" : [ "Update depth" , "Pellet sinking" ], 
  "Predator" : [ "Pellet sinking" ],
  "Basal_predator" : [ "Pellet sinking" ]
}

filename = sys.argv[1]
out_filename = filename.split(".")[0].strip() + '.py'
f = open(out_filename, "w")
f.write("import math\n")

dom = parse(filename)
fgroups = dom.getElementsByTagName("functionalgroup")[0:1]
species = dom.getElementsByTagName("species")[1:]

for fg in fgroups:
  fgroup = FGroup(fg)
  for s in species:
    #print s.toxml()
    if s.getAttribute("fg") == fgroup.name: 
      sname = s.getAttribute("name").replace(' ', '_')
      print "Adding species: " + sname
      fgroup.add_species(sname , s )
      fgroup.write_parameters(f, sname)
  fgroup.write_update_kernels(f)


