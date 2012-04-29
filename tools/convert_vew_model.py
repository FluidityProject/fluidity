#!/usr/bin/env python
import sys
import re
from pyparsing import Literal, Word, Group, OneOrMore, ZeroOrMore, alphanums, nums, ParseException, Forward, replaceWith
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
fgroup = None

### Metamodel definitions
class Stage:
  def __init__(self, dom_element):
    global id_counter
    self.dom = dom_element
    self.name = self.dom.getElementsByTagName("name")[0].firstChild.data
    self.id = id_counter
    self.functions = []
    self.function_eqns = {}
    self.assigned_locals = []
    id_counter = id_counter + 1

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

    self.stages = {}
    for s in self.dom.getElementsByTagName("stage"):
      sname = s.getElementsByTagName("name")[0].firstChild.data
      self.stages[sname] = Stage(s)

    for function in fg.getElementsByTagName("function"):
      fname = function.getElementsByTagName("name")[0].firstChild.data
      for stage in function.getElementsByTagName("calledin"):
        sname = stage.firstChild.data
        self.stages[sname].add_function(fname, function)

  def eval_var(self, t, hist_ind=None):
    v = str(t[1])
    if hist_ind != None:
      return "vars['" + v + "_" + str(hist_ind) + "']"

    if v == "TimeStep":
      v = "dt_in_hours"
    elif v == "PI":
      v = "math.pi"
    elif v == "Temp":
      v = "env['Temperature']"
    elif v == "Vis_Irrad":
      v = "env['Irradiance']"
    elif v == "Density":  # Use Pade approximation in Fluidity
      v = "(1000.*(env['Density'])-1000.)"
    elif v.endswith("$Pool"):
      v = "vars['" + v.split("$")[0] + "']"
    elif v.endswith("$Ingested"):
      v = "vars['" + v.split("$")[0] + "Ingested" + "']"
    elif v.endswith("$Conc"):
      v = "env['Dissolved" + v.split("$")[0] + "']"
    elif v == 'V_m':
      v = "vars['V_m']"

    elif v in self.parameters.keys():
      v = "param['" + v + "']"
    elif v in self.state_vars:
      v = "vars['" + v + "']"
    elif v in self.local_vars:
      v = v
      self.used_vars.append(v)
    elif v in self.variety_local:
      v = v
    elif v in self.variety_param.keys():
      v = "param['" + v + "']"
    # External values that we need to re-create VEW hacks
    elif v in ['MLDepth', 'Max_MLD', 'd_year']:
      v = "param['" + v + "']"
    # Food sets and ingested cells
    elif v in ['P']:
      v = "env['" + self.name + v + "Concentration']"
    elif v in ['IngestedCells']:
      v = "vars['P" + v + "']"

    # Dev exceptions...
    elif v in ['z', 'S_t']:
      v = "vars['" + v + "']"
    else:
      raise Exception("Unkown variable: " + v)
    return v

  ### The big eval function....
  def eval_token(self, t):
    global indent

    if t[0] == var:
      return self.eval_var(t)

    elif t[0] == assign:
      if t[1][0] == var:
        # When assigning state or pool variables, 
        # we create a buffer var and write the assignment at the end
        v = str(t[1][1])
        self.initialised_vars.append(v)
        self.used_vars = []
        if v.endswith("$Pool"):
          v = v.split("$")[0] + "_new"
          self.pool_update_vars[v] = t[1]
        elif v in self.state_vars:
          v = v + "_new"
          self.pool_update_vars[v] = t[1]
        elif v == 'V_m':
          v = "vars['V_m']"

        term = self.eval_token(t[2])
        code = ""
        for term_var in self.used_vars:
          if term_var not in self.initialised_vars:
            code = code + term_var + " = "+ str(self.local_vars[term_var]) + "\n"
            for i in range(indent):
              code = code + " "
            self.initialised_vars.append(term_var)
        return code + v + " = " + term
      else:
        raise Exception("Assigning a non-variable..?")

    elif t[0] == val:
      return str(float(t[2]))
    elif t[0] == add:
      c = "(" + self.eval_token(t[1])
      for summand in t[2:]:
        c = c + " + " + self.eval_token(summand)
      return c + ")"
    elif t[0] == sub:
      return "(" + self.eval_token(t[1]) + " - " + self.eval_token(t[2]) + ")"
    elif t[0] == mul:
      c = "(" + self.eval_token(t[1])
      for factor in t[2:]:
        c = c + " * " + self.eval_token(factor)
      return c + ")"
    elif t[0] == div:
      return "(" + self.eval_token(t[1]) + " / " + self.eval_token(t[2]) + ")"
    elif t[0] == exp:
      return "math.exp(" + self.eval_token(t[1]) + ")"
    elif t[0] == log10:
      return "math.log10(" + self.eval_token(t[1]) + ")"
    elif t[0] == pow:
      return "math.pow(" + self.eval_token(t[1]) + ", " + self.eval_token(t[2]) + ")"
    elif t[0] == max:
      return "max(" + self.eval_token(t[1]) + ", " + self.eval_token(t[2]) + ")"
    elif t[0] == min:
      return "min(" + self.eval_token(t[1]) + ", " + self.eval_token(t[2]) + ")"
    elif t[0] == minus:
      return "-" + self.eval_token(t[1])
    elif t[0] == abs:
      return "abs( " + self.eval_token(t[1]) + " )"

    elif t[0] == cond:
      return "((" + self.eval_token(t[2]) + ") if (" + self.eval_token(t[1]) + ") else (" + self.eval_token(t[3]) + "))"
    elif t[0] == boolor:
      return "(" + self.eval_token(t[1]) + ") or (" + self.eval_token(t[2]) + ")"
    elif t[0] == booland:
      return "(" + self.eval_token(t[1]) + ") and (" + self.eval_token(t[2]) + ")"
    elif t[0] == eq:
      return "(" + self.eval_token(t[1]) + " == " + self.eval_token(t[2]) + ")"
    elif t[0] == neq:
      return "(" + self.eval_token(t[1]) + " != " + self.eval_token(t[2]) + ")"
    elif t[0] == gt:
      return "(" + self.eval_token(t[1]) + " > " + self.eval_token(t[2]) + ")"
    elif t[0] == geq:
      return "(" + self.eval_token(t[1]) + " >= " + self.eval_token(t[2]) + ")"
    elif t[0] == lt:
      return "(" + self.eval_token(t[1]) + " < " + self.eval_token(t[2]) + ")"
    elif t[0] == leq:
      return "(" + self.eval_token(t[1]) + " <= " + self.eval_token(t[2]) + ")"

    elif t[0] == ifthen:
      c = "if " + self.eval_token(t[1]) + ":\n"
      indent = indent + 2
      for i in range(indent):
        c = c + " "
      c = c + self.eval_token(t[2])
      indent = indent - 2
      return c

    elif t[0] == divide:
      return "vars['Size'] = vars['Size'] * " + self.eval_token(t[1])
    elif t[0] == change:
      s = str(t[2])
      return "vars['Stage'] = stage_id('" + self.name + "', '" + s + "')"
    elif t[0] == uptake:
      chem = str(t[2][1]).split("$")[0]
      return "vars['" + chem + "Uptake'] = " + self.eval_token(t[1])
    elif t[0] == release:
      chem = str(t[2][1]).split("$")[0]
      return "vars['" + chem + "Release'] = " + self.eval_token(t[1])

    elif t[0] == varietysum:
      return self.eval_token(t[1])
    elif t[0] == varhist:
      hist_ind = int(float(self.eval_token(t[2])))
      return self.eval_var(t[1], hist_ind)
    elif t[0] == visIrradAt:
      return "param['surface_irradiance']"

    # TODO
    elif t[0] == integrate:
      return self.eval_token(t[1])
    elif t[0] == ingest:
      species_conc = self.eval_token(t[1])
      ing_threshold = self.eval_token(t[2])
      ing_amount = self.eval_token(t[3])
      # PRequest should not be hardcoded, but I'm lazy today...
      return "vars['PRequest'] = " + ing_amount + " if (" + species_conc + " > " + ing_threshold + ") else 0.0"
    elif t[0] == create:
      s = str(t[2])
      indent_str = ""
      for i in range(indent):
        indent_str = indent_str + " "
      c = "new_agent_vars = {}"
      c = c + "\n" + indent_str + "new_agent_vars['Stage'] = stage_id('" + self.name + "', '" + s + "')"
      c = c + "\n" + indent_str + "new_agent_vars['Size'] = " + self.eval_token(t[3])
      for i in range(4,len(t)):
        c = c + "\n" + indent_str + "new_agent_" + self.eval_token(t[i][1]) + " = " + self.eval_token(t[i][2]) 
      c = c + "\n" + indent_str + "lebiology.add_agent(new_agent_vars)"
      return c
    elif t[0] == pchange:
      s = str(t[2])
      return "#TODO PCHANGE( " + s + ", " + self.eval_token(t[3]) + " )" 

    elif t[0] == rnd:
      return "#TODO RND( " + self.eval_token(t[1]) + " )"

    else:
      raise Exception("Unkown token: " + str(t))
  ### End of big eval function

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
    if len(self.variety_param) > 0:
      file.write("### Variety parameters ###\n")
      for vp in sorted_nicely( self.variety_param.keys() ):
        file.write("    '" + vp + "' : " + str(self.variety_param[vp]) + ",\n")
    file.write("}\n")

  def write_update_kernel(self, file, stage):
    global fg_motion_functions

    print "  Writing Stage: " + stage.name
    file.write("\ndef update_" + stage.name + "_" + self.name + "(param, vars, env, dt):\n")
    file.write('  """ FGroup:  ' + self.name + '\n')
    file.write('      Stage:   ' + stage.name + '\n  """\n')
    file.write("  dt_in_hours = dt / 3600.0\n")

    self.pool_update_vars = {}
    self.initialised_vars = []    
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

          for t in tokens:
            eq_code = eq_code + self.eval_token(t)
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
    for temp_var, token in self.pool_update_vars.iteritems():
      eq_code = self.eval_token(token) + " = " + temp_var
      file.write("  " + eq_code + "\n")

### Main model parsing ###

fg_motion_functions = {
  "Diatom" : [ "Motion" ],
  "Copepod" : [ "Update depth" , "Pellet sinking" ], 
  "Predator" : [ "Pellet sinking" ],
  "Basal_predator" : [ "Pellet sinking" ]
}

fg_write_stages = {
  "Diatom" : [ "Living", "Dead" ],
  "Copepod" : [ "Dead", "OW5", "OWA5" ]
}

filename = sys.argv[1]
out_filename = filename.split(".")[0].strip() + '.py'
f = open(out_filename, "w")
f.write("import math\n")
f.write("from lebiology import stage_id\n")

dom = parse(filename)
fgroups = dom.getElementsByTagName("functionalgroup")[0:2]
species = dom.getElementsByTagName("species")[1:]

for fg in fgroups:
  fgroup = FGroup(fg)

  # Add and print species parameters
  for s in species:
    if s.getAttribute("fg") == fgroup.name: 
      sname = s.getAttribute("name").replace(' ', '_')
      print "Adding species: " + sname
      fgroup.add_species(sname , s )

  fgroup.write_parameters(f, sname)

  sorted_stages = sorted(fgroup.stages.iteritems(), key=lambda stage: stage[1].id)
  for (sname, stage) in sorted_stages:
    if sname in fg_write_stages[fgroup.name]:
      fgroup.write_update_kernel(f, stage)


