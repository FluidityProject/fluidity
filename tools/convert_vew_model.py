#!/usr/bin/env python
import sys
import re
from pyparsing import Literal, Word, Group, OneOrMore, alphanums, nums, ParseException, Forward
from xml.dom.minidom import parse

### taken from: http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python
def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


### VEW grammar ###
ident = Word( alphanums + '_' + '$' )
number = Word( nums + '.' + '-' + 'e' + 'E' )
comma = Literal( ',').suppress()
lpar = Literal( '{' ).suppress()
rpar = Literal( '}' ).suppress()

assign = Literal( r'\assign' )
var = Literal( r'\var' )
val = Literal( r'\val' )
sival = Literal( r'\sival' )
unit = Literal( r'\unit' )
add = Literal( r'\add' )
sub = Literal( r'\sub' )
mul = Literal( r'\mul' )
div = Literal( r'\div' )
exp = Literal( r'\exp' )
pow = Literal( r'\pow' )
max = Literal( r'\max' )
min = Literal( r'\min' )
log10 = Literal( r'\log10' )
rnd = Literal( r'\rnd' )
cond = Literal( r'\conditional' )
ifthen = Literal( r'\ifthen' )
eq = Literal( r'\equal' )
geq = Literal( r'\greaterequal' )
gt = Literal( r'\greater' )
leq = Literal( r'\lessequal' )
lt = Literal( r'\less' )
boolor = Literal( r'\or' )
booland = Literal( r'\and' )
divide = Literal( r'\divide' )
uptake = Literal( r'\uptake' )
release = Literal( r'\release' )
change = Literal( r'\change' )
stage = Literal( r'\stage' )

expr = Forward()
bool_expr = Forward()
stmt = Forward()

variable = Group( var + lpar + ident + rpar )
si_value = sival + lpar + number + comma + number + rpar
unit_def = unit + lpar + number + comma + Word(alphanums) + comma + number + rpar
value = Group( val + lpar + si_value + comma + unit_def + rpar )
assignment = Group( assign + lpar + variable + comma + expr + rpar )

addition = Group( add + lpar + expr + OneOrMore( comma + expr ) + rpar )
subtraction = Group( sub + lpar + expr + comma + expr + rpar )
multiplication = Group( mul + lpar + expr + OneOrMore( comma + expr ) + rpar )
division = Group( div + lpar + expr + comma + expr + rpar )
exponential = Group( exp + lpar + expr + rpar )
power = Group( pow + lpar + expr + comma + expr + rpar )
maximum = Group( max + lpar + expr + comma + expr + rpar )
minimum = Group( min + lpar + expr + comma + expr + rpar )
logarithm = Group( log10 + lpar + expr + rpar )
random = Group( rnd + lpar + expr + rpar )

bool_op = ( geq | leq | eq | gt | lt )

bool_expr << ( Group( bool_op + lpar + expr + comma + expr + rpar ) 
             | Group( boolor + lpar + bool_expr + comma + bool_expr + rpar ) 
             | Group( booland + lpar + bool_expr + comma + bool_expr + rpar ) 
             )
conditional = Group( cond + lpar + bool_expr + comma + expr + comma + expr + rpar )
cond_stmt = Group( ifthen + lpar + bool_expr + comma + stmt + rpar )

cell_division = Group( divide + lpar + value + rpar )
chemical_uptake = Group( uptake + lpar + expr + comma + variable + rpar )
chemical_release = Group( release + lpar + expr + comma + variable + rpar )
change_stage =  Group( change + lpar + stage + lpar + ident + rpar + rpar )

expr << ( variable | value | addition | subtraction | multiplication | division | 
          exponential | conditional | power | maximum | minimum | random | logarithm )
stmt << ( assignment | cond_stmt | cell_division | chemical_uptake | chemical_release | change_stage )
planktonica = ( stmt )

indent = 2
id_counter = 0
stage_id = {}
stage_functions = {}
pool_update_vars = {}

functions_to_exclude = [ 'Motion' ]

def eval_token( t ):
  global indent, stage_id, pool_update_vars

  if t[0] == var:
    v = str(t[1])
    if v == 'TimeStep':
      v = 'dt_in_hours'
    if v == 'Temp':
      v = 'env["Temperature"]'
    if v == 'Vis_Irrad':
      v = 'env["Irradiance"]'
    if v.endswith('$Pool'):
      v = 'vars["' + v.split('$')[0] + '"]'
    if v.endswith('$Ingested'):
      v = 'vars["' + v.split('$')[0] + 'Ingested' + '"]'
    if v.endswith('$Conc'):
      v = 'env["Dissolved' + v.split('$')[0] + '"]'
    return v

  elif t[0] == assign:
    if t[1][0] == var:
      v = str(t[1][1])
      if v.endswith('$Pool'):
        # record the pool variable, so we can set it later
        v = v.split('$')[0]
        pool_update_vars[v + '_new'] = t[1]
        return v + '_new = ' + eval_token(t[2])
      else:
        return eval_token(t[1]) + ' = ' + eval_token(t[2])
    else:
      return eval_token(t[1]) + ' = ' + eval_token(t[2])

  elif t[0] == val:
    return str(float(t[2]))
  elif t[0] == add:
    c = '(' + eval_token(t[1])
    for summand in t[2:]:
      c = c + ' + ' + eval_token(summand)
    return c + ')'
  elif t[0] == sub:
    return '(' + eval_token(t[1]) + ' - ' + eval_token(t[2]) + ')'
  elif t[0] == mul:
    c = '(' + eval_token(t[1])
    for factor in t[2:]:
      c = c + ' * ' + eval_token(factor)
    return c + ')'
  elif t[0] == div:
    return '(' + eval_token(t[1]) + ' / ' + eval_token(t[2]) + ')'
  elif t[0] == exp:
    return 'math.exp(' + eval_token(t[1]) + ')'
  elif t[0] == log10:
    return 'math.log10(' + eval_token(t[1]) + ')'
  elif t[0] == pow:
    return 'math.pow(' + eval_token(t[1]) + ', ' + eval_token(t[2]) + ')'
  elif t[0] == max:
    return 'max(' + eval_token(t[1]) + ', ' + eval_token(t[2]) + ')'
  elif t[0] == min:
    return 'min(' + eval_token(t[1]) + ', ' + eval_token(t[2]) + ')'

  elif t[0] == cond:
    return '((' + eval_token(t[2]) + ') if (' + eval_token(t[1]) + ') else (' + eval_token(t[3]) + '))'
  elif t[0] == boolor:
    return '(' + eval_token(t[1]) + ') or (' + eval_token(t[2]) + ')'
  elif t[0] == booland:
    return '(' + eval_token(t[1]) + ') and (' + eval_token(t[2]) + ')'
  elif t[0] == eq:   # not sure why this is necessary here..?
    return '(' + eval_token(t[1]) + ' == ' + eval_token(t[2]) + ')'
  elif t[0] == gt:
    return '(' + eval_token(t[1]) + ' > ' + eval_token(t[2]) + ')'
  elif t[0] == geq:
    return '(' + eval_token(t[1]) + ' >= ' + eval_token(t[2]) + ')'
  elif t[0] == lt:
    return '(' + eval_token(t[1]) + ' < ' + eval_token(t[2]) + ')'
  elif t[0] == leq:
    return '(' + eval_token(t[1]) + ' <= ' + eval_token(t[2]) + ')'

  elif t[0] == ifthen:
    c = 'if ' + eval_token(t[1]) + ':\n'
    indent = indent + 2
    for i in range(indent):
      c = c + ' '
    c = c + eval_token(t[2])
    indent = indent - 2
    return c

  elif t[0] == divide:
    return 'vars["Size"] = vars["Size"] * ' + eval_token(t[1])
  elif t[0] == change:
    s = str(t[2])
    return 'vars["Stage"] = ' + str(float(stage_id[s])) + '  # ' + s
  elif t[0] == uptake:
    chem = str(t[2][1]).split('$')[0]
    return 'vars["' + chem + 'Uptake"] = ' + eval_token(t[1])
  elif t[0] == release:
    chem = str(t[2][1]).split('$')[0]
    return 'vars["' + chem + 'Release"] = ' + eval_token(t[1])
  else:
    raise Exception('Unkown token: ' + str(t))

### Main model parsing ###

filename = sys.argv[1]
out_filename = filename.split('.')[0].strip() + ".py"
f = open(out_filename, 'w')
f.write('import math\n')

dom = parse(filename)
fgroups = dom.getElementsByTagName('functionalgroup')[0:1]
for fg in fgroups:
  fg_name = fg.getElementsByTagName('name')[0].firstChild.data
  print "FG: " + fg_name

  params = {}
  for parameter in fg.getElementsByTagName('parameter'):
    param_name = parameter.getElementsByTagName('name')[0].firstChild.data
    val_str = parameter.getElementsByTagName('value')[0].firstChild.data
    params[param_name] =  float(val_str)

  stage_functions = {}
  stage_id = {}
  stage_functions = {}
  stages = fg.getElementsByTagName('stage')
  for stage in stages:
    stage_name = stage.getElementsByTagName('name')[0].firstChild.data
    stage_functions[stage_name] = []
    stage_id[stage_name] = id_counter
    id_counter = id_counter + 1

  func_equations = {}
  for function in fg.getElementsByTagName('function'):
    func_name = function.getElementsByTagName('name')[0].firstChild.data
    func_equations[func_name] = function.getElementsByTagName('equation')

    for call_stage in function.getElementsByTagName('calledin'):
      stage_functions[call_stage.firstChild.data].append(func_name)

  f.write('\n# Parameters for FG ' + fg_name + '\n')
  for p in sorted_nicely( params.keys() ):
    f.write(p + ' = ' + str(params[p]) + '\n')

  # Now print to file
  for stage in stage_functions.keys():
    print "Writing Stage: " + stage
    f.write('\n## FG: ' + fg_name + ';   Stage: ' + stage + ';   ID: ' + str(stage_id[stage]) )
    f.write('\ndef update_' + stage + '_' + fg_name + '(vars, env, dt):\n')
    f.write('  dt_in_hours = dt / 3600.0\n')

    pool_update_vars = {}
    for func in stage_functions[stage]:
      if func in functions_to_exclude:
        print "!!Excluding function: " + func
        continue

      print "  Function: " + func
      f.write('\n  ### ' + func + ' ###\n')
      for equation in func_equations[func]:
        eq_name = equation.getElementsByTagName('name')[0].firstChild.data
        f.write('  # ' + eq_name + '\n')
        eq_string = equation.getElementsByTagName('eq')[0].firstChild.data
        eq_code = ''

        try:
          tokens = planktonica.parseString( eq_string )

          for t in tokens:
            eq_code = eq_code + eval_token(t)
          print "    Eq: " + eq_name + " ... OK"
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
          f.write('  ' + eq_code + '\n')

    # Add the housekeeping for pool updates
    f.write('\n  ### Setting pool variables\n')
    for temp_var, token in pool_update_vars.iteritems():
      eq_code = eval_token(token) + ' = ' + temp_var
      f.write('  ' + eq_code + '\n')
