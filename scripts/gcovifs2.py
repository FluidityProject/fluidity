#!/usr/bin/env python

def number(s):
  res = True
  for i in range(len(s)):
    if not s[i].isdigit():
      res = False
      break
  return res

def elimIF(s):
  s=s.strip()
  s=s.lstrip('if')
  s=s.lstrip('IF')
  s=s.strip()
  s=s.lstrip('(')
  res=s.strip()
  return res

def elimTHEN(s):
  s=s.strip()
  s=s.rstrip('then')
  s=s.rstrip('THEN')
  s=s.rstrip()
  res = s[:len(s)-1]
  return res

def elimELSEIF(s):
  s=s.strip()
  s=s.lstrip('else')
  s=s.lstrip()
  s=s.lstrip('if')
  s=s.lstrip('ELSEIF')
  s=s.strip()
  s=s.lstrip('(')
  res=s.strip()
  return res

def elimBR(s):
  s=s.strip()
  s=s.rstrip('{')
  s=s.rstrip()
  res = s[:len(s)-1]
  return res
  

def gettokensF90(filename):
  f = open(filename,'r')
  line = f.readline()
  text = None
  while line != '':
    tokens=line.split()
    first = tokens[0]
    if tokens[0] == '#####:' and text != None:
      if text.find('have_option')!=-1:
        text = elimTHEN(text)
        print text,'\n'
        text = None
    elif len(tokens)>=3 and (tokens[2].startswith('if') or 
           tokens[2].startswith('IF') or 
           tokens[2].startswith('ELSEIF') or
           tokens[2].startswith('else if')):
      text=''
      cond = line.split(None,2)
      s = elimIF(cond[2])
      s = elimELSEIF(cond[2])
      text=filename+':'+text+cond[1]+s+'\n'
      if not number(first[:len(first)-1]):
        while line!='' and not tokens[len(tokens)-1].endswith('then') and not tokens[len(tokens)-1].endswith('THEN'):
          line = f.readline()
          cond = line.split(None,1)
          if len(cond)<2:
            text = text+filename+':'+line
          else:
            text=text+filename+':'+cond[1]
          tokens=line.split()
    elif number(first[:len(first)-1]):
      #either the else is not executed or one of the other cases in the if block
      text = None      
        
    
    line=f.readline()    
  f.close()

def gettokensCPP(filename):
  f = open(filename,'r')
  line = f.readline()
  text = None
  while line != '':
    tokens=line.split()
    first = tokens[0]
    if tokens[0] == '#####:' and text != None:
      if text.find('have_option')!=-1:
        text = elimBR(text)
        print text,'\n'
        text = None
    elif len(tokens)>=3 and (tokens[2].startswith('if') or 
           tokens[2].startswith('else if')):
      text=''
      cond = line.split(None,2)
      s = elimIF(cond[2])
      s = elimELSEIF(cond[2])
      text=filename+':'+text+cond[1]+s+'\n'
      if not number(first[:len(first)-1]):
        while line!='' and not tokens[len(tokens)-1].endswith('{') and not tokens[len(tokens)-1].endswith(';'):
          line = f.readline()
          cond = line.split(None,1)
          tok = cond[1].split()
          if tok[len(tok)-1].endswith('{'):
            break
          if len(cond)<2:
            text = text+filename+':'+line
          else:
            text=text+filename+':'+cond[1]
          tokens=line.split()
        else:
          if tokens[len(tokens)-1].endswith(';'):
            if text.find('have_option')!=-1:
              text = elimBR(text)
              print text,'\n'
            text=None
    elif number(first[:len(first)-1]):
      #either the else is not executed or one of the other cases in the if block
      text = None      
        
    
    line=f.readline()    
  f.close()


if __name__ == '__main__':
  import sys
  filen = sys.argv[1]
  textf = open(filen, 'r')
  l = textf.readline()
  while l != '':
    filen = l.strip();
    if filen.find('.F90.gcov')!=-1:
      gettokensF90(filen)
    elif filen.find('.cpp.gcov')!=-1:
      gettokensCPP(filen)
    l = textf.readline()
  textf.close()





