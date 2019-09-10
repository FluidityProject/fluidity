### Copyright (C) 2005 Thomas M. Hinkle
### Copyright (C) 2007 Imperial College London and others.
### Please see the AUTHORS file for a full list of contributors.
###
### This library is free software; you can redistribute it and/or
### modify it under the terms of the GNU General Public License as
### published by the Free Software Foundation; either version 2 of the
### License, or (at your option) any later version.
###
### This library is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
### General Public License for more details.
###
### You should have received a copy of the GNU General Public License
### along with this library; if not, write to the Free Software
### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
### 02111-1307, USA.

from gi.repository import Pango as pango
from gi.repository import Gtk as gtk
from gi.repository import Gdk as gdk
import xml.sax.saxutils

from . import debug

class AttrIterator():
    def __init__(self, attributes=[]):
        self.attributes = attributes
        self.attribute_stack = []
        self.start_index = 0
        self.end_index = 0
        if not self.next():
            self.end_index = 2 ** 32 - 1
  
    def next(self):
        if len(self.attributes) == 0 and len(self.attribute_stack) == 0:
            return False
        self.start_index = self.end_index
        self.end_index = 2 ** 32 - 1
  
        to_remove = []
        for attr in self.attribute_stack:
            if attr.end_index == self.start_index:
                to_remove.append(attr)
            else:
                self.end_index = min(self.end_index, attr.end_index)
  
        while len(to_remove) > 0:
            attr = to_remove[0]
            self.attribute_stack.remove(to_remove[0])
            try:
                to_remove.remove(attr)
            except:
                pass
  
        while len(self.attributes) != 0 and self.attributes[0].start_index == self.start_index:
            if self.attributes[0].end_index > self.start_index:
                self.attribute_stack.append(self.attributes[0])
                self.end_index = min(self.end_index, self.attributes[0].end_index)
            self.attributes = self.attributes[1:]
        if len(self.attributes) > 0:
            self.end_index = min(self.end_index, self.attributes[0].start_index)
        return True
  
    def range(self):
        return (self.start_index, self.end_index)
  
    def get_font(self):
        tmp_list1 = self.attribute_stack
        fontdesc = pango.FontDescription()
        for attr in self.attribute_stack:
            if attr.klass.type == pango.AttrType.FONT_DESC:
                tmp_list1.remove(attr)
                attr.__class__ = pango.AttrFontDesc
                fontdesc = attr.desc
        return (fontdesc, None, self.attribute_stack)
  
  
def get_iterator(self):
    tmplist = []
  
    def fil(val, unused_data):
        tmplist.append(val)
        return False
  
    self.filter(fil, None)
    return AttrIterator(tmplist)
  
  
setattr(pango.AttrList, 'get_iterator', get_iterator)

class PangoBuffer (gtk.TextBuffer):
    desc_to_attr_table = {
        'family':[pango.AttrType.FAMILY,""],
        'style':[pango.AttrType.STYLE,pango.Style.NORMAL],
        'variant':[pango.AttrType.VARIANT,pango.Style.NORMAL],
        'weight':[pango.AttrType.WEIGHT,pango.Weight.NORMAL],
        'stretch':[pango.AttrType.STRETCH,pango.Stretch.NORMAL],
        }
    desc_to_attr_type = {
        'family':pango.AttrString,
        'size':pango.AttrInt,
        'weight':pango.AttrInt,
        'underline':pango.AttrInt,
        'stretch':pango.AttrInt,
        'variant':pango.AttrInt,
        'style':pango.AttrInt,
        'scale':pango.AttrFloat,
        'strikethrough':pango.AttrInt,
        'rise':pango.AttrInt,
        }
    attval_to_markup={
            'underline':{pango.Underline.SINGLE:'single',
                         pango.Underline.DOUBLE:'double',
                         pango.Underline.LOW:'low',
                         pango.Underline.NONE:'none'},
            'stretch':{pango.Stretch.ULTRA_EXPANDED:'ultraexpanded',
                       pango.Stretch.EXPANDED:'expanded',
                       pango.Stretch.EXTRA_EXPANDED:'extraexpanded',
                       pango.Stretch.EXTRA_CONDENSED:'extracondensed',
                       pango.Stretch.ULTRA_CONDENSED:'ultracondensed',
                       pango.Stretch.CONDENSED:'condensed',
                       pango.Stretch.NORMAL:'normal',
                       },
            'variant':{pango.Variant.NORMAL:'normal',
                       pango.Variant.SMALL_CAPS:'smallcaps',
                       },
            'style':{pango.Style.NORMAL:'normal',
                     pango.Style.OBLIQUE:'oblique',
                     pango.Style.ITALIC:'italic',
                     },
            'stikethrough':{1:'true',
                            True:'true',
                            0:'false',
                            False:'false'},
            }
    def __init__ (self):
        self.tagdict = {}
        self.tags = {}
        #self.buf = buf
        #self.set_text(txt)
        gtk.TextBuffer.__init__(self)

    def set_text (self, txt):
        gtk.TextBuffer.set_text(self,"")
        try:
            success,self.parsed,self.txt,self.separator = pango.parse_markup(txt,-1,'\x00')
        except:
            debug.deprint('Escaping text, we seem to have a problem here!', 2)
            txt=xml.sax.saxutils.escape(txt)
            success,self.parsed,self.txt,self.separator = pango.parse_markup(txt,-1,'\x00')
        self.attrIter = self.parsed.get_iterator()
        self.add_iter_to_buffer()
        while self.attrIter.next():
            self.add_iter_to_buffer()

    def add_iter_to_buffer (self):
        range=self.attrIter.range()
        font,lang,attrs = self.attrIter.get_font()
        tags = self.get_tags_from_attrs(font,lang,attrs)
        text = self.txt[range[0]:range[1]]
        if tags: self.insert_with_tags(self.get_end_iter(),text,*tags)
        else: self.insert_with_tags(self.get_end_iter(),text)

    def get_tags_from_attrs (self, font,lang,attrs):
        tags = []
        if font:
            font,fontattrs = self.fontdesc_to_attrs(font)
            fontdesc = font.to_string()
            if fontattrs:
                attrs.extend(fontattrs)
            if fontdesc and fontdesc!='Normal':
                if font.to_string() not in self.tags:
                    tag=self.create_tag()
                    tag.set_property('font-desc',font)
                    if tag not in self.tagdict: self.tagdict[tag]={}
                    self.tagdict[tag]['font_desc']=font.to_string()
                    self.tags[font.to_string()]=tag
                tags.append(self.tags[font.to_string()])
        if lang:
            if lang not in self.tags:
                tag = self.create_tag()
                tag.set_property('language',lang)
                self.tags[lang]=tag
            tags.append(self.tags[lang])
        if attrs:
            for a in attrs:
                if a.klass.type == pango.AttrType.FOREGROUND:
                    a.__class__ = pango.AttrColor
                    gdkcolor = self.pango_color_to_gdk(a.color)
                    key = 'foreground%s'%self.color_to_hex(gdkcolor)
                    if key not in self.tags:
                        self.tags[key]=self.create_tag()
                        self.tags[key].set_property('foreground-gdk',gdkcolor)
                        self.tagdict[self.tags[key]]={}
                        self.tagdict[self.tags[key]]['foreground']="#%s"%self.color_to_hex(gdkcolor)
                    tags.append(self.tags[key])
                    a.__class__ = pango.Attribute
                if a.klass.type == pango.AttrType.BACKGROUND:
                    a.__class__ = pango.AttrColor
                    gdkcolor = self.pango_color_to_gdk(a.color)
                    key = 'background%s'%self.color_to_hex(gdkcolor)
                    if key not in self.tags:
                        self.tags[key]=self.create_tag()
                        self.tags[key].set_property('background-gdk',gdkcolor)
                        self.tagdict[self.tags[key]]={}
                        self.tagdict[self.tags[key]]['background']="#%s"%self.color_to_hex(gdkcolor)
                    tags.append(self.tags[key])
                    a.__class__ = pango.Attribute
                if a.klass.type.value_nick in self.desc_to_attr_type:
                    prop=a.klass.type.value_nick
                    a.__class__ = self.desc_to_attr_type[prop]
                    val=getattr(a,'value')
                    mval = val
                    if prop in self.attval_to_markup:
                        if val in self.attval_to_markup[prop]:
                            mval = self.attval_to_markup[prop][val]
                        else:
                            debug.deprint("hmmm, didn't know what to do with value %s"%val, 2)
                    key="%s%s"%(prop,val)
                    if key not in self.tags:
                        self.tags[key]=self.create_tag()
                        self.tags[key].set_property(prop,val)
                        self.tagdict[self.tags[key]]={}
                        self.tagdict[self.tags[key]][prop]=mval
                    tags.append(self.tags[key])
                    a.__class__ = pango.Attribute
        return tags

    def get_tags (self):
        tagdict = {}
        for pos in range(self.get_char_count()):
            iter=self.get_iter_at_offset(pos)
            for tag in iter.get_tags():
                if tag in tagdict:
                    if tagdict[tag][-1][1] == pos - 1:
                        tagdict[tag][-1] = (tagdict[tag][-1][0],pos)
                    else:
                        tagdict[tag].append((pos,pos))
                else:
                    tagdict[tag]=[(pos,pos)]
        return tagdict

    def get_text (self, start=None, end=None, include_hidden_chars=True):
        tagdict=self.get_tags()
        if not start: start=self.get_start_iter()
        if not end: end=self.get_end_iter()
        txt = str(gtk.TextBuffer.get_text(self,start,end))
        cuts = {}
        for k,v in list(tagdict.items()):
            stag,etag = self.tag_to_markup(k)
            for st,e in v:
                if st in cuts: cuts[st].append(stag) #add start tags second
                else: cuts[st]=[stag]
                if e+1 in cuts: cuts[e+1]=[etag]+cuts[e+1] #add end tags first
                else: cuts[e+1]=[etag]
        last_pos = 0
        outbuff = ""
        cut_indices = list(cuts.keys())
        cut_indices.sort()
        soffset = start.get_offset()
        eoffset = end.get_offset()
        cut_indices = [i for i in cut_indices if eoffset >= i >= soffset]
        for c in cut_indices:
            if not last_pos==c:
                outbuff += xml.sax.saxutils.escape(txt[last_pos:c])
                last_pos = c
            for tag in cuts[c]:
                outbuff += tag
        outbuff += xml.sax.saxutils.escape(txt[last_pos:])
        return outbuff

    def tag_to_markup (self, tag):
        stag = "<span"
        for k,v in list(self.tagdict[tag].items()):
            stag += ' %s="%s"'%(k,v)
        stag += ">"
        return stag,"</span>"

    def fontdesc_to_attrs (self,font):
        nicks = font.get_set_fields().value_nicks
        attrs = []
        for n in nicks:
            if n in self.desc_to_attr_table:
                attrType, norm = self.desc_to_attr_table[n]
                # create an attribute with our current value
                a = pango.Attribute()
                a.init(pango.AttrClass())
                a.klass.type = attrType
                a.__class__ = self.desc_to_attr_type[n]
                a.value = getattr(font,'get_%s'%n)()
                a.__class__ = pango.Attribute
                attrs.append(a)
                # unset our font's value
                getattr(font,'set_%s'%n)(norm)
        return font,attrs

    def pango_color_to_gdk (self, pc):
        return gdk.Color(pc.red,pc.green,pc.blue)

    def color_to_hex (self, color):
        hexstring = ""
        for col in 'red','green','blue':
            hexfrag = hex(getattr(color,col)//(16*16)).split("x")[1]
            if len(hexfrag)<2: hexfrag = "0" + hexfrag
            hexstring += hexfrag
        return hexstring

    def apply_font_and_attrs (self, font, attrs):
        tags = self.get_tags_from_attrs(font,None,attrs)
        for t in tags: self.apply_tag_to_selection(t)

    def remove_font_and_attrs (self, font, attrs):
        tags = self.get_tags_from_attrs(font,None,attrs)
        for t in tags: self.remove_tag_from_selection(t)

    def setup_default_tags (self):
        self.italics = self.get_tags_from_attrs(None,None,[pango.AttrStyle('italic')])[0]
        self.bold = self.get_tags_from_attrs(None,None,[pango.AttrWeight('bold')])[0]
    def get_selection (self):
        bounds = self.get_selection_bounds()
        if not bounds:
            iter=self.get_iter_at_mark(self.get_insert())
            if iter.inside_word():
                start_pos = iter.get_offset()
                iter.forward_word_end()
                word_end = iter.get_offset()
                iter.backward_word_start()
                word_start = iter.get_offset()
                iter.set_offset(start_pos)
                bounds = (self.get_iter_at_offset(word_start),
                          self.get_iter_at_offset(word_end+1))
            else:
                bounds = (iter,self.get_iter_at_offset(iter.get_offset()+1))
        return bounds

    def apply_tag_to_selection (self, tag):
        selection = self.get_selection()
        if selection:
            self.apply_tag(tag,*selection)

    def remove_tag_from_selection (self, tag):
        selection = self.get_selection()
        if selection:
            self.remove_tag(tag,*selection)

    def remove_all_tags (self):
        selection = self.get_selection()
        if selection:
            for t in list(self.tags.values()):
                self.remove_tag(t,*selection)

class InteractivePangoBuffer (PangoBuffer):
    def __init__ (self,
                  normal_button=None,
                  toggle_widget_alist=[]):
        """An interactive interface to allow marking up a gtk.TextBuffer.
        txt is initial text, with markup.
        buf is the gtk.TextBuffer
        normal_button is a widget whose clicked signal will make us normal
        toggle_widget_alist is a list that looks like this:
        [(widget, (font,attr)),
         (widget2, (font,attr))]
         """
        PangoBuffer.__init__(self)
        if normal_button: normal_button.connect('clicked',lambda *args: self.remove_all_tags())
        self.tag_widgets = {}
        self.internal_toggle = False
        self.insert = self.get_insert()
        self.connect('mark-set',self._mark_set_cb)
        self.connect('changed',self._changed_cb)
        for w,tup in toggle_widget_alist:
            self.setup_widget(w,*tup)

    def setup_widget_from_pango (self, widg, markupstring):
        """setup widget from a pango markup string"""
        #font = pango.FontDescription(fontstring)
        success,a,t,s = pango.parse_markup(markupstring,-1,'\x00')
        ai=a.get_iterator()
        font,lang,attrs=ai.get_font()
        return self.setup_widget(widg,font,attrs)

    def setup_widget (self, widg, font, attr):
        tags=self.get_tags_from_attrs(font,None,attr)
        self.tag_widgets[tuple(tags)]=widg
        return widg.connect('toggled',self._toggle,tags)

    def _toggle (self, widget, tags):
        if self.internal_toggle: return
        if widget.get_active():
            for t in tags: self.apply_tag_to_selection(t)
        else:
            for t in tags: self.remove_tag_from_selection(t)

    def _mark_set_cb (self, buffer, iter, mark, *params):
        # Every time the cursor moves, update our widgets that reflect
        # the state of the text.
        if hasattr(self,'_in_mark_set') and self._in_mark_set: return
        self._in_mark_set = True
        if mark.get_name()=='insert':
            for tags,widg in list(self.tag_widgets.items()):
                active = True
                for t in tags:
                    if not iter.has_tag(t):
                        active=False
                self.internal_toggle=True
                widg.set_active(active)
                self.internal_toggle=False
        if hasattr(self,'last_mark'):
            self.move_mark(self.last_mark,iter)
        else:
            self.last_mark = self.create_mark('last',iter,left_gravity=True)
        self._in_mark_set = False

    def _changed_cb (self, tb):
        if not hasattr(self,'last_mark'): return
        # If our insertion point has a mark, we want to apply the tag
        # each time the user types...
        old_itr = self.get_iter_at_mark(self.last_mark)
        insert_itr = self.get_iter_at_mark(self.insert)
        if old_itr!=insert_itr:
            # Use the state of our widgets to determine what
            # properties to apply...
            for tags,w in list(self.tag_widgets.items()):
                if w.get_active():
                    for t in tags: self.apply_tag(t,old_itr,insert_itr)





class SimpleEditor:
    def __init__ (self):
        self.w = gtk.Window()
        self.vb = gtk.VBox()
        self.editBox = gtk.HButtonBox()
        self.nb = gtk.Button('Normal')
        self.editBox.add(self.nb)
        self.sw = gtk.ScrolledWindow()
        self.tv = gtk.TextView()
        self.sw.add(self.tv)
        self.ipb = InteractivePangoBuffer(
            normal_button=self.nb)
        self.ipb.set_text("""<b>This is bold</b>. <i>This is italic</i>
            <b><i>This is bold, italic, and <u>underlined!</u></i></b>
            <span background="blue">This is a test of bg color</span>
            <span foreground="blue">This is a test of fg color</span>
            <span foreground="white" background="blue">This is a test of fg and bg color</span>
            """)
        #    Here are some more: 1-2, 2-3, 3-4, 10-20, 30-40, 50-60
        #    This is <span color="blue">blue</span>, <span color="red">red</span> and <span color="green">green</span>""")
        #self.ipb.set_text("""This is a numerical range (three hundred and fifty to four hundred) 350-400 which may get messed up.
        #Here are some more: 1-2, 2-3, 3-4, 10-20, 30-40, 50-60""")

        self.tv.set_buffer(self.ipb)
        for lab,stock,font in [('gtk-italic',True,'<i>italic</i>'),
                               ('gtk-bold',True,'<b>bold</b>'),
                               ('gtk-underline',True,'<u>underline</u>'),
                               ('Blue',True,'<span foreground="blue">blue</span>'),
                               ('Red',False,'<span foreground="red">smallcaps</span>'),
                               ]:
            button = gtk.ToggleButton(lab)
            self.editBox.add(button)
            if stock: button.set_use_stock(True)
            self.ipb.setup_widget_from_pango(button,font)
        self.vb.add(self.editBox)
        self.vb.add(self.sw)
        self.actionBox = gtk.HButtonBox()
        self.qb = gtk.Button(stock='quit')
        self.pmbut = gtk.Button('Print markup')
        self.pmbut.connect('clicked',self.print_markup)
        self.qb.connect('clicked',lambda *args: self.w.destroy() or gtk.main_quit())
        self.actionBox.add(self.pmbut)
        self.actionBox.add(self.qb)
        self.vb.add(self.actionBox)
        self.w.add(self.vb)
        self.w.show_all()

    def print_markup (self,*args):
        debug.dprint(self.ipb.get_text(), 0)
