#!/usr/bin/python
import textwrap; 
s=open('input_large.pgm', 'rb').read().encode("hex").upper(); 
t="".join(["\\x"+x+y for (x,y) in zip(s[0::2], s[1::2])]) ; 
print "static const char *input = \\\n\t\"%s\";"%"\" \\\n\t\"".join(textwrap.wrap(t,80))
