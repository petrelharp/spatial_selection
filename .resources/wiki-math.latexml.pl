#!/usr/bin/perl
# ProjectEngine compute server
# Copyright (C) 2010 Lee Worden <worden.lee@gmail.com>
# http://lalashan.mcmaster.ca/theobio/projects/index.php/ProjectEngine
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
# http://www.gnu.org/copyleft/gpl.html

# this is called after using latexmlmath to render an inline latex expression.
# it expects a single <math> element, which it wraps in
# <span class="latexml-document">.

my($dbg) = 0;

print "<span class=\"latexml\">";
while(<>)
{ chomp;
  next if /^<\?/; # skip the opening XML line
  s/^\s+//;
  # this is shameful - learn to use XSLT right!
  s/(<|<\/)m:/$1/g;
  s|<math\s|<math xmlns="http://www.w3.org/1998/Math/MathML" |ig;
  s/display="block"/display="inline"/ig;
  s/(<|<\/)svg:/$1/g;
  s|<svg\s|<svg xmlns="http://www.w3.org/2000/svg" |ig;
  print;
}
print "</span>";
exit;
