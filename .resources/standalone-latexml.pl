#!/usr/bin/perl
# ProjectEngine compute server
# Copyright (C) 2010 Lee Worden <worden.lee@gmail.com>
# http://lalashan.mcmaster.ca/theobio/projects/ProjectEngine
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

#my($WAITING,$BODYTAG,$INBODY,$PASTBODY) = (1 .. 4);
#my($state) = $WAITING;
my($dbg) = 0;

my($latexml_css) = my($base) = $ENV{'PUBLISH_URL'};
$latexml_css =~ s/project=.*?&amp;//;
$latexml_css .= "css\/latexml.css";

while(<>)
{ # this is shameful - learn to use XSLT right!
  s/(<|<\/)m:/$1/g;
  s|<math\s|<math xmlns="http://www.w3.org/1998/Math/MathML" |ig;
  s/(<|<\/)svg:/$1/g;
  s|<svg\s|<svg xmlns="http://www.w3.org/2000/svg" |ig;
  my($add_to_head) = '
    <script type="text/javascript" '
             .'src="http://cdn.mathjax.org/mathjax/latest/MathJax.js'
                      .'?config=TeX-AMS-MML_HTMLorMML"></script>';
  #<base href=\"http://lalashan.mcmaster.ca$publish_url\"/>";
  s/<head>/<head>$add_to_head/i;
  s/core.css/$latexml_css/;
  s/<body/<body class="latexml"/i;
  s/src="([^\/"]+(\/[^\/"]+)*)"/src=\"$base$1\"/i;
  print;
}
exit;
