<?xml version="1.0" encoding="ISO-8859-1"?>
<!--
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
-->

<xsl:stylesheet version="1.0"
  xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:m="http://www.w3.org/1998/Math/MathML"
  xmlns:svg="http://www.w3.org/2000/svg">

<xsl:output indent="no" omit-xml-declaration="yes"/>
<xsl:strip-space elements="*"/>

<xsl:template match="/">
<span class="latexml">
  <xsl:apply-templates select=".//*[@class='document']"/>
</span>
</xsl:template>

<xsl:template match="div|p">
  <xsl:apply-templates/>
</xsl:template>

<xsl:template match="@* | node()">
  <xsl:copy>
    <xsl:apply-templates select="@* | node()"/>
  </xsl:copy>
</xsl:template>

</xsl:stylesheet>
