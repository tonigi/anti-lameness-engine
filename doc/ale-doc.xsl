<!--  
  -  Copyright 2006 David Hilvert
  - 
  -  ALE documentation stylesheet definition.
  -
  -  This file is part of the Anti-Lamenessing Engine documentation.
  -
  -  The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
  -  it under the terms of the GNU General Public License as published by
  -  the Free Software Foundation; either version 2 of the License, or
  -  (at your option) any later version.
  -
  -  The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
  -  but WITHOUT ANY WARRANTY; without even the implied warranty of
  -  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  -  GNU General Public License for more details.
  -
  -  You should have received a copy of the GNU General Public License
  -  along with Anti-Lamenessing Engine; if not, write to the Free Software
  -  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  -->

<!DOCTYPE xsl:stylesheet [

	<!--
	  -  Entities
	  -
	  -
	  -  NOTE: all of the below entities operate on edit records, which
	  -  are expected to conform to exactly one of the following patterns:
	  -
	  -  <edit by="GIVEN-NAME SURNAME" on="2006-Sep-24">
	  -  <edit by="GIVEN-NAME SURNAME" in-month="2006-Sep">
	  -  <edit by="GIVEN-NAME SURNAME" in-year="2006">
	  -
	  -  The first form is preferred for occasional edits; the latter two
	  -  are provided as shorthand for edits occurring over a number of
	  -  days in succesion, or for cases where the exact date of the edit
	  -  is not known (e.g., when basing edit history on copyright notices
	  -  from older files).
	  -->


	<!--
	  -  Obtain an editor name from an edit record.
	  -->

	<!ENTITY editor 'normalize-space(@by)'>

	<!--
	  -  Obtain an edit year from an edit record.
	  -->

	<!ENTITY year 'concat(@in-year, substring-before(@in-month, "-"), substring-before(@on, "-"))'>

	<!--
	  -  Generate a string unique for each (editor, year) combination.
	  -->

	<!ENTITY editor-year 'concat(&editor;, &year;)'>

	<!--
	  -  Boolean test to determine whether the current node is within the scope of
	  -  $objectroot.
	  -->

	<!ENTITY scope 'count(ancestor::node()|$objectroot) = count(ancestor::node())'>

	<!--
	  -  Editor first name.
	  -->

	<!ENTITY editor-firstname 'substring-before(&editor;, " ")'>

	<!--
	  -  Editor last name
	  -->

	<!ENTITY editor-surname 'substring-after(&editor;, " ")'>

	<!--
	  -  Sort order for outputs including editor names.
	  -->

	<!ENTITY editor-sort-order '&editor-surname;'>

]>

<xsl:stylesheet id="style1"
                version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
		xmlns:xi="http://www.w3.org/2001/XInclude">

	<!-- 
	  -  Product information
          -->

	<xsl:param name="product-name" select="'ALE'"/>

	<xsl:param name="product-number" select="'0.8.5-prerelease'"/>

	<xsl:param name="copyright-info">
	  <copyright id="copyright-info">
	    <year> 2002, 2003, 2004, 2005, 2006 </year>
	    <holder> David Hilvert </holder>
	  </copyright>
	</xsl:param>

	<!-- 
	  -  Editing information
          -->

	<xsl:key name="editors" match="edit" use="&editor;"/>
	<xsl:key name="editor-years" match="edit" use="&editor-year;"/>

	<!-- 
	  -  If no other template matches, we copy the tree structure unchanged.
          -->

	<xsl:template match="@*|node()">
	  <xsl:copy>
	      <xsl:apply-templates select="@*|node()"/>
          </xsl:copy>
	</xsl:template>

	<!-- 
	  -  Do not pass edit tags to DocBook, as it doesn't understand them.
          -->

	<xsl:template match="edit"/>

	<!-- 
	  -  Add package information to titles of articles, books and sets.
          -->

	<xsl:template match="setinfo/title|bookinfo/title|articleinfo/title">
	  <title>
		  <xsl:copy-of select="$product-name"/>
		  <xsl:text> </xsl:text>
		  <xsl:copy-of select="$product-number"/>
		  <xsl:text> </xsl:text>
	          <xsl:apply-templates/>
	  </title>
	</xsl:template>

	<!-- 
	  -  Generate title, author, and copyright information for articles,
	  -  books and sets.
          -->

	<xsl:template match="setinfo|bookinfo|articleinfo">
	<xsl:copy>
	  <xsl:param name="objectroot" select=".."/>
	  <xsl:param name="editors-unique" select="$objectroot//edit[count(.|key('editors', &editor;)[&scope;][1]) = 1]"/>
	  <xsl:param name="editor-years-unique" select="$objectroot//edit[count(.|key('editor-years', &editor-year;)[&scope;][1]) = 1]"/>

  	  <!-- 
	    -  Generate the title, if available.
            -->

	    <xsl:apply-templates select="title"/>

	  <!--
	    -  Preserve abstracts.
	    -->

	    <xsl:apply-templates select="abstract"/>

  	  <!-- 
	    -  Add author information.
            -->

	    <xsl:for-each select="$editors-unique">
	      <xsl:sort select="&editor-sort-order;"/>
	      <author>
		<firstname>
		  <xsl:value-of select="&editor-firstname;"/>
		</firstname>
		<surname>
		  <xsl:value-of select="&editor-surname;"/>
		</surname>
	      </author>
	    </xsl:for-each>

  	  <!-- 
	    -  Add copyright information
            -->

	    <xsl:for-each select="$editors-unique">
	      <xsl:sort select="&editor-sort-order;"/>
	      <xsl:param name="this-editor" select="&editor;"/>
	      <copyright>
	      <holder>
	        <xsl:value-of select="&editor;"/>
	      </holder>
	      <xsl:for-each select="$editor-years-unique">
	        <xsl:if test="&editor; = $this-editor">
		  <year>
		    <xsl:value-of select="&year;"/>
		  </year>
		</xsl:if>
	      </xsl:for-each>
	      </copyright>
	    </xsl:for-each>

	</xsl:copy>
	</xsl:template>


</xsl:stylesheet>
