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
	  -  <edit by="David Hilvert" on="2006-Sep-24">
	  -  <edit by="David Hilvert" in-month="2006-Sep">
	  -  <edit by="David Hilvert" in-year="2006">
	  -
	  -  The first form is preferred for occasional edits; the latter two
	  -  are provided as shorthand for edits occurring over a number of
	  -  days in succession, or for cases where the exact date of the edit
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

	<!--
	  -  License information
	  -->

	<xsl:template match="*|/" mode="license-terms">
	<para>
	  This <xsl:apply-templates select="." mode="document-type"/> is free
	  documentation; you can redistribute it and/or modify it under the
	  terms of the GNU General Public License as published by the Free
	  Software Foundation; either version 2 of the License, or (at your
	  option) any later version.
	</para>
	</xsl:template>

	<xsl:template match="*|/" mode="license-object">
	<para>
	  Unless otherwise defined within the GNU General Public License, any
	  references to "object code" within the license are to be interpreted
	  to refer to any non-source version of this <xsl:apply-templates
	  select="." mode="document-type"/> (or of any work based on this
	  <xsl:apply-templates select="." mode="document-type"/>).
	</para>
	</xsl:template>

	<xsl:template match="*|/" mode="license-warranty">
	<para>
	   This <xsl:apply-templates select="." mode="document-type"/> is
	   distributed in the hope that it will be useful, but WITHOUT ANY
	   WARRANTY; without even the implied warranty of MERCHANTABILITY or
	   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
	   for more details.
	 </para>
	</xsl:template>

	<xsl:template match="*|/" mode="license-availability">
	<para>
	   You should have received a copy of the GNU General Public License
	   along with this <xsl:apply-templates select="."
	   mode="document-type"/>; if not, write to the Free Software
	   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
	</para>
	</xsl:template>

	<xsl:template match="*|/" mode="license">
	<legalnotice>
	  <xsl:apply-templates select="." mode="license-terms"/>
	  <xsl:apply-templates select="." mode="license-warranty"/>
	  <xsl:apply-templates select="." mode="license-availability"/>
	  <xsl:apply-templates select="." mode="license-object"/>
	</legalnotice>
	</xsl:template>

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
	  -  Generate title, author, copyright, and license information for
	  -  articles, books and sets.
          -->

	<xsl:template match="setinfo" mode="document-type">
	  <xsl:text>set</xsl:text>
	</xsl:template>

	<xsl:template match="bookinfo" mode="document-type">
	  <xsl:text>book</xsl:text>
	</xsl:template>

	<xsl:template match="articleinfo" mode="document-type">
	  <xsl:text>article</xsl:text>
	</xsl:template>

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


	  <!--
	    -  Add a license notice
	    -->

	    <xsl:apply-templates select="." mode="license"/>

	</xsl:copy>
	</xsl:template>

	<!-- 
	  -  Abbreviations for DocBook elements.
          -->

	<xsl:template match="ul">
	<itemizedlist>
	  <xsl:apply-templates/>
	</itemizedlist>
	</xsl:template>

	<xsl:template match="li">
	<listitem>
	  <xsl:apply-templates/>
	</listitem>
	</xsl:template>


</xsl:stylesheet>
