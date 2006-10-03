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
	  -->

	<!--
	  -  Character classes
	  -->

	<!ENTITY uppercase "'ABCDEFGHIJKLMNOPQRSTUVWXYZ'">
	<!ENTITY lowercase "'abcdefghijklmnopqrstuvwxyz'">
	
	<!--
	  -
	  -  NOTE: the below entities operate on edit records, which
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

	<!--
	  -  Sort order for years
	  -->

	<!ENTITY year-sort-order '&year;'>

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

	<xsl:param name="product-version" select="'0.8.5-prerelease'"/>

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
	  -  Capitalize the initial letter of titles.
	  -->

	<xsl:template name="write_title">
	   <xsl:param name="product-name"/>
	   <xsl:param name="product-version"/>
           <xsl:param name="title" select="."/>

	   <xsl:variable name="space-stripped-title" select="normalize-space($title)"/>
	   <xsl:variable name="initial" select="substring($space-stripped-title, 1, 1)"/>
	   <xsl:variable name="sequel" select="substring($space-stripped-title, 2)"/>

	   <title>

	     <xsl:copy-of select="$product-name"/>
	     <xsl:text> </xsl:text>
	     <xsl:copy-of select="$product-version"/>
	     <xsl:text> </xsl:text>
	     <xsl:copy-of select="concat(translate($initial, &lowercase;, &uppercase;), $sequel)"/>

           </title>
	</xsl:template>

	<xsl:template match="text()" mode="title">
	  <xsl:call-template name="write_title">
	    <xsl:with-param name="title" select="."/>
	  </xsl:call-template>
	</xsl:template>

	<!--
	  - Titles other than article, book and set titles.
	  -->

	<xsl:template match="title|t">
	  <xsl:apply-templates select="text()" mode="title"/>
	</xsl:template>

	<!-- 
	  -  Add package information to titles of articles, books and sets.
          -->

	<xsl:template match="setinfo/title|bookinfo/title|articleinfo/title|setinfo/t|bookinfo/t|articleinfo/t">
	   <xsl:call-template name="write_title">
	     <xsl:with-param name="product-name" select="$product-name"/>
	     <xsl:with-param name="product-version" select="$product-version"/>
	     <xsl:with-param name="title" select="text()"/>
	   </xsl:call-template>
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

	    <xsl:apply-templates select="title|t"/>

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
	        <xsl:sort select="&year-sort-order;"/>
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
	  -  Changelogs
	  -->

	<!--
	  -  Changelog item
	  -->

	<xsl:template match="entry">
	  <listitem>
	    <xsl:apply-templates/>
	  </listitem>
	</xsl:template>

	<!--
	  - Keys
	  -->

	<xsl:key name="taxonomy-names" match="taxonomy//*[not(keyword)]" use="name()"/>
	<xsl:key name="taxonomy-keywords" match="taxonomy//*" use="./keyword/@keyword"/>

	<!--
	  - Taxonomy
	  -->

	<xsl:variable name="taxonomy" select="document('taxonomy.xml')"/>

	<!--
	  - Attempt to map all changes to sublevels
	  -->

	<xsl:template name="map-changes-to-sublevels">
	  <xsl:param name="taxonomy"/>
	  <xsl:param name="changes"/>
	  <xsl:param name="fragment"/>
	  <xsl:param name="iteration" select="1"/>

	  <xsl:if test="count($taxonomy/*[$iteration]) = 0">
	    <xsl:copy-of select="$fragment">
	  </xsl:if>

	  <xsl:variable name="new-fragment">
	  </xsl:variable>

	  <xsl:call-template name="map-changes-to-sublevels">
	    <xsl:with-param name="taxonomy" select="$taxonomy"/>
	    <xsl:with-param name="changes" select="$changes"/>
	    <xsl:with-param name="fragment" select="$new-fragment"/>
	    <xsl:with-param name="iteration" select="$iteration + 1"/>
	  </xsl:call-template>

	</xsl:template name="map-changes-to-sublevels">

	<!--
	  - Write changes
	  -->

	<xsl:template name="write-changes">
	  <xsl:param name="taxonomy"/>
	  <xsl:param name="changes"/>

	  <!--
	    - Attempt to map changes to sublevels
	    -->

	  <xsl:variable name="map">
	    <xsl:call-template name="map-changes-to-sublevels">
	      <xsl:with-param name="taxonomy" select="$taxonomy">
	      <xsl:with-param name="changes" select="$changes">
	    </xsl:call-template>
	  </xsl:variable>

	  <xsl:choose>
	    <xsl:when test="contains($map, 'fail-to-insert-error'">
	      
	      <!--
	        - In case of failure, insert all nodes at the current level. 
		-->

	      <itemizedlist>
	      <xsl:apply-template select=".//entry">
	      </itemizedlist>

	    </xsl:when>
	    <xsl:otherwise>

	      <!--
	        - In case of success, insert all nodes within sublevels.
		-->

	      <xsl:for-each select="$taxonomy/*">
	      </xsl:for-each>

	    </xsl:otherwise>
	  </xsl:choose>

	</xsl:template>

	<!--
	  - Release
	  -->

	<xsl:template match="release">
	    <section>
		<title><xsl:value-of select="@version"/></title>

		<!--
		  - Notes
		  -->

		<xsl:apply-templates select=".//note"/>

		<!--
		  - Write changes according to the change taxonomy
		  -->

		<xsl:call-template name="write-changes">
		  <xsl:with-param name="taxonomy" select="document('taxonomy.xml')"/>
		  <xsl:with-param name="changes" select="."/>
		</xsl:call-template>
            </section>
	</xsl:template>

	<xsl:template match="changelog">
	  <xsl:apply-templates select="release"/>
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

	<xsl:template match="s">
	<section>
	  <xsl:apply-templates/>
	</section>
	</xsl:template>


</xsl:stylesheet>
