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
	<xsl:param name="product-version" select="'0.8.5'"/>
	<xsl:param name="site-URL" select="
		concat('http://auricle.dyndns.org/', $product-name, '/')"/>
	<xsl:param name="download-URL" select="concat($site-URL, 'download/')"/>
	<xsl:param name="windows-binary-package-name" select="
		concat(translate($product-name, &uppercase;, &lowercase;), '-', 
		       translate($product-version, '.', '_'), '-win32.zip')"/>
	<xsl:param name="source-package-name" select="
		concat(translate($product-name, &uppercase;, &lowercase;), '-', $product-version)"/>
	<xsl:param name="source-package-name-tar-gz" select="
		concat($source-package-name, '.tar.gz')"/>
	<xsl:param name="windows-binary-URL" select="concat($download-URL, $windows-binary-package-name)"/>
	<xsl:param name="source-URL" select="concat($download-URL, $source-package-name-tar-gz)"/>
	<xsl:param name="mailing-list-address" select="'ale@ventricle.dyndns.org'"/>

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
	  Unless explicitly defined within the applicable version of the GNU
	  General Public License, any references to "object code" within the
	  license shall refer to any non-source version of this
	  <xsl:apply-templates select="." mode="document-type"/> (or of any
	  work based on this <xsl:apply-templates select="."
	  mode="document-type"/>).
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
	   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
	   02110-1301, USA.
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
	  <xsl:variable name="objectroot" select=".."/>
	  <xsl:variable name="editors-unique" select="$objectroot//edit[count(.|key('editors', &editor;)[&scope;][1]) = 1]"/>
	  <xsl:variable name="editor-years-unique" select="$objectroot//edit[count(.|key('editor-years', &editor-year;)[&scope;][1]) = 1]"/>

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
	      <xsl:variable name="this-editor" select="&editor;"/>
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
	  -  Changelogs and news files.
	  -->

	<!--
	  -  Changelog item
	  -->

	<xsl:template match="entry">
	  <listitem>
	    <xsl:apply-templates select="text/text()"/>
	  </listitem>
	</xsl:template>

	<!--
	  - Keys
	  -->

	<xsl:key name="word-map" match="entry" use="word/child::text()"/>

	<!--
	  - Taxonomy
	  -->

	<xsl:variable name="taxonomy-root" select="document('taxonomy.xmli')"/>

	<!--
	  - Write a fragment substring for a (change, category) pair.
	  -->

	<xsl:template name="write-change-category">
	  <xsl:param name="change"/>
	  <xsl:param name="category"/>
	  <xsl:value-of select="concat(generate-id($change), ':', generate-id($category), ' ')"/>
	</xsl:template>

	<!--
	  - Get a change category from a fragment result
	  -->

	<xsl:template name="get-change-category">
	  <xsl:param name="change"/>
	  <xsl:param name="fragment"/>
	  <xsl:value-of select="substring-before(substring-after(substring-after($fragment, generate-id($change)), ':'), ' ')"/>
	</xsl:template>

	<!--
	  - Get a tree fragment with categories set to zero.
	  -->

	<xsl:template name="get-zero-categories">
	  <xsl:param name="changes"/>
	    <xsl:for-each select="$changes">
	      <xsl:call-template name="write-change-category">
	        <xsl:with-param name="change" select="."/>
	        <xsl:with-param name="category" select="$taxonomy-root"/>
	      </xsl:call-template>
	    </xsl:for-each>
	</xsl:template>

	<!--
	  - Determine whether a category tree is valid (no undefineds) for all
	  - changes.  Returns a tree fragment including the string 'fail'
	  - upon failure.
	  -->

	<xsl:template name="categories-valid">
	  <xsl:param name="changes"/>
	  <xsl:param name="fragment"/>
	    <xsl:for-each select="$changes">
	      <xsl:variable name="category">
	        <xsl:call-template name="get-change-category">
	          <xsl:with-param name="change" select="."/>
	          <xsl:with-param name="fragment" select="$fragment"/>
	        </xsl:call-template>
	      </xsl:variable>

	      <xsl:if test="$category = ''">
	        <xsl:text>fail: </xsl:text>
		<apply-templates select="."/>
              </xsl:if>
	    </xsl:for-each>
	</xsl:template>

	<!--
	  - make a map unique
	  -->
	<xsl:template name="make-map-unique">
	  <xsl:param name="map"/>
	  <xsl:param name="map-prefix"/>

	  <xsl:choose>
	    <xsl:when test="$map = ''">
	    </xsl:when>
	    <xsl:when test="contains($map-prefix, substring-before($map, ':'))">
	      <xsl:call-template name="make-map-unique">
	        <xsl:with-param name="map" select="substring-after($map, ' ')"/>
	        <xsl:with-param name="map-prefix" select="$map-prefix"/>
	      </xsl:call-template>
	    </xsl:when>
	    <xsl:otherwise>
	      <xsl:value-of select="substring-before($map, ' ')"/>
	      <xsl:text> </xsl:text>
	      <xsl:call-template name="make-map-unique">
	        <xsl:with-param name="map" select="substring-after($map, ' ')"/>
	        <xsl:with-param name="map-prefix" select="concat($map-prefix, substring-before($map, ' '))"/>
	      </xsl:call-template>
	    </xsl:otherwise>
	  </xsl:choose>
	</xsl:template>

	<!--
	  - Attempt to map all changes to sublevels
	  -->

	<xsl:template name="map-changes-to-sublevels">
	  <xsl:param name="taxonomy"/>
	  <xsl:param name="changes"/>

          <xsl:variable name="non-unique-map">
	  <xsl:for-each select="$taxonomy">
	    <xsl:sort select="@match-priority" order="descending"/>

	      <xsl:variable name="category" select="."/>

	      <xsl:variable name="keywords" select=".//keyword"/>

	        <xsl:for-each select="$changes">

		  <xsl:variable name="change" select="."/>

		  <xsl:choose>
		    <xsl:when test="name($category) = name(.)">
	              <xsl:call-template name="write-change-category">
	                <xsl:with-param name="change" select="."/>
	                <xsl:with-param name="category" select="$category"/>
	              </xsl:call-template>
		    </xsl:when>
		    <xsl:when test="count(key('word-map', $keywords/child::text()))
		                  = count(key('word-map', $keywords/child::text())|.)">
	              <xsl:call-template name="write-change-category">
	                <xsl:with-param name="change" select="."/>
	                <xsl:with-param name="category" select="$category"/>
	              </xsl:call-template>
		    </xsl:when>
		  </xsl:choose>
	        </xsl:for-each>
	      </xsl:for-each>
	    </xsl:variable>

	    <xsl:call-template name="make-map-unique">
	      <xsl:with-param name="map" select="$non-unique-map"/>
	    </xsl:call-template>
	</xsl:template>

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
	      <xsl:with-param name="taxonomy" select="$taxonomy"/>
	      <xsl:with-param name="changes" select="$changes"/>
	    </xsl:call-template>
	  </xsl:variable>

	  <xsl:variable name="map-failures">
	    <xsl:call-template name="categories-valid">
	      <xsl:with-param name="changes" select="$changes"/>
	      <xsl:with-param name="fragment" select="$map"/>
	    </xsl:call-template>
	  </xsl:variable>

	  <xsl:choose>
	    <xsl:when test="contains($map-failures, 'fail')">
	      
	      <!--
	        - In case of failure, insert all nodes at the current level. 
		-->

	      <itemizedlist>
	        <xsl:apply-templates select="$changes"/>
	      </itemizedlist>

	    </xsl:when>
	    <xsl:otherwise>

	      <!--
	        - In case of success, insert all nodes within sublevels.
		-->

	      <xsl:for-each select="$taxonomy">
	        <xsl:variable name="category" select="."/>
		<xsl:if test="contains($map, concat(generate-id($category), ' '))">
	          <section tocexclude="1">
		    <xsl:call-template name="write_title">
		      <xsl:with-param name="title">
		        <xsl:choose>
		        <xsl:when test="@t != ''">
		          <xsl:value-of select="@t"/>
			</xsl:when>
			<xsl:otherwise>
		          <xsl:value-of select="translate(name($category), '-', ' ')"/>
			</xsl:otherwise>
			</xsl:choose>
		      </xsl:with-param>
		    </xsl:call-template>
		    <xsl:call-template name="write-changes">
		      <xsl:with-param name="taxonomy" select="$category/*"/>
		      <xsl:with-param name="changes" select="$changes[name(.) = name($category)]/*|$changes[not(name(.) = name($category))][contains($map, concat(generate-id(.), ':', generate-id($category), ' '))]"/>
		    </xsl:call-template>
		  </section>
		</xsl:if>
	      </xsl:for-each>

	    </xsl:otherwise>

	  </xsl:choose>

	</xsl:template>

	<!--
	  - News entries.
	  -->

	<xsl:template match="fm">
	  <para>
	    <xsl:choose>
	      <xsl:when test="@nh">
	        <xsl:apply-templates/><xsl:text> (Freshmeat announcement via Neohapsis)</xsl:text>
	      </xsl:when>
	      <xsl:otherwise>
	        <xsl:apply-templates/><xsl:text> (Freshmeat announcement)</xsl:text>
	      </xsl:otherwise>
	    </xsl:choose>
	  </para>
	</xsl:template>

	<xsl:template match="sum">
	  <section tocexclude="1"><title>Program summary</title>
	    <xsl:if test="@revised">
	      <para>This release includes a revised summary:</para>
		<!-- This release is accompanied by a revised summary: -->
	    </xsl:if>
	    <para>
	      <xsl:apply-templates/>
	    </para>
	  </section>
	</xsl:template>

	<xsl:template match="ed-note">
	  <xsl:apply-templates/>
	</xsl:template>

	<xsl:template match="notes">
	  <section tocexclude="1"><title>Notes</title>
	  <para>
	    <xsl:apply-templates/>
	  </para>
	  </section>
	</xsl:template>

	<xsl:template match="ch">
	  <section tocexclude="1"><title>Changelog summary</title>
	    <xsl:apply-templates/>
	  </section>
	</xsl:template>

	<xsl:template match="ml">
	  <section tocexclude="1"><title>Mailing list announcement</title>
	    <xsl:apply-templates/>
	  </section>
	</xsl:template>

	<!--
	  - Releases
	  -->

	<xsl:template match="release" name="release-news" mode="news">
	    <xsl:param name="label" select="1"/>
	    <section label="{$label}" tocexclude="1">
	    <!-- <section label="{@version}"> -->
	    	<xsl:choose>
	    	  <xsl:when test="@date">
		    <title>
		      <xsl:text>Version </xsl:text>
		      <xsl:value-of select="@version"/>
		      <xsl:text>, </xsl:text>
		      <xsl:value-of select="@date"/>
		    </title>
		  </xsl:when>
	    	  <xsl:otherwise>
		    <title>
		      <xsl:text>Version </xsl:text>
		      <xsl:value-of select="@version"/>
		    </title>
		  </xsl:otherwise>
		</xsl:choose>

		<!--
		  - Notes
		  -->

		<xsl:apply-templates select=".//note"/>

		<!--
		  - Write subsections.
		  -->

		<xsl:apply-templates/>

            </section>
	</xsl:template>

	<xsl:template match="release" name="release-changelog" mode="changelog">
	    <xsl:param name="label" select="1"/>
	    <section label="{$label}" tocexclude="1">
	    <!-- <section label="{@version}"> -->
	    	<xsl:choose>
	    	  <xsl:when test="@date">
		    <title>
		      <xsl:text>Version </xsl:text>
		      <xsl:value-of select="@version"/>
		      <xsl:text>, </xsl:text>
		      <xsl:value-of select="@date"/>
		    </title>
		  </xsl:when>
	    	  <xsl:otherwise>
		    <title>
		      <xsl:text>Version </xsl:text>
		      <xsl:value-of select="@version"/>
		    </title>
		  </xsl:otherwise>
		</xsl:choose>

		<!--
		  - Notes
		  -->

		<xsl:apply-templates select=".//note"/>

		<!--
		  - Write changes according to the change taxonomy
		  -->

		<xsl:call-template name="write-changes">
		  <xsl:with-param name="taxonomy" select="$taxonomy-root/taxonomy/*"/>
		  <xsl:with-param name="changes" select="./*[name() != 'edit']"/>
		</xsl:call-template>
            </section>
	</xsl:template>

	<xsl:template match="news">
	  <xsl:for-each select="release">
	    <xsl:call-template name="release-news">
	      <xsl:with-param name="label" select="count(../release) - position()"/>
	    </xsl:call-template>
	  </xsl:for-each>
	</xsl:template>

	<xsl:template match="changelog">
	  <xsl:for-each select="release">
	    <xsl:call-template name="release-changelog" mode="changelog">
	      <xsl:with-param name="label" select="count(../release) - position()"/>
	    </xsl:call-template>
	  </xsl:for-each>
	</xsl:template>

	<!--
	  - Inline product information
	  -->

	<xsl:template match="winpack">
	  <xsl:value-of select="$windows-binary-package-name"/>
	</xsl:template>

	<xsl:template match="ver">
	  <xsl:value-of select="$product-version"/>
	</xsl:template>

	<xsl:template match="sourcepack">
	  <xsl:value-of select="$source-package-name"/>
	</xsl:template>

	<xsl:template match="sourcepacktargz">
	  <xsl:value-of select="$source-package-name-tar-gz"/>
	</xsl:template>

	<xsl:template match="winurl">
	  <xsl:choose>
	    <xsl:when test="contains($product-version, 'pre')">
	      <note>
	        Windows binaries are not available for this prerelease version.
	      </note>
	    </xsl:when>
	    <xsl:otherwise>
	      <ulink url="{$windows-binary-URL}"/>
	    </xsl:otherwise>
	  </xsl:choose>
	</xsl:template>

	<xsl:template match="sourceurl">
	  <ulink url="{$source-URL}"/>
	</xsl:template>

	<xsl:template match="mailinglist">
	  <xsl:value-of select="$mailing-list-address"/>
	</xsl:template>

	<!-- 
	  -  Abbreviations for DocBook elements.
          -->

	<xsl:template match="p">
	<para>
	  <xsl:apply-templates/>
	</para>
	</xsl:template>

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

	<xsl:template match="code">
	<literal>
	  <xsl:apply-templates/>
	</literal>
	</xsl:template>

	<xsl:template match="l">
	<literal>
	  <xsl:apply-templates/>
	</literal>
	</xsl:template>

	<xsl:template match="pre">
	<literallayout class="monospaced">
	  <xsl:apply-templates/>
	</literallayout>
	</xsl:template>

	<xsl:template match="ll">
	<literallayout class="monospaced">
	  <xsl:apply-templates/>
	</literallayout>
	</xsl:template>

	<xsl:template match="ui">
	<userinput>
	  <xsl:apply-templates/>
	</userinput>
	</xsl:template>

	<xsl:template match="meta">
	<emphasis>
	  <xsl:apply-templates/>
	</emphasis>
	</xsl:template>


</xsl:stylesheet>
