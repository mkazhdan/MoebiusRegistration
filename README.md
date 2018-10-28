<CENTER><H1>M&ouml;bius Registration (V2.5)</H1></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#EXECUTABLE">executables</A>
<A HREF="#USAGE">usage</A>
<A HREF="#NOTES">notes</A>
<A HREF="#CHANGES">changes</A>
</CENTER>
<HR>
This distribution contains code for constructing and registering conformal spherical parametrizations of water-tight, genus-zero surfaces. Specifically, it provides implementations for:
<UL>
<LI>Computing a conformal parametrization over the sphere
<LI>Centering the parametrization with respect to M&ouml;bius inversions
<LI>Tessellating the conformal parametrization to a regular equirectangular grid
<LI>Performing fast spherical correlation to find the rotation/reflections that best aligning two centered parametrizations
<LI>Using the registered parametrizations to compute dense correspondences from a source mesh to a target
</UL>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<UL>
<B>Papers</B>: <A HREF="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">[Kazhdan, Solomon, and Ben-Chen, 2012]</A> <A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP18.pdf">[Baden, Crane, and Kazhdan, 2018]</A><br>
<B>Binaries</B>: <A href="MoebiusRegistration.x64.zip">Windows Executables</A> <A href="MoebiusRegistration.x64.lib.zip">Windows libraries and DLLs</A><br>
<B>Source Code</B>: <A href="MoebiusRegistration.zip">ZIP</A> <A HREF="https://github.com/mkazhdan/MoebiusRegistration">GitHub</A><br>
<B>Older Versions</B>: <A HREF="../Version1.0">V1.0</A><BR>
<B>License</B>: <A href="../license.txt">BSD</A><br>
</UL>

<HR>
<A NAME="EXECUTABLES"><B>EXECUTABLES</B></A><br>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>SphereMap</b></font>:
(1) Computes a conformal parametrization of a water-tight, genus-zero surface to the sphere using Conformalized Mean Curvature Flow <A HREF="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">[Kazhdan, Solomon, and Ben-Chen, 2012]</A>; (2) cannonically centers the parametrization relative to M&ouml;bius inversions (or, more generally, the low-frequency spherical harmonics); and (3) tessellates the spherical mapping over a regular equirectangular grid.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input mesh</i>&gt;</dt>
<dd> This string is the name of the file from which the mesh will be read.<br>
The file is assumed to be in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
</dd>

<dt>[<b>--out</b> &lt;<i>output triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the spherically parametrized mesh will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If the input contains colors, they will be copied to the output. Otherwise, colors are assigned using the surface normals.
</dd>

<dt>[<b>--outT</b> &lt;<i>output tessellated triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the mesh obtained by tessellating against a regular equirectangular grid will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the original vertex positions), "px", "py", "pz" (for the associated positions on the unit sphere), and "red", "green", "blue" (for the per-vertex colors). If the input contains colors, they will be copied to the output. Otherwise, colors are assigned using the surface normals.
</dd>

<dt>[<b>--outG</b> &lt;<i>output spherical grid</i>&gt;]</dt>
<dd> This string is the name of the file to which the spherical grid of conformal factors will be written.<br>
If the file extension is ".sgrid", the grid will be written as a 2D array of values (in binary). If the file extension is ".ply", the grid will be visualized as a triangle mesh obtained by scaling points on the unit sphere in proportion to their value.
</dd>

<dt>[<b>--mesh</b> &lt;<i>mesh type</i>&gt;]</dt>
<dd> This integer values specifies the type of mesh used (and how the associated mass and stiffness matrices are defined).<br>
A value of "1" indicates a triangle mesh with the standard cotangent Laplacian.<br>
A value of "2" indicates a polygon mesh using the <A HREF="https://dl.acm.org/citation.cfm?doid=1964921.1964997">Discrete Laplacians on General Polygonal Meshes</A> is used.<br>
A value of "3" indicates that a polygon mesh is used, but the mass and stiffness matrices are defined over non-triangle polygons by adding the center, triangulating by connecting the edges of the polygon to the center, using the standard cotan. Laplacian, but adding the constraint that the value at the center has to be the average of the values at the polygon corners. <br>
The default value for this parameter is 1.
</dd>

<dt>[<b>--fill</b> &lt;<i>hole fill type</i>&gt;]</dt>
<dd> This integer values specifies how holes should be handled.<br>
A value of "0" indicates that holes should be left as is.<br>
A value of "1" indicates that each hole should be filled by adding the center of the hole to the list of mesh vertices and adding the triangles defined by connecting the hole edges to the center vertex to the list of polygons.<br>
A value of "2" indicates that each hole should be filled by using the minimal area triangulation of the hole boundary. <br>
A value of "3" indicates that each hole should be filled by adding the polygon whose edges are made up of the hole boundary. (This is not supported when the <b>--mesh</b> type is set to 1.)
The default value for this parameter is 0.
</dd>


<dt>[<b>--iters</b> &lt;<i>number of CMCF iterations</i>&gt;]</dt>
<dd> This integer values specifies the number of Conformalized Mean Curvature Flow iterations to be used to obtain the conformal spherical parametrization.<br>
The default value for this parameter is 100.
</dd>

<dt>[<b>--stepSize</b> &lt;<i>the temporal size of each CMCF step</i>&gt;]</dt>
<dd> This floating point values specifies the units for the temporal discretization of the Conformalized Mean Curvature Flow.<br>
The default value for this parameter is 0.1.
</dd>

<dt>[<b>--cutOff</b> &lt;<i>M&ouml;bius centering cut-off</i>&gt;]</dt>
<dd> This floating point value specifies the threshold for terminating the M&ouml;bius centering iterations.<br>
The default value for this parameter is 10^(-10).
</dd>

<dt>[<b>--degree</b> &lt;<i>spherical harmonics degree</i>&gt;]</dt>
<dd> This integer value specifies the degrees of the spherical harmonics that should be centered out using explicit advection.<BR>
If this parameter is not specified, the code reverts to centering with respect to M&ouml;bius inversions.<BR>
Only degrees 1, 2, 3, and 4 are supported at this point.
</dd>

<dt>[<b>--aSteps</b> &lt;<i>advection steps</i>&gt;]</dt>
<dd> If a spherical harmonic degree is specified, this integer value specifies the number of advection steps to be performed within each centering step.<BR>
The default value for this parameter is 4.
</dd>

<dt>[<b>--aStepSize</b> &lt;<i>advection step size</i>&gt;]</dt>
<dd> If a spherical harmonic degree is specified, this floating point value specifies the size of each advection step.<BR>
The default value for this parameter is 0.25.<BR>
<I>If the resulting spherical parameterization exhibits triangle flips, it is likely that the advection step size should be reduced.</I>
</dd>

<dt>[<b>--res</b> &lt;<i>equirectangular grid resolution</i>&gt;]</dt>
<dd> This integer value specifies the resolution of the equirectangular grid used to tessellate the spherical parametrization.<br>
The default value for this parameter is 256.
</dd>

<dt>[<b>--smooth</b> &lt;<i>spherical diffusion time</i>&gt;]</dt>
<dd> This floating point value specifies the temporal duration for the heat diffusion used to antialias the sampled spherical function.<br>
The default value for this parameter is 0.0005.
</dd>

<dt>[<b>--c2i</b> &lt;<i>center to inversion type</i>&gt;]</dt>
<dd> This integer value specifies how the gradient descent value is to be interpreted as a center of inversion. A value of <B>0</B> indicates that a trivial interpretation is to be used. A value of <B>1</B> indicates that a golden section search should be performed along the descent direction. A value of <B>2</B> indicates that the length of centering transformation should be rescaled using the metric for the Poincar&eacute; disk model.<BR>
The default value for this parameter is 2.
</dd>

<dt>[<b>--gssTolerance</b> &lt;<i>golden section search tolerance</i>&gt;]</dt>
<dd> This floating point value specifies the tolerance for the golden section search.<BR>
The default value for this parameter is 0.000001.
</dd>

<dt>[<b>--random</b>]</dt>
<dd> If enabled, this flag specifies that the vertices of the input mesh should be assigned random positions within the unit ball before performing the Conformalized Mean Curvature Flow. (But after extracting the stiffness matrix.)
</dd>

<dt>[<b>--noCenter</b>]</dt>
<dd> If enabled, no M&ouml;bius centering is performed after computing the conformal spherical parametrization.
</dd>

<dt>[<b>--collapse</b>]</dt>
<dd> If enabled, the triangles falling into a single equirectangular cell are collapsed into a single quad before extracting the spherical tessellation.
</dd>

<dt>[<b>--ascii</b>]</dt>
<dd> If enabled, all PLY files are output in ASCII mode.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>

<dt>[<b>--fullVerbose</b>]</dt>
<dd> If enabled, even more details regarding the running times of the different stages of processing are output.
</dd>

</DETAILS>
</dl>
</ul>


<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Parametrization</b></font>:
Replaces the positions of the vertices of the triangle mesh with the parameteric positions on the unit sphere.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input mesh</i>&gt;</dt>
<dd> This string is the name of the file from which the parametrized triangle mesh will be read.<br>
The file is assumed to be in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and should contain fields "x", "y", "z" (for the original vertex positions), and "px", "py", "pz" (for the associated positions on the unit sphere).
</dd>

<dt>[<b>--out</b> &lt;<i>output triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the triangle mesh will be written.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and will contain vertices with fields "x", "y", "z" (for the positions of the parameterization on the unit sphere). If the input mesh contains per-vertex colors, these will be preserved in the output.
</dd>

</DETAILS>
</dl>
</ul>


<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>Register</b></font>:
Performs fast spherical correlation to align the rotational component of the M&ouml;bius transformation and outputs the transformed source. (The transformed source is obtained by rotating the parametric positions and, optionally, by replacing the source vertex positions with the corresponding target vertex positions.)
</SUMMARY>
<dt><b>--in</b> &lt;<i>input source/target</i>&gt;</dt>
<dd> These pair of strings are the names of the source and target file from which the spherical parameterizations will be read.<br>
The files are either both in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and should contain fields "x", "y", "z" (for the original vertex positions), and "px", "py", "pz" (for the associated positions on the unit sphere), or they should be the ".sgrid" files output by <B>SphereMap</B>.
</dd>

<dt>[<b>--out</b> &lt;<i>output triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file containing the source mesh with vertex positions on the target.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
[This output is only supported in the case that the input is in PLY format.]
</dd>

<dt>[<b>--res</b> &lt;<i>equirectangular grid resolution</i>&gt;]</dt>
<dd> This integer value specifies the resolution of the equirectangular grid used to tessellate the spherical parametrization.<br>
The default value for this parameter is 256.<br>
[This is only used in the case that the intput is in PLY format.]
</dd>

<dt>[<b>--smooth</b> &lt;<i>spherical diffusion time</i>&gt;]</dt>
<dd> This floating point value specifies the temporal duration for the heat diffusion used to antialias the sampled spherical function.<br>
The default value for this parameter is 0.0005.<br>
[This is only used in the case that the input is in PLY format.]
</dd>

<dt>[<b>--cType</b> &lt;<i>correlation type</i>&gt;]</dt>
<dd> This integer value specifies how to perform correlation. A value of <b>1</B> indicates that correlation should only be performed over the orthogonal transformations with determinant 1. A value of <B>2</B> indicates that the correlation should only be performed over the orthogonal transformations with determinant -1. A value of <B>3</B> indicates that the correlation should be performed over all orthogonal transformations.<br>
The default value for this parameter is 1.
</dd>

<dt>[<b>--correspondence</b>]</dt>
<dd> If enabled, the output mesh is defined by using the triangulation of the source and setting the vertex positions to the corresponding positions on the target. Otherwise, the parameteric coordinates of the source are rotated.<BR>
[This is only used in the case that the input is in PLY format.]
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
</dd>


</DETAILS>
</dl>
</ul>


<HR>
<A NAME="NOTES"><B>NOTES</B></A><br>
<UL>
<LI> The implementation of this code relies on the <A HREF="http://eigen.tuxfamily.org/">Eigen</A>, <A HREF="https://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A>, and <A HREF="http://www.fftw.org/">FFTW</A> libraries. The source for Eigen and SOFT are included and should compile under both Windows and Linux. For the FFTW, Windows .lib and .dll files can be found <A href="MoebiusRegistration.x64.lib.zip">here</A>.</UL>

<HR>
<DETAILS>
<SUMMARY>
<A NAME="CHANGES"><B>CHANGES</B></A><br>
</SUMMARY>
<UL>
<A HREF="../Version2.0/"><B>Version 2.0</B></A>:
<OL>
<LI> Added support for general polygonal meshes using <A HREF="https://dl.acm.org/citation.cfm?doid=1964921.1964997">Discrete Laplacians on General Polygonal Meshes</A> via the <b>--poly</b> flag.
<LI> Added support for genus-zero surfaces with boundaries by filling in the holes, via the <b>--fill</b> flag.
</OL>
<A HREF="../Version2.5/"><B>Version 2.5</B></A>:
<OL>
<LI> Expanded support for polygon mesh types with mass/stiffness matrices obtained by adding the center, triangulating to the center, and constraining the system so that the value at the center is the average value at the polygon vertices.
<LI> Added support for filling holes by adding the minimal area triangulation and add the polygonal hole to the polygon mesh.
</OL>
</UL>
</DETAILS>



<HR>
<A HREF="http://www.cs.jhu.edu/~misha">HOME</A>
