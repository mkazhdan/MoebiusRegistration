<CENTER><H1>Moebius Registration</H1></CENTER>
<CENTER>
<A HREF="#LINKS">links</A>
<A HREF="#EXECUTABLE">executables</A>
<A HREF="#USAGE">usage</A>
<A HREF="#NOTES">notes</A>
</CENTER>
<HR>
This distribution contains code for constructing and registering conformal spherical parametrizations of water-tight, genus-zero surfaces. Specifically, it provides implementations for:
<UL>
<LI>Computing a conformal parametrization over the sphere
<LI>Centering the parametrization with respect to Moebius inversions
<LI>Tessellating the conformal parametrization to a regular equirectangular grid
<LI>Performing fast spherical correlation to find the rotation/reflections that best aligning two centered parametrizations
<LI>Using the registered parameerizations to compute dense correspondences from a source mesh to a target
</UL>
<HR>
<A NAME="LINKS"><B>LINKS</B></A><br>
<A href="http://www.cs.jhu.edu/~misha/MyPapers/SGP18.pdf">SGP 2018 Paper</A><br>
<A href="MoebiusRegistration.x64.zip">Windows Executables</A><br>
<A href="MoebiuseRegistration.x64.lib.zip">Windows libraries and DLLs</A><br>
<A href="MoebiusRegistration.zip">Source Code</A><br>
<A href="license.txt">License</A><br>

<HR>
<A NAME="EXECUTABLES"><B>EXECUTABLES</B></A><br>

<ul>
<dl>
<DETAILS>
<SUMMARY>
<font size="+1"><b>SphereMap</b></font>:
(1) Computes a conformal parametrization of a water-tight, genus-zero surface to the sphere using Conformalized Mean Curvature Flow <A HREF="http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf">[Kazhdan, Solomon, and Ben-Chen, 2012]</A>; (2) cannoically centers the parametrization relative to Moebius inversions; and (3) tessellates the spherical mapping over a regular equirectangular grid.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input mesh</i>&gt;</dt>
<dd> This string is the name of the file from which the point set will be read.<br>
The file is assumed to be in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
</dd>

<dt>[<b>--out</b> &lt;<i>output triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file to which the triangle mesh will be written.<br>
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

<dt>[<b>--iters</b> &lt;<i>number of CMCF iterations</i>&gt;]</dt>
<dd> This integer values specifies the number of Conformalized Mean Curvature Flow iterations to be used to obtain the conformal spherical parametrization.<br>
The default value for this parameter is 100.
</dd>

<dt>[<b>--stepSize</b> &lt;<i>the temporal size of each CMCF step</i>&gt;]</dt>
<dd> This floating point values specifies the units for the temporal discretization of the Conformalized Mean Curvature Flow.<br>
The default value for this parameter is 0.1.
</dd>

<dt>[<b>--cutOff</b> &lt;<i>Moebius centering cut-off</i>&gt;]</dt>
<dd> This floating point value specifies the threshold for terminating the Moebius centering iterations.<br>
The default value for this parameter is 10^(-10).
</dd>

<dt>[<b>--res</b> &lt;<i>equirectangular grid resolution</i>&gt;]</dt>
<dd> This integer value specifies the resolution of the equirectangular grid used to tessellate the spherical parametrization.<br>
The default value for this parameter is 256.
</dd>

<dt>[<b>--smooth</b> &lt;<i>spherical diffusion time</i>&gt;]</dt>
<dd> This floating point value specifies the temporal duration for the heat diffusion used to antialias the sampled spherical function.<br>
The default value for this parameter is 0.0005.
</dd>

<dt>[<b>--random</b>]</dt>
<dd> If enabled, this flag specifies that the vertices of the input mesh should be assigned random positions within the unit ball before performing the Conformalized Mean Curvature Flow. (But after extracting the stiffness matrix.)
</dd>

<dt>[<b>--noCenter</b>]</dt>
<dd> If enabled, no Moebius centering is performed after computing the conformal spherical parametrization.
</dd>

<dt>[<b>--collapse</b>]</dt>
<dd> If enabled, the triangles are falling into a single equirectangular cell are collapsed into a single quad before extracting the spherical tessellation.
</dd>

<dt>[<b>--ascii</b>]</dt>
<dd> If enabled, all PLY files are output in ASCII mode.
</dd>

<dt>[<b>--verbose</b>]</dt>
<dd> If enabled, details regarding the running times of the different stages of processing are output.
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
Performs fast spherical correlation to align the rotational component of the Moebius transformation and computes dense point-to-point correspondences between the source and target mesh.
</SUMMARY>
<dt><b>--in</b> &lt;<i>input source/target</i>&gt;</dt>
<dd> These pair of strings are the names of the source and target file from which the spherical parameterizations will be read.<br>
The files are either both in the <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format and should contain fields "x", "y", "z" (for the original vertex positions), and "px", "py", "pz" (for the associated positions on the unit sphere), or they should be the ".sgrid" files output by <B>SphereMap</B>.
</dd>

<dt>[<b>--out</b> &lt;<i>output triangle mesh</i>&gt;]</dt>
<dd> This string is the name of the file containing the source mesh with vertex positions on the target.<br>
The file is written in <a href="http://www.cc.gatech.edu/projects/large_models/ply.html">PLY</a> format.<br>
Note that this output is only supported in the case that the input is in PLY format.
</dd>

<dt>[<b>--res</b> &lt;<i>equirectangular grid resolution</i>&gt;]</dt>
<dd> This integer value specifies the resolution of the equirectangular grid used to tessellate the spherical parametrization.<br>
The default value for this parameter is 256.<br>
Note that this is only used in the case that the intput is in PLY formate.
</dd>

<dt>[<b>--smooth</b> &lt;<i>spherical diffusion time</i>&gt;]</dt>
<dd> This floating point value specifies the temporal duration for the heat diffusion used to antialias the sampled spherical function.<br>
The default value for this parameter is 0.0005.<br>
Note that this is only used in the case that the input is in PLY format.
</dd>

<dt>[<b>--cType</b> &lt;<i>correlation type</i>&gt;]</dt>
<dd> This integer value specifies how to perform correlation. A vlaue of <b>1</B> indicates that correlation should only be performed over the orthogonal transformations with determinant 1. A value of <B>2</B> indicates that the correlation should only be performed over the orthogonal transformations with determinant -1. A value of <B>3</B> indicates that the correlation should be performed over all orthogonal transformations.
The default value for this parameter is 1.
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
<LI> The implementation of this code relies on the <A HREF="http://eigen.tuxfamily.org/">Eigen</A>, <A HREF="https://www.cs.dartmouth.edu/~geelong/soft/">SOFT</A>, and <A HREF="http://www.fftw.org/">FFTW</A>. The source for Eigen and SOFT are included and should compile uncer both Windows and Linux. For the FFTW, Windows .lib and .dll files can be found <A href="MoebiuseRegistration.x64.lib.zip">here</A>.
</UL>


<HR>
<A HREF="http://www.cs.jhu.edu/~misha">HOME</A>