===========================================================================
         ________    _________   _________     __     _________
        |        |\ |        |\ |        |\   |  |\  /         \\
        |    ____|| |    ____|| |    __  ||   |__|| |    ___    ||
        |   |\___\| |   |\___\| |   |\_| ||    \_\| |   //  \   ||
        |   ||      |   ||      |   || | ||    __   |  ||    |  ||
        |   ||__    |   ||__    |   ||_| ||   |  |\ |  ||    |  ||
        |       |\  |       |\  |         \\  |  || |  ||    |  ||
        |    ___||  |    ___||  |    ___   || |  || |  ||    |  ||
        |   |\__\|  |   |\__\|  |   |\__|  || |  || |  ||    |  ||
        |   ||      |   ||      |   ||  |  || |  || |  ||    |  ||
        |   ||      |   ||___   |   ||__|  || |  || |   \\__/   ||
        |   ||      |        |\ |          || |  || |           ||
        |___||      |________|| |__________|| |__||  \_________//

      F I N I T E   E L E M E N T S   F O R   B I O M E C H A N I C S

                   --- v e r s i o n - 1 . 2 . 2 ---


  Musculoskeletal Research Laboratory
  University of Utah
  http://mrl.sci.utah.edu

  copyright (c) 2006-2010 - All rights reserved

===========================================================================

R E L E A S E   N O T E S

===========================================================================

This document describes the features that have been added to the latest
FEBio release as well as the more important bug fixes. Please consult
the FEBio User’s Manual for a more detailed description of the new features.

===========================================================================

 R E L E A S E   1.2.2

===========================================================================

1. FEATURE: micro-material
--------------------------
Using the new micro-material, FEBio has taken a first step in the solution
of computational homogenization problems. For such problems, the macro-scale
problem's material response is described by the iterative solution of a 
micro-scale problem. 
In FEBio a new material, called "micro-material" is used to solve computational
homogenization problems. The material references the micro model which is 
specified in a separate input file. 


2. FEATURE: relative kinematics
-------------------------------
An option was added to describe rigid body kinematics in a local coordinate system
and connect the rigid bodies in a child-parent hierarchy. This facilitates the
description of complex rigid body prescribed kinematics. 


3. FEATURE: set minimum residual
--------------------------------
Whenever the residual drops below a minimal value, FEBio assumes that there are no
force acting on the system and will consider the time step as converged. The user can
now set the value of this minimum residual using the <min_residual> parameter in
the <Control> section of the input file. The default value is 1e-20. 


4. FEATURE: new constraints section
-----------------------------------
A new section was added that will deal with constraints. Although all constraints will
eventually be moved to this section, for now it was added to allow the user to describe
rigid kinematics (see item 5).


5. FEATURE: redefine rigid kinematics in multi-step analysis
------------------------------------------------------------
The user can now redefine the rigid body constraints using the Constraints section in
a multi-step analysis. The Constraints section may be repeated in each step, and the user
may redefine the rigid body constraints in each step.


6. FEATURE: parameter optimization
----------------------------------
A new parameter optimization algorithm was implemented. Although parameter optimization
has been implemented for a while, it was not available since it used the NAG library. 
The NAG routines have now been replaced with a custom implementation of the Levenberg-
Marquardt method, making the parameter optimization module available for all FEBio users.


7. FEATURE: ehanced loadcurves
------------------------------
Two new features were added to the load curves. First, the user can define the load curve
type, which can be set to step, linear or smooth. For the smooth type, FEBio will inter-
polate the load curve values with a cubic polynomial creating smooth curves (C1 continuous).
Second, the user can now also define the extend mode which defines the values of the load
curve outside its specified domain. The new options are set as attributes to the <loadcurve>
element. For example, the following defines a smooth load curve with linear extrapolation
outside the domain.

<loadcurve type="smooth" extend="extrapolate">

The extend mode can take on the following values: "constant", "extrapolate", "repeat",
"repeat offset". The default value is "constant". 


8. FEATURE: nodally integrated tet
----------------------------------
A new nodally integrated tet element is available which performs better in incompressible
and bending problems. Performs better in the sense that it suffers less from locking 
problems that often cause problems using the default linear tet elements in FEBio.


9. FEATURE: constant surface traction boundary condition
--------------------------------------------------------
A new boundary condition was implemented which allows the user to apply a constant traction
force to a surface. Unlike the pressure boundary force, the surface traction will remain
independent of the deformation.


10. BUG: crash when using triangular shells
-------------------------------------------
A bug was found and fixed in the input routine for triangular shells. The bug was responsible
for unexpected crashes that could occur anywhere during the run. 



===========================================================================

 R E L E A S E   1.2.1

===========================================================================


1. FEATURE: mat_axis
--------------------
A mat_axis option was added for the ElementData section. The user can now 
override the material axis for each element individually.

===========================================================================

 R E L E A S E   1.2.0

===========================================================================


1. FEATURE: sparse linear solvers
---------------------------------
FEBio now adds support for several fast sparse linear solvers, such as Pardiso
and SuperLU. The preferred linear solver can either be specified in the 
configuration file or in the FEBio input file by adding the following line
to the Control section (e.g. for Pardiso):

<linear_solver type="pardiso"></linear_solver>

These new sparse linear solvers are faster than the default skyline solver
since they take better advantage of the sparsity of the stiffness matrix and
can use multiple processors. Note that if you want to use multiple processors
you may need to set specific environment variables. We refer to the documentation
of the linear solvers for more details. 


2. FEATURE: contact improvements
--------------------------------
Two new frictionless sliding contact algorithms were implemented which add more 
stability to contact analysis. The first is referred to as the facet-to-facet
sliding, and the second one is reffered to as sliding2. To select a particular
algorithm, set the type attribute of the contact xml-element to the appropriate
value (e.g. for facet-to-facet sliding). 

<contact type="facet-to-facet sliding"> ... </contact>

For sliding2, simply enter the value "sliding2". The sliding2 algorithm is
inherently unsymmetric, which means you'll need to use a non-symmetric linear
solver (both pardiso and superlu can handle non-symmetric matrices). However,
this algoritm does have a symmetric mode, although this mode was shown not to
be as robust as the non-symmetric. You can choose between the two modes by
setting the following flag in the contact section.

<symmetric_stiffness>1</symmetric_stiffness>

A value of 0 will use the non-symmetric mode. Please see the user's manual for
more details on these new contact algorithms.


3. FEATURE: Biphasic contact
----------------------------
The new sliding contact algorithm sliding2 can also deal with biphasic contact.
Here two poroelastic materials are brought in contact and can exchange fluid
across the contact interface. The contact algorithm also automatically enforces
a free-draining boundary condition on the part of the surface that is not in
contact. To use biphasic contact you'll need to define two contacting poroelastic
surfaces and you need to set the contact type to "sliding2".

<contact type="sliding2"> ... </contact>

In addition, you'll need to define several control parameters which are described
in detail in the user's manual.


4.FEATURE: auto-penalty
-----------------------
A new auto-penalty algorithm was implemented for all sliding contact algorithms. 
This algorithm derives the penalty factor from the initial elasticity tensor of
the material. The advantage is that the algorithm now works with any material. 
To use the new algorithm, add the following xml-element to the contact section
of the input file.

<auto_penalty>1</auto_penalty>

Note that the old algorithm is still available for now for backward compatibility
but is considered obselete and may be removed in future releases.


5. FEATURE: isotropic elastic
-----------------------------
A new material was implemented to replace the linear elastic and the St.Venant-
Kirchhoff material. This material is valid for both small deformations and 
large deformations and reduces to linear elasticity for small strains. An
example illustrating the use of this material:

<material id="1" type="isotropic elastic">
	<E>1.0</E>
	<v>0.45</v>
</material>

Note that the linear elastic and the St.Venant-Kirchhoff materials are still
available for backward compatibility, but it is highly recommended to use
this new material instead.


6. FEATURE: convergence tolerances
----------------------------------
A simple mechanism was implemented to turn off particular convergence tolerances.
As of this version, setting a value of 0 (zero) for the displacement (dtol),
energy (etol) or residual (rtol) will deactivate the corresponding convergence 
criteria. Similarly, setting a value of zero for the line search tolerance (lstol),
will deactivate the line search.


7. FEATURE: time command
------------------------
The time command will display the time that has passed so far and also gives an
estimate of the remaining time. The command can be entered at the FEBio command 
prompt which is activated after interrupting the run with ctrl+c.


8. BUG: running Linux input files on Windows
--------------------------------------------
An issue was identified that prevented files created on a linux machine to be run
on windows. The problem had to do with way the different platforms deal with line
endings. The issue has been resolved, and both platforms should now be able to read
files that were created on the other.


===========================================================================

 R E L E A S E   1.1.7

===========================================================================

1. FEATURE: Initial nodal velocities
------------------------------------
Support fot initial nodal velocities was added to FEBio and to the FEBio
input file. To set initial nodal velocities, add the Initial section to the
input file. 

<febio_spec>
        ...
	<Initial>
		<velocity>
			<node id="1">1,0,0</node>
			...
			<node id="n">1,0,0</node>
		</velocity>
	</Initial>
</febio_spec>

Note that you need to run the problem dynamically for this to have any effect.
To run a dynamic analysis, add the following element to the Control section
of your input file.

<analysis type="dynamic"></analysis>


2. FEATURE: strain-dependant poroelasticity
------------------------------------------
A new poroelastic material was added that has a strain dependant permeability
tensor. This material uses the constitutive model proposed by Holmes and Mow.


3. FEATURE: customizing plot file output
----------------------------------------
You can now replace the predefined LSDYNA data fields with FEBio datafields. 
To do this, define the plotfile element in the Output section of your input file. 
You can then define maps to redefined the meaning of an output field in the plot
file. For example, 

<Output>
  <plotfile>
    <map field="velocity">CONTACT_TRACTION</map>
  </plotfile>
</Output>

will replace the velocity data with contact traction data.


4. FEATURE: frictional contact
------------------------------
The sliding interface is augmented with a frictional algorithm. To use friction
simply add the following two parameters to the your contact section.

<fric_coeff>0.1</fric_coeff>
<fric_penalty>1</fric_penalty>

The first parameter sets the friction coefficient. The second parameter is the
frictional penalty parameter.


5. FEATURE: contact stiffness multipliers
-----------------------------------------
The user can set a scalar multiplier for both the normal as well as the 
tangential stiffness. The "ktmult" parameter defines the multiplier for the
tangential component (default = 0), and "knmult" the multiplier for the 
normal component. These parameters need to be set in the sliding interface
contact section of your input file.

<contact type="sliding_with_gaps">
  ...
  <knmult>1</knmult>
  <ktmult>0</ktmult>
  <surface>
  ...
</contact>


6. FEATURE: configuration file
------------------------------
FEBio now use a configuration file to store platform dependant parameters. 
Currently it only stores the default linear solver you wishes to use which
is set to Pardiso on most platforms. The configuration file needs to be
stored in the same location of the FEBio executable. See the manual for more
options about using this configuration file. 

7. BUG: node numbering
----------------------
When a negative jacobian occurs the node numbering were reported as zero-based
in stead of one-based. This has been fixed.


===========================================================================

 R E L E A S E   1.1.6

===========================================================================

1.Multi-step Analysis
---------------------
The user can now split the analysis in separate steps where in each step
he can define different control parameters and boundary conditions. This
is useful for instance for defining time dependant boundary conditions or
for switching between (quasi-) static and dynamic analysis. Multi-step
problems are defined slightly different. For each step you need to define
a <Step> section in the input file. The overall structure of the input file
now looks like this.

<febio_spec>
	<Control>
		<!-- global control parameters (e.g. title) --!>
	</Control>
	<Material>
		<!-- material data --!>
	</Material>
	<Boundary>
		<!-- global boundary constraints --!>
	</Boundary>
	<Step>
		<Control>
			<!-- local control parameters --!>
		</Control>
		<Boudary>
			<!-- local boundary constraints --!>
		</Boundary>
	</Step>
	<Step>
		...
	</Step>
</febio_spec>

In each step section, you can redefine the control parameters, as well as 
specify boundary conditions that will only be enforced during this step.
Note that the boundary constraints defined outside the step section are
considered global and will remain active during all steps.


2.FEATURE: Ogden material
-------------------------
The Ogden material was added to the FEBio's material library. To use this
material set the material type to "Ogden". Please see the user manual for
a description of the material parameters.

<material id="1" type="Ogden">
	<c1>1</c1>
	<m1>2</m1>
	<k>100</k>
</material>

Note that this material is formulated using an uncoupled formulation. Incom-
pressibility can be enforced by either increasing the bulk modulus (k) and
by setting the augmentation flag (<laugon>) to "1".


3.BUG FIX: contact
------------------
A bug was fixed in the frictionless contact algorithm. Improved convergence
can be expected for some contact problems. 

===========================================================================

 R E L E A S E   1.1.5

===========================================================================

1. FEATURE: visco-elasticity
----------------------------
Visco-elasticity has been added to the material library. To use the visco-elastic
material you need to add two materials: (1) the visco-elastic material and (2)
a material describing the elastic equilibrium repsonse.

<material id="1" type="neo-Hookean">
	<E>1</E>
	<v>0.33</v>
</material>
<material id="2" type="visco-elastic">
	<t1>0.01</t1>
	<g1>1</g1>
	<solid_id>1</solid_id>
</material>

The solid_id element references the material id of the elastic material. 
The visco-elastic material is described using a Prony series. Each term
in this series is parameterized by a relaxation time (ti) and a scalar 
modulus gi. Up to six terms can be prescribed.


2. FEATURE: continuous fiber material
-------------------------------------
A new material model was added that has a 2D random fiber distribution. 
This model replaces the "random" fiber type that could be used to describe
random fiber materials. Currently, this material uses a Mooney-Rivlin matrix.

3.FEATURE: body forces for shells
----------------------------------
Body forces can now also be allowed to shell geometry.

4.FEATURE: linear orthotropic material
------------------------------------
A linear orthotropic material was added to FEBio's material library.

5.FEATURE: linear constraints
-----------------------------
The user can now apply linear constraints to the problem. A linear constraint
defines a degree of freedom in function of the other dofs. It is always of 
the form.

	ui = A1*u1+ ... + A(i-1)*u(i-1) + A(i+1)*u(i+1) + ... + An*un

As many linear constraints can be applied as desired. Only the non-zero
coefficients need to be specified.

6.FEATURE: SuperLU solver
-------------------------
The SuperLU solver has been implemented. To select it, add the following
line to the Control section of the input file.

<linear_solver type="superlu"></linear_solver>

7.BUG FIX: off-diagonal strains in logfile
------------------------------------------
The off-diagonal strians reported in the logfile were off by a constant
term of 0.5.


===========================================================================

 R E L E A S E   1.1.4

===========================================================================

1. FEATURE: new muscle and tendon material models
-------------------------------------------------
The muscle and tendon material models due to Silvia Blemker have been
implemented. These models can be used to simulate the passive and active
response of muscles and the passive response of tendons. See the user
manual for more details on these material models.

2. FEATURE: strain output in log file
-------------------------------------
The element strain values can now also be outputted to the log file. See
the user manual for more details.

2. BUGFIX: plot file data for triangular shells
-----------------------------------------------
The stress and strain data for triangular shells was not stored correctly
in the plot database file. This has now been fixed. 

3. BUGFIX: divergence logic
---------------------------
A small bug was fixed in the divergence logic.

===========================================================================

 R E L E A S E   1.1.3

===========================================================================

1.FEATURE: iterative linear solver
----------------------------------
An iterative linear solver has been added that can be used as an alternative
to the default skyline solver to solve the system of linear equations Ku = R.
This iterative solver is attractive for large problems that either require too 
much memory or take too long to factorize the stiffness matrix. 
To select the solver, add the following section to the Control section
of your input file.

<linear_solver type="conjugate gradient">
	<tolerance>0.01</tolerance>
	<max_iterations>50</max_iterations>
	<print_level>0</print_level>
</linear_solver>

This solver takes three parameters (shown in the example above with their 
default values). The "tolerance" parameter specifies the convergence tolerance
on the residual norm of the linear system, where the residual is defined here
as r = R - Ku. The "max_iterations" parameter sets a maximum on the number of
iterations. The "print_level" parameter sets the amount of output that the
linear solver will generate. When set to zero (=default), no output is generated.
When set to a non-zero value, convergence information is printed only to the
screen (not to the logfile!).


2.FEATURE: shell strains
------------------------
FEBio can now output shell strains. To output shell strains set the "shell_strain"
parameter to 1 in the "plotfile" section of the "Output" section of your input file.

<Output>
	<plotfile>
		<shell_strain>1</shell_strain>
	</plotfile>
</Output>

The default value for this flag is 0 (=off).


3.FEATURE: progress indicator in window title
---------------------------------------------
A progress indicator is displayed in the title of the console window
or terminal for linux/mac). That way the user can monitor the progress of a
run even when the window is minimized.


4.BUG FIX: NIKE3D wedge element
-------------------------------
FEBio read the node numbering incorrectly of the wedge element for NIKE3D
input files.


5.BUG FIX: linesearch
---------------------
A small bug was fixed in the linesearch algorithm.


===========================================================================

 R E L E A S E   1.1.2

===========================================================================

1.FEATURE: Augmented Lagrangian for incompressible materials
----------------------------------------------------------
An augmented Lagrangian algorithm was added for the incompressibility
enforcement. To use this feature, add the following two parameters
to any material that uses a decoupled formulation (e.g. Mooney-Rivlin)

<laugon>1</laugon>
<atol>0.05</atol>

The first parameter activates the augmented Lagrangian for incompres-
sibility and the second parameter specifies the augmented Lagrangian
convergence tolerance. 


2.FEATURE: TC nonlinear orthotropic
---------------------------------
A tension-compression nonlinear orthotropic material was added. To use
this material set the material type to "TC nonlinear orthotropic" and
define the material parameters. See the user's manual for a description
of the material parameters.


3.FEATURE: contact auto penalty calculation
-----------------------------------------
The initial penalty factor can now be calculated automatically. Simply
specify the "auto" attribute and as value "on" or "off" to use the 
auto penalty estimation. Default is off.

<penalty auto="on">1</penalty>

In this case the specified value of the penalty element is a scale factor.


4.FEATURE: max condition number
-----------------------------
When the condition number of the stiffness matrix gets too large, the matrix
becomes ill-conditioned. When this happens, FEBio will reform the stiffness
matrix. The maximum value for the condition number (previously hardcoded at 1e5)
can now be set by the user. To set this value, add the "cmax" parameter to
the control section of your input file.

<cmax>1e5</cmax>


5.BUG FIX: contact with tetrahedral elements
------------------------------------------
There were a few bugs found and fixed in the contact algorithm 
for tetrahedral elements. One bug was rather severe since it could
crash FEBio.


6.BUG FIX: Dilatational stiffness
-------------------------------
A bug was found and corrected in the routine that calculates the dilational
stiffness. This may improve convergence behavior for all materials using a 
decoupled formulation. 


7.BUG FIX: auto timestepper
-------------------------
A bug was found and fixed in the auto timestepper routine. Specifically,
the max time step size was sometimes not respected when it was specified
through a loadcurve.


===========================================================================

 R E L E A S E   1.1.1

===========================================================================

1.FEATURE: data check option
----------------------------
A command line option was added that requests FEBio to check the file for 
potential problems, but without running the file. To use the feature add
the -c option on the command line. For instance, the following command line
will perform a data check on the file without running the problem.

>febio -i in.feb -c [enter]


2. FEATURE: tied contact interface
----------------------------------
A tied contact interface was added. This interface allows you to connect two
non-conforming meshes. A tied contact is defined similarly as a sliding 
contact interface, except the contact type is now set to "tied" and there 
is two pass option. For, instance the following example defines a simple
tied interface.

<contact type="tied">
	<tolerance>0.01</tolerance>
	<penalty>50</penalty>
	<surface type="master">
		<quad4>1,2,3,4</quad4>
	</surface>
	<surface type="slave">
		<quad4>5,6,7,8</quad4>
	</surface>
</contact>


3.BUG FIX: divergence logic
----------------------------
A small bug was fixed in the algorithm that determines whether an iteration
is diverging or not.


4.BUG FIX: rigid body reaction forces
--------------------------------------
When using prescribed rigid body forces, the rigid body reaction forces where
not reported accurately.


5.BUG FIX: rigid body data logging
-----------------------------------
A bug was fixed in the data logging feature of rigid body data when the user
did not explicitly define the rigid material id(s).


6.BUG FIX: bug in bandwidth optimizer
-------------------------------------
A bug was found and fixed in the bandwidth optimization algorithm.


7.BUG FIX: St.Venant-Kirchhoff
------------------------------
The St.Venant-Kirchhoff material was not recognized by FEBio when spelled
like this. Instead it expected a hyphen instead of the first dot. This is 
fixed and now the correct type for this material is spelled as follow:
"St.Venant-Kirchhoff"


===========================================================================

 R E L E A S E   1.1.0

===========================================================================

1. FEATURE: Bandwidth reduction algorithm
-----------------------------------------
A bandwidth reduction algorithm has been implemented for the default Skyline
solver. This algorithm attempts to minimize the bandwidth of the global
stiffness matrix. The result can be a dramatic reduction in memory
requirements and computing time. To use the algorithm, simply add the
following option to the Control section of your input file.

<optimize_bw>1</optimize_bw>


2. FEATURE: Data Logging
-----------------------
The user can request specific data to be stored in the log file. Data can be
requested for nodes, elements and rigid bodies. To request additional data you
need to define the new Output section in your input file. In the Output
section, define a logfile section and in this section define the data that you
wish to add to the log file. 

<Output>
	<logfile>
		<node_data data="x;y;z"></node_data>
	</logfile>
</Output>


3. FEATURE: Rigid wall interface
--------------------------------
The rigid wall interface is a simplification of the sliding contact interface
where the master surface is a movable wall which can be defined implicitly
using a plane equation. No geometry needs to be defined for the wall. The
position of the wall can be controlled with a load curve.

<contact type="rigid_wall">
	<plane lc="1">0,0,1,0</plane>
	<tolerance>0.1</tolerance>
	<penalty>1</penalty>
</contact>


4. FEATURE: User fiber distribution
-----------------------------------
An option has been added for transversely isotropic materials to input a user
defined fiber distribution. The local fiber direction can now be defined on an
element per element basis. To use this option set the fiber type to "user" in
the material section. In addition, define the new ElementData section in the
Geometry section of your input file. In this section specify for each element
the local fiber direction.

<ElementData>
	<element id="1">
		<fiber>1,0,0</fiber>
	</element> 
</ElementData>


5. FEATURE: Shell formulation
-----------------------------
FEBio now supports shell elements. The current formulation allows either four
node quadrilateral or three node triangular shells with possible thickness
variations over the shell’s surface. Shells are defined in the Geometry
section of your input file. Use quad4 for a quadrilateral element and tri3
for a triangular element.

<Elements>
	<quad4 id="1" mat="1">1,2,3,4</quad4>
	<tri3 id="2" mat="1">5,6,7</tri3>
</Elements>

The shell thickness is defined in the ElementData section and needs to be
speficied for each node.

<ElementData>
	<element id="1">
		<thickness>0.02,0.02,0.02,0.02</thickness>
	</element>
</ElementData>


6. FEATURE: Restart File
------------------------
The user can define a restart input file that can be used to alter some
control parameters when restarting a previously terminated run. To request the
creation of a restart file, add the following option to the Control section.

<restart>1</restart>

This option requests the creation of a dump file that can be used to restart
the analysis. When restarting the analysis you can specify a restart input file
that may override some of the original control options.


7. FEATURE: Rigid Body Forces
-----------------------------
Prescribed forces can now be applied directly onto the rigid body’s center of
mass. To specify a force simply set the type attribute of the corresponding
rigid body’s degree of freedom to force and set the value to the force’s
magnitude. Applying a force to a rotational degree of freedom defines a torque
around the center of mass.

<material id="1" type="rigid body">
	<trans_x type="force" lc="1">5.0</trans_x>
	...
</material>

8. FEATURE: Plot Level
----------------------
The user has now more control over when exactly a state is added to the
plotfile. Choices include: never, all major iterations, all minor iterations,
only for must points, .... Enter the following option in the Control section of
your input file to set the plot level.

<plot_level>PLOT_MAJOR_ITRS</plot_level>


9. BUG FIX: Improved Contact
----------------------------
Several bugs were found and fixed in the contact implementation. Improved
convergence can be expected for contact analyses.


10. BUG FIX: Linesearch
-----------------------
A bug was found and fixed in the linesearch algorithm. Problems that depend
heavily on a linesearch will show improved convergence.


11. BUG FIX: Memory allocation failure
--------------------------------------
FEBio now prints an error message when memory allocation for the global
stiffness matrix fails.
