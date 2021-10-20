Glider_ADCP_Toolbox
---------------------------

The glider-ADCP toolbox is a set of MATLAB (version >= 2018b) scripts and functions developed at University of Perpignan to manage the data collected by autonomous underwater gliders. They cover the main stages of the data management process only in delayed time mode: data concatenation, data processing, and generation of data products and figures. 

A focus is maide on the processing of Acoustic Doppler Current Profile (ADCP) data, which is a complex sensor and where a need for processing exists.

FEATURES
---------------------------

The following features are already implemented in the glider-ADCP toolbox (privat folder):

	* Two main scripts to perform delayed time data processing:
		- Glider_ADCP_define_param.m
		- Glider_ADCP_main_program.m
	* Support Slocum Glider G1 model.
	* Concatenation of glider data segments.
	* Data processing, including:
		- unit conversions
		- Quality control with optional parameter estimation
		- parameters derivations, such as absolute velocities
	* Data interpolation over a 2D grid and data mapping.
	* Generation of figure products.
	* Generation of matrix data products (.mat)
	* Configuration of every processing stage using configuration file (Glider_ADCP_define_param.m)
		- parameters for conversions, corrections, derivations, interpolations and filters
		- parameters for data gridding
		- customizable figure outputs
		
The following features are planned or in development:
	* raw data retrieval from glider, storage and load options.
	* Implementation of inverse method for absolute velocities estimation.
	
DOCUMENTATION
--------------------------

The toolbox is documented through a Chapter II of the Mathieu Gentil's thesis "Glider-ADCP toolbox: A MATLAB toolbox for processing active acoustic data onto underwater glider". The documentation describes all the steps of the treatment with figure outputs as well as the parameters that can be control by the user. In addition, the limitations and recommendations of the toolbox are described in the documentation, in order to provide a framework to process acoustic glider data (e.g., to as best possible avoid “black box” usability of the tools that leads to the user not knowing the way the data has been handled). 

The documentation folder contains:
	* Documentation: "Glider-ADCP toolbox: A MATLAB toolbox for processing active acoustic data onto underwater glider".
	* The glider toolbox diagram: from deployment database to figures.
	
GLIDER DATASET
-------------------------

A dataset of an acoustic-glider deployment is provided to run the toolbox (glider_data folder). This is a deployment carried out in the NW Mediterranean in 2017, as part of the Matugli project.

	* The raw glider data are stored in the "Raw" folder and includes: .dbd files (navigation data), .ebd files (sciences data), .PD0 files (acoustic data), .mat file (Suspended particulates matter calibration).
	* The processed folder includes: figure outputs and the 5 levels (L0 to L4) of matlab files (.mat) outputs (see the documentation).
	
RUN
-------------------------

	- Run the Script Glider_ADCP_main_program.m;
	- Select the path of the toolbox folder in the dialog box;
	- Follow the instructions in the dialog box and messages dispalyed in the command window.

	
NOTE
-------------------------

It is a toolbox in the development stage. The codes are "home-made" and deserve to be optimised. Furthermore, the addition of external data would be interesting to compare the methods used to derive absolute velocities from glider-ADCP.

Contributions and criticism are welcome. Feel free to contribute ! :)


Authors
--------------------------

This toolbox was mainly developed in the framework of Mathieu Gentil's thesis. But it is the result of many contributions from Xavier Durrieu de Madron, Pierre Cauchy, François Bourrin and especially Gaël Many.



	


