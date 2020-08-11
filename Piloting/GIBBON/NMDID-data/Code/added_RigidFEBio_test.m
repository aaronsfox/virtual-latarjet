
%Get a template with default settings 
[febio_spec] = febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version = '2.5'; 

%Module section
febio_spec.Module.ATTR.type = 'solid'; 

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type = 'dynamic';
stepStruct.Control.time_steps = 10;
stepStruct.Control.step_size = 1/10;
stepStruct.Control.time_stepper.dtmin = (1/10)/100;
stepStruct.Control.time_stepper.dtmax = 1/10; 
stepStruct.Control.time_stepper.max_retries = 5;
stepStruct.Control.time_stepper.opt_iter = 10;
stepStruct.Control.max_refs = 25;
stepStruct.Control.max_ups = 0;
stepStruct.Control.symmetric_stiffness = 0;

%Add template based default settings to proposed control section
[stepStruct.Control] = structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec = rmfield(febio_spec,'Control'); 

%Step specific control section
febio_spec.Step{1}.Control = stepStruct.Control;
febio_spec.Step{1}.ATTR.id = 1;
% % % febio_spec.Step{2}.Control = stepStruct.Control;
% % % febio_spec.Step{2}.ATTR.id = 2;
    
%Material section
%%%%% TODO: update from neoHookean material if appropriate

% % % %Glenoid material
% % % febio_spec.Material.material{1}.ATTR.type = 'neo-Hookean';
% % % febio_spec.Material.material{1}.ATTR.name = 'elasticGlenoid';
% % % febio_spec.Material.material{1}.ATTR.id = 1;
% % % febio_spec.Material.material{1}.E = 10; %10*1e+6; %Youngs modulus; MPa to Pa conversion for mm to m
% % % febio_spec.Material.material{1}.v = 0.4; %Poisson ratio
% % % febio_spec.Material.material{1}.density = 1850*1000; %1850 kg/m3 of bone...I think

% % % %Glenoid material
% % % febio_spec.Material.material{1}.ATTR.type = 'incomp neo-Hookean';
% % % febio_spec.Material.material{1}.ATTR.name = 'elasticGlenoid';
% % % febio_spec.Material.material{1}.ATTR.id = 1;
% % % febio_spec.Material.material{1}.G = 1000; %10*1e+6; %Youngs modulus; MPa to Pa conversion for mm to m
% % % febio_spec.Material.material{1}.k = 5000; %Poisson ratio
% % % % febio_spec.Material.material{1}.density = 1850*1000; %1850 kg/m3 of bone...I think
% % % febio_spec.Material.material{1}.ATTR.type='neo-Hookean';
% % % febio_spec.Material.material{1}.ATTR.id=1;
% % % febio_spec.Material.material{1}.E=17000;
% % % febio_spec.Material.material{1}.v=0.40;

%%%%% standard neo-Hookean is still too squishy for bone--- and also
%%%%% crashes with the same simulation parameters as the Ogden material

% % % % % % %Sphere sliding material
% % % febio_spec.Material.material{1}.ATTR.type='Ogden';
% % % febio_spec.Material.material{1}.ATTR.id=1;
% % % febio_spec.Material.material{1}.c1=1e-3;
% % % febio_spec.Material.material{1}.m1=8;
% % % febio_spec.Material.material{1}.c2=1e-3;
% % % febio_spec.Material.material{1}.m2=-8;
% % % febio_spec.Material.material{1}.k=1e-3*1e2;

%Glenoid material
febio_spec.Material.material{1}.ATTR.type = 'rigid body';
febio_spec.Material.material{1}.ATTR.name = 'rigidGlenoid';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.E = 1.8e+10; %only really needed for calculating contact penalty?
febio_spec.Material.material{1}.v = 0.4; %only really needed for calculating contact penalty?
febio_spec.Material.material{1}.density = 1850; %may need to change
febio_spec.Material.material{1}.center_of_mass = mean(glenoidV)/1000;
% febio_spec.Material.material{2}.center_of_mass = center_of_mass;

%Head material
febio_spec.Material.material{2}.ATTR.type = 'rigid body';
febio_spec.Material.material{2}.ATTR.name = 'rigidHead';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.E = 1.8e+10; %only really needed for calculating contact penalty?
febio_spec.Material.material{2}.v = 0.4; %only really needed for calculating contact penalty?
febio_spec.Material.material{2}.density = 1850;
febio_spec.Material.material{2}.center_of_mass = mean(headV)/1000;

%Geometry section

% % % V = V / 1000;

% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name = 'nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = V; %The nodel coordinates

% -> Elements

%Glenoid
febio_spec.Geometry.Elements{1}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set 
febio_spec.Geometry.Elements{1}.ATTR.name = 'glenoid'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(glenoidVolE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = glenoidVolE;

%Head
febio_spec.Geometry.Elements{2}.ATTR.type = 'tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material index for this set 
febio_spec.Geometry.Elements{2}.ATTR.name = 'Sphere'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id = size(glenoidVolE,1)+(1:1:size(headF,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL = headF;

% % % % -> NodeSets
% % % %Nodes to support back of glenoid
% % % febio_spec.Geometry.NodeSet{1}.ATTR.name = 'bcSupportList';
% % % febio_spec.Geometry.NodeSet{1}.node.ATTR.id = bcSupportList(:);

% -> Surfaces

%Head surface
febio_spec.Geometry.Surface{1}.ATTR.name = 'contactHead';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid = (1:1:size(F_contact_master,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL = F_contact_master;

%Glenoid surface
febio_spec.Geometry.Surface{2}.ATTR.name = 'contactGlenoid';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid = (1:1:size(F_contact_slave,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL = F_contact_slave;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name = 'headToGlenoid';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface = febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface = febio_spec.Geometry.Surface{2}.ATTR.name;

%Boundary condition section 

% -> Fix boundary conditions

% % % %Fix the glenoid support nodes
% % % febio_spec.Boundary.fix{1}.ATTR.bc = 'x';
% % % febio_spec.Boundary.fix{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
% % % febio_spec.Boundary.fix{2}.ATTR.bc = 'y';
% % % febio_spec.Boundary.fix{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
% % % febio_spec.Boundary.fix{3}.ATTR.bc = 'z';
% % % febio_spec.Boundary.fix{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body

febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{6}.ATTR.bc = 'Rz';

%Initial displacement
%Note that distance to contact surface is 0.01m from the glenoid, so we add
%a tiny little bit to ensure contact
febio_spec.Step{1}.Boundary.rigid_body{2}.ATTR.mat = 2;
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{5}.ATTR.bc = 'Rz';
% % % febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.ATTR.bc = 'z';
% % % febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.ATTR.lc = 1;
% % % febio_spec.Step{1}.Boundary.rigid_body{2}.prescribed.VAL = 2.5; %10-0.5;

% % % %Compression to establish contact
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat = 2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'z';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL = 2.5;

% % % %Translation
% % % %%%% TODO: fix arbitrary distance and step size etc.
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat = 2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'z';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'x';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL = 0.02;

%Loads
febio_spec.Loads.body_load{1}.ATTR.type='linear';
febio_spec.Loads.body_load{1}.z.VAL=-50;
febio_spec.Loads.body_load{1}.z.ATTR.lc=1;

%Contact section
pointSpacing = mean(patchEdgeLengths(glenoidVolFb,V)); %value used earlier for basic remesh
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.two_pass=0;
febio_spec.Contact.contact{1}.laugon=0;
febio_spec.Contact.contact{1}.tolerance=0.1;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=pointSpacing/10;
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
%From FEBio discussion post: https://forums.febio.org/showthread.php?2086-Rigid-Body-Contact-Issues
%Penalty can be estimated from the Young's Modulus E of the materials
%idealising as rigid bodies, along with the maximum penetration g that
%would be willing to tolerate under the loads applied. If E = 10 GPa and g =
%0.2, then penalty should be set at 50GPa/mm (i.e. here it is 10/0.2)
febio_spec.Contact.contact{1}.auto_penalty=0;
febio_spec.Contact.contact{1}.penalty=1.8e+10 / (1/1000);
febio_spec.Contact.contact{1}.fric_coeff=0;


%LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
febio_spec.LoadData.loadcurve{1}.point.VAL = [0 0; 1 1; 2 1]; %increase third y value if wanting to keep compressing

% % % febio_spec.LoadData.loadcurve{2}.ATTR.id = 2;
% % % febio_spec.LoadData.loadcurve{2}.ATTR.type = 'linear';
% % % febio_spec.LoadData.loadcurve{2}.point.VAL = [0 0; 1 0; 2 1];

%Output section 

%%%%%%% TODO: edit these

% Defining file names
febioFebFileNamePart='sphereSliding_adaptation_rigid';
febioFebFileName=[febioFebFileNamePart,'.feb']; %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);


%% Exporting the FEBio input file
% Exporting the febio_spec structure to an FEBio input file is done using
% the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Running the FEBio analysis
% To run the analysis defined by the created FEBio input file the
% |runMonitorFEBio| function is used. The input for this function is a
% structure defining job settings e.g. the FEBio input file name. The
% optional output runFlag informs the user if the analysis was run
% succesfully. 

febioAnalysis.run_filename=febioFebFileName; %The input file name
febioAnalysis.run_logname=febioLogFileName; %The name for the log file
febioAnalysis.disp_on=1; %Display information on the command window
febioAnalysis.disp_log_on=1; %Display convergence information in the command window
febioAnalysis.runMode='external';%'internal';
febioAnalysis.t_check=0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi=1e99; %Max analysis time
febioAnalysis.maxLogCheckTime=10; %Max log file checking time

clc
[runFlag]=runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
