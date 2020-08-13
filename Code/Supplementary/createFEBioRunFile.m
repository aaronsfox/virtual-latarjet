function createFEBioRunFile(glenoidMesh,headMesh,febioFebFileNamePart,[],generatePlots)

%% This function serves to import in and create the base surface system for
%  running the FEA analysis of the humeral head against the glenoid.
%
%  Inputs:
%
%   ...                 
%   ...                 ...
%
%  Outputs
%
%   ...                 ...
%
%  This function, like parts of the main code, uses elements of the GIBBON
%  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
%  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).

%% TODO: currently running sim with a neutral position (i.e. 0,0,-1) for compression etc.

    %% Define febio structure

    %Get a template with default settings 
    [febio_spec] = febioStructTemplate;

    %febio_spec version 
    febio_spec.ATTR.version = '2.5'; 

    %Module section
    febio_spec.Module.ATTR.type = 'solid'; 

    %% Material section
    
    % Rigid materials are given a Young's modulus high enough to be rigid,
    % but are also estimated off the cortical bone values used in Edwards
    % et al. (2010), Clin Biomech, 25: 372-377. Admittedly this was a study
    % on the tibia, but nonetheless it represents an estimate high enough
    % that represents rigid bone.

    %Rigid glenoid
    febio_spec.Material.material{1}.ATTR.type = 'rigid body';
    febio_spec.Material.material{1}.ATTR.name = 'rigidGlenoid';
    febio_spec.Material.material{1}.ATTR.id = 1;
    febio_spec.Material.material{1}.E = 18600; %MPa; large value for rigidity
    febio_spec.Material.material{1}.v = 0.4;
    febio_spec.Material.material{1}.density = 0.00000186; %from estimate of bone density at 1,860kg/m^3

    %Rigid head
    febio_spec.Material.material{2}.ATTR.type = 'rigid body';
    febio_spec.Material.material{2}.ATTR.name = 'rigidHead';
    febio_spec.Material.material{2}.ATTR.id = 2;
    febio_spec.Material.material{2}.E = 18600; %MPa; large value for rigidity
    febio_spec.Material.material{2}.v = 0.4;
    febio_spec.Material.material{2}.density = 0.00000186; %from estimate of bone density at 1,860kg/m^3

    %% Geometry section

    % % % %Join node sets
    % % % V = [glenoidVolV; headVolV; ]; %Combined node sets
    % % % headVolFb = headVolFb + size(glenoidVolV,1); %Fixed element indices for head surface
    
    %%%%% TODO: set mesh aspects as input structres
    glenoidVolV = glenoidMesh.glenoidVolV;
    headVolV = headMesh.headVolV;
    glenoidVolE = glenoidMesh.glenoidVolE;
    headVolE = headMesh.headVolE;
    
    %Shift head nodes to be 0.5 mm from glenoid surface
    %%%%% TODO: this will need to be adapted based on the way the head has
    %%%%% been rotated
    headVolV(:,3) = headVolV(:,3) - 9.5;
    
    %Add nodes
    febio_spec.Geometry.Nodes{1}.ATTR.name = 'glenoid'; %The node set name
    febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(glenoidVolV,1))'; %The node id's
    febio_spec.Geometry.Nodes{1}.node.VAL = glenoidVolV; %The nodel coordinates

    febio_spec.Geometry.Nodes{2}.ATTR.name = 'head'; %The node set name
    febio_spec.Geometry.Nodes{2}.node.ATTR.id = size(glenoidVolV,1)+(1:size(headVolV,1))'; %The node id's
    febio_spec.Geometry.Nodes{2}.node.VAL = headVolV; %The nodel coordinates

    %Add elements
    febio_spec.Geometry.Elements{1}.ATTR.type = 'tet4'; %Element type of this set
    febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
    febio_spec.Geometry.Elements{1}.ATTR.name = 'glenoidPart'; %Name of the element set
    febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(glenoidVolE,1))'; %Element id's
    febio_spec.Geometry.Elements{1}.elem.VAL = glenoidVolE;

    febio_spec.Geometry.Elements{2}.ATTR.type = 'hex8'; %Element type of this set
    febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material index for this set
    febio_spec.Geometry.Elements{2}.ATTR.name = 'headPart'; %Name of the element set
    febio_spec.Geometry.Elements{2}.elem.ATTR.id = size(glenoidVolE,1)+(1:1:size(headVolE,1))'; %Element id's
    febio_spec.Geometry.Elements{2}.elem.VAL = headVolE + size(glenoidVolV,1);

    %% Surface section
    
    %%%%% TODO: set mesh aspects as input structres
    glenoidF = glenoidMesh.glenoidF;
    headF = headMesh.headVolFb;
    
    %%%%% TODO: hemisphere does seem kind of 'un-centred' but is based on
    %%%%% alignment of humerus and scapula axes?????
    
    % Plotting surface models and contact faces
    if generatePlots
        cFigure; hold on;
        title('Contact sets and normal directions');
        gpatch(glenoidF,glenoidVolV,'kw','none',0.3); 
        hl(1) = gpatch(glenoidF,glenoidVolV,'g','k',1); 
        patchNormPlot(glenoidF,glenoidVolV);
        hl(2) = gpatch(headF,headVolV,'b','k',1);
        patchNormPlot(headF,headVolV);
        legend(hl,{'Master','Slave'});
        axisGeom;
        camlight headlight;
    end
    
    %Add surfaces
    febio_spec.Geometry.Surface{1}.ATTR.name = 'headToGlenoid_master';
    febio_spec.Geometry.Surface{1}.tri3.ATTR.lid = (1:1:size(glenoidF,1))';
    febio_spec.Geometry.Surface{1}.tri3.VAL = glenoidF;
    
    febio_spec.Geometry.Surface{2}.ATTR.name = 'headToGlenoid_slave';
    febio_spec.Geometry.Surface{2}.quad4.ATTR.lid = (1:1:size(headF,1))';
    febio_spec.Geometry.Surface{2}.quad4.VAL = headF + size(glenoidVolV,1);

    %Surface pairs
    febio_spec.Geometry.SurfacePair{1}.ATTR.name = 'headToGlenoid';
    febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface = febio_spec.Geometry.Surface{1}.ATTR.name;
    febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface = febio_spec.Geometry.Surface{2}.ATTR.name;

    %% Rigid body boundaries for initial step (and throughout)

    %Rigid glenoid
    febio_spec.Boundary.rigid_body{1}.ATTR.mat = 1;
    febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
    febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
    febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
    febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
    febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
    febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

    %% Contact section
    
    %%%%% TODO: refine idea of inputting pointSpacing in here
    pointSpacing = 0.1;
    
    % The contact penalty here represents the amount of penetration allowed
    % with respect to the Young's modulus of the rigid bodies. This is
    % currently set at 18,600 --- hence if we wanted to allow 1mm of
    % penetration we would set this at 18,600MPa/mm. Here we allow a half
    % mm of penetration, hence 18600*(1/0.5) is the penalty.
    
    %Contact section
    febio_spec.Contact.contact{1}.ATTR.surface_pair = febio_spec.Geometry.SurfacePair{1}.ATTR.name;
    febio_spec.Contact.contact{1}.ATTR.type = 'sliding-elastic';
    febio_spec.Contact.contact{1}.ATTR.name = 'headToGlenoid';
    febio_spec.Contact.contact{1}.two_pass = 1;
    febio_spec.Contact.contact{1}.laugon = 0;
    febio_spec.Contact.contact{1}.tolerance = 0.1;
    febio_spec.Contact.contact{1}.gaptol = 0;
    febio_spec.Contact.contact{1}.search_tol = 0.01;
    febio_spec.Contact.contact{1}.search_radius = pointSpacing/10;
    febio_spec.Contact.contact{1}.symmetric_stiffness = 0;
    febio_spec.Contact.contact{1}.auto_penalty = 0;
    febio_spec.Contact.contact{1}.penalty = 18600*(1/0.5); %25,000 worked in FEBio Studio
    febio_spec.Contact.contact{1}.fric_coeff = 0;

    %% Load curves

    %%%%% TODO: refine these, don't need them all, but just copying what comes
    %%%%% with FEBio studio output --- also edit to appropriately apply force
    %%%%% steps...
    
    %%%%% TODO: load curves will need to be altered to firstly represent
    %%%%% different translation values (i.e. not just pure X-translation),
    %%%%% as well as the increasing magnitude of the load so that it
    %%%%% increases by a certain amount each time step
    
    %%%%% It's likely we will need load curves for every x,y and z
    %%%%% direction, and just set to zero in certain cases...
    
    %LoadData section
    
    %Compression load curve
    febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
    febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
    febio_spec.LoadData.loadcurve{1}.point.VAL = [0 0; 3 1];

    %Translation load curve
    febio_spec.LoadData.loadcurve{2}.ATTR.id = 2;
    febio_spec.LoadData.loadcurve{2}.ATTR.type = 'linear';
    febio_spec.LoadData.loadcurve{2}.point.VAL = [2 0; 5 1];

    %% Output section

    %Clear out original
    febio_spec = rmfield(febio_spec,'Output');

    %Repopulate
    febio_spec.Output.plotfile.ATTR.type = 'febio';
    febio_spec.Output.plotfile.var{1}.ATTR.type = 'displacement';

    %% Step section
    
    %%%%% TODO: refine step section to be more appropriate for the force
    %%%%% control being applied
    
    %Remove the control section to use step control
    febio_spec = rmfield(febio_spec,'Control');

    %Add first step

    %Control section
    %Compression step
    febio_spec.Step{1}.ATTR.name = 'compress';
    febio_spec.Step{1}.Control.time_steps = 30;
    febio_spec.Step{1}.Control.step_size = 0.1;
    febio_spec.Step{1}.Control.max_refs = 25;
    febio_spec.Step{1}.Control.max_ups = 0;
    febio_spec.Step{1}.Control.diverge_reform = 1;
    febio_spec.Step{1}.Control.reform_each_time_step = 1;
    febio_spec.Step{1}.Control.dtol = 0.001;
    febio_spec.Step{1}.Control.etol = 0.01;
    febio_spec.Step{1}.Control.rtol = 0;
    febio_spec.Step{1}.Control.lstol = 0.9;
    febio_spec.Step{1}.Control.min_residual = 1e-20;
    febio_spec.Step{1}.Control.qnmethod = 0;
    febio_spec.Step{1}.Control.rhoi = 0;
    febio_spec.Step{1}.Control.time_stepper.dtmin = 0.01;
    febio_spec.Step{1}.Control.time_stepper.dtmax = 0.1;
    febio_spec.Step{1}.Control.time_stepper.max_retries = 10;
    febio_spec.Step{1}.Control.time_stepper.opt_iter = 10;
    febio_spec.Step{1}.Control.analysis.ATTR.type = 'dynamic';
    febio_spec.Step{1}.Control.symmetric_stiffness = 0;

    %Boundary section
    %Only allow the compression motion in first step
    febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 2;
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';

    %Load section
    %Compression load
    febio_spec.Step{1}.Loads.body_load{1}.ATTR.type = 'const';
    febio_spec.Step{1}.Loads.body_load{1}.ATTR.elem_set = 'headPart';
    febio_spec.Step{1}.Loads.body_load{1}.z.VAL = 100;
    febio_spec.Step{1}.Loads.body_load{1}.z.ATTR.lc = 1;

    %Add second step

    %Control section
    %Translation step
    febio_spec.Step{2}.ATTR.name = 'translate';
    febio_spec.Step{2}.Control.time_steps = 20;
    febio_spec.Step{2}.Control.step_size = 0.1;
    febio_spec.Step{2}.Control.max_refs = 25;
    febio_spec.Step{2}.Control.max_ups = 0;
    febio_spec.Step{2}.Control.diverge_reform = 1;
    febio_spec.Step{2}.Control.reform_each_time_step = 1;
    febio_spec.Step{2}.Control.dtol = 0.001;
    febio_spec.Step{2}.Control.etol = 0.01;
    febio_spec.Step{2}.Control.rtol = 0;
    febio_spec.Step{2}.Control.lstol = 0.9;
    febio_spec.Step{2}.Control.min_residual = 1e-20;
    febio_spec.Step{2}.Control.qnmethod = 0;
    febio_spec.Step{2}.Control.rhoi = 0;
    febio_spec.Step{2}.Control.time_stepper.dtmin = 0.01;
    febio_spec.Step{2}.Control.time_stepper.dtmax = 0.1;
    febio_spec.Step{2}.Control.time_stepper.max_retries = 10;
    febio_spec.Step{2}.Control.time_stepper.opt_iter = 10;
    febio_spec.Step{2}.Control.analysis.ATTR.type = 'dynamic';
    febio_spec.Step{2}.Control.symmetric_stiffness = 0;

    %Boundary section
    %Allow compression and translation movement
    febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat = 2;
    febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'y';
    febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'Rx';
    febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Ry';
    febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Rz';

    %Load section
    %Translational load
    febio_spec.Step{2}.Loads.body_load{1}.ATTR.type = 'const';
    febio_spec.Step{2}.Loads.body_load{1}.ATTR.elem_set = 'headPart';
    febio_spec.Step{2}.Loads.body_load{1}.x.VAL = -100;
    febio_spec.Step{2}.Loads.body_load{1}.x.ATTR.lc = 2;

    %% Output section 
    
    %%%%% TODO: set febioFileNamePart as function input
    febioFebFileNamePart = 'baselineSim';
    
    %Defining file names
    febioFebFileName = [febioFebFileNamePart,'.feb']; %FEB file name
    febioLogFileName = [febioFebFileNamePart,'.txt']; %FEBio log file name
    febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
    % febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

    %Log file
    febio_spec.Output.logfile.ATTR.file = febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
    febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
    %%% gets all nodes, which may be useful for animating
    febio_spec.Output.logfile.node_data{1}.VAL = 1:size([glenoidVolV;headVolV],1);

    %% Save to file 

    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode