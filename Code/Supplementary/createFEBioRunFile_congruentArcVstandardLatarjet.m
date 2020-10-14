function [glenoidMeshOutput,headMeshOutput] = createFEBioRunFile_congruentArcVstandardLatarjet(glenoidMesh,febioFebFileNamePart,...
    glenoidRotation,humerusCS,shapes,simType,generatePlots)
    

    %%%% Potentially need to mesh the head sphere volumetrically and to
    %%%% apply the appropriate element types for this and the surface...
    
    %% This function serves to import in and create the base surface system for
    %  running the FEA analysis of the humeral head against the glenoid.
    %
    %  Inputs:
    %
    %   glenoidMesh         structure containing surface and volumetric aspects
    %                       of the current glenoid mesh
    %   febioFileNamePart   name to store FEBio file outputs under
    %   glenoidRotation     amount to rotate the glenoid counter-clockwise
    %                       to change the direction of translation
    %   humerusCS           structure containing details of humerus coordinate
    %                       system so that the glenoid can be rotated
    %   shapes              structure containing imported humeral head
    %                       sphere shape
    %   simType             string of 'translation' or 'force' for the two
    %                       different driven simulations
    %   generatePlots       flag whether to generate figures from the
    %                       processing throughout the function (default = false)
    %   
    %  Outputs
    %
    %   None                
    %
    %  This function, like parts of the main code, uses elements of the GIBBON
    %  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
    %  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).
    %
    %  References
    %
    %  Moroder et al. (2019). Challenging the current concept of critical
    %  glenoid bone loss in shoulder instability: Does the size measurement
    %  really tell all? Am J Sports Med, 47: 688-694.
    %
    %  Walia et al. (2013). Theoretical model of the effect of combined
    %  glenohumeral bone defects on anterior shoulder instability: A finite
    %  element approach. J Orthop Res, 31: 601-607.
    %
    %  Walia et al. (2015). Influence of combined Hill-Sachs and bony
    %  Bankart defects on range of motion in anterior instability of the
    %  shoulder in a finite element model. Arthroscopy, 31: 2119-2127.

    %% Check inputs

    %Need all inputs
    if nargin < 5
        error('Need all function inputs except for generatePlots')
    end
    if nargin < 6
        generatePlots = false;
    end
    
    %% Set-up
    
    %Set the starting deep glenoid origin and humeral head points
    originPt = [0 0 0];
    headPt = humerusCS.Xc(1:3);
    
    %Set the starting compression and translation vectors that will be
    %updated based on the rotation of the objects. The starting compression
    %is directly along the Z-axis, while the translation relates to the
    %different translation clock face inputs
    compressionVec = [0 0 1];
    
% % %     %Create a string for the coordinate system labels
% % %     csLabels = [{'Xc'}, {'Yc'}, {'Zc'}];
% % %     csColours = [{'r'}; {'g'}; {'b'}];
    
    %Extract the nodes, elements and surfaces from the mesh structures
    glenoidVolV = glenoidMesh.glenoidVolV;
% % %     headVolV = headMesh.headVolV;
    glenoidVolE = glenoidMesh.glenoidVolE;
% % %     headVolE = headMesh.headVolE;
    glenoidF = glenoidMesh.glenoidF;
% % %     headF = headMesh.headVolFb;
    
    %Calculate point spacing if not input
    pointSpacing = mean(patchEdgeLengths(glenoidF,glenoidVolV));        
    
    %Create a list humeral landmarks
% % %     humerusLandmarks = [{'GHJC'},{'EL'},{'EM'},{'EJC'}];
% % %     humerusLandmarks = [{'GH'},{'EL'},{'EM'},{'EJC'}];
% % %     
% % %     %Create a list of scapula landmarks
% % %     scapulaLandmarks = [{'AA'},{'AI'},{'TS'},{'DeepGlenoid'},{'SGT'},{'IGT'},{'AntGlenoid'},{'PostGlenoid'},{'SupGlenoid'},{'InfGlenoid'}];
    
    %% Rotate glenoid head considering clock face translation
    
    %The glenoid is now aligned that the 3 o clock direction represents the
    %positive X-axis (i.e. absolute anterior translation). The
    %glenoidRotation input argument specifies the amount to
    %counter-clockwise rotate the glenoid so that the translational X-axis
    %shifts in the antero-inferior direction. 
    
    %Check for non-zero value
    if glenoidRotation ~= 0

        %Create the rotation transform about the pure Z-axis
        zRot = createRotationOz(deg2rad(glenoidRotation));

        %Apply the rotation to the surfaces
        for pp = 1:length(glenoidVolV)
            glenoidVolV(pp,:) = transformPoint3d(glenoidVolV(pp,:),zRot);    
        end
        clear pp

    end
    
    %% Create and position a humeral head sphere above glenoid centre
    
    %A simple sphere is created here to represent the humeral head. It is
    %placed at a point so that it's edge is 0.5mm away from the glenoid
    %origin point. Given that no rotation or translation has occurred with
    %respect to the original head point, no additional movement is required
    %for the starting point.
    
    %Create the sphere centering around the specified humeral head point
    %and with the specified head radius.
    [headVolE,headVolV,~] = geoSphere(3,shapes.headSphere.radius);
    
    %Shift the created sphere to the starting head point
    headVolV(:,3) = headVolV(:,3) + (headPt(3) - 9.5);
    
    %Create a volumetric mesh of the sphere
    headInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
    headInputStruct.Faces = headVolE; %Boundary faces
    headInputStruct.Nodes = headVolV; %Nodes of boundary
    headInputStruct.regionPoints= getInnerPoint(headVolE,headVolV); %Interior points for regions
    headInputStruct.holePoints = []; %Interior points for holes
    headInputStruct.regionA = tetVolMeanEst(headVolE,headVolV); %Desired tetrahedral volume for each region

    %Mesh model
    [headMesh] = runTetGen(headInputStruct); %Run tetGen

    %Get outputs of mesh structure
    headVolE = headMesh.elements; %The elements
    headVolV = headMesh.nodes; %The vertices or nodes
    headVolFb = headMesh.facesBoundary; %The boundary faces
    headVolCb = headMesh.boundaryMarker; %The boundary markers
    
    %Need to invert faces for some reason
    headVolFb = fliplr(headVolFb);
       
    %% Set force vector values
    
    %Set original force values
    compressionForce = 100;
    translationForce = 100;
    
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
    % that represents rigid bone. These values are also only used in the
    % contact model calculations, given that rigid materials are used.

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

    febio_spec.Geometry.Elements{2}.ATTR.type = 'tet4'; %Element type of this set
    febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material index for this set
    febio_spec.Geometry.Elements{2}.ATTR.name = 'headPart'; %Name of the element set
    febio_spec.Geometry.Elements{2}.elem.ATTR.id = size(glenoidVolE,1)+(1:1:size(headVolE,1))'; %Element id's
    febio_spec.Geometry.Elements{2}.elem.VAL = headVolE + size(glenoidVolV,1);

    %% Surface section
    
    %Having aspects on the cut flat surface as contact faces seems to slow
    %down the simulation, potentially due to the shering sort of contact
    %they experience. We can remove these from the contact surface by
    %extracting the patch norms that are directly points along the X-axis.
    %We can also remove those on the back surface too for brevity based on
    %those that are exactly -Z
    
    %%%%% TODO: sphere contacts coracoid as well with at least a 2 o clock
    %%%%% translation, can't have those as contact surfaces...
    
    %Extract the glenoid face normals
    glenoidNormals = patchNormal(glenoidF,glenoidVolV);
    
    %Identify those that are either directly pointing along the +X or -Z
    notFaces = glenoidNormals(:,1) == 1 | glenoidNormals(:,3) == -1;
    
    %Get the glenoid surface faces that we want to keep
    keepFaces = glenoidF(~notFaces,:);
    
    % Plotting surface models and contact faces
    if generatePlots
        cFigure; hold on;
        title('Contact sets and normal directions');
        gpatch(glenoidF,glenoidVolV,'kw','none',0.3); 
        hl(1) = gpatch(glenoidF,glenoidVolV,'g','k',1); 
        patchNormPlot(keepFaces,glenoidVolV);
        hl(2) = gpatch(headVolFb,headVolV,'b','k',1);
        patchNormPlot(headVolFb,headVolV);
        legend(hl,{'Master','Slave'});
        axisGeom;
        camlight headlight;
    end
    
    %Add surfaces
    febio_spec.Geometry.Surface{1}.ATTR.name = 'headToGlenoid_master';
    febio_spec.Geometry.Surface{1}.tri3.ATTR.lid = (1:1:size(keepFaces,1))';
    febio_spec.Geometry.Surface{1}.tri3.VAL = keepFaces;
    
    febio_spec.Geometry.Surface{2}.ATTR.name = 'headToGlenoid_slave';
    febio_spec.Geometry.Surface{2}.tri3.ATTR.lid = (1:1:size(headVolFb,1))';
    febio_spec.Geometry.Surface{2}.tri3.VAL = headVolFb + size(glenoidVolV,1);

    %Surface pairs
    febio_spec.Geometry.SurfacePair{1}.ATTR.name = 'headToGlenoid';
    febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface = febio_spec.Geometry.Surface{1}.ATTR.name;
    febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface = febio_spec.Geometry.Surface{2}.ATTR.name;

    %% Rigid body boundaries for initial step (and throughout)

    %Glenoid remains rigidly fixed across entire simulation
    febio_spec.Boundary.rigid_body{1}.ATTR.mat = 1;
    febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
    febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
    febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
    febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
    febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
    febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

    %% Contact section
    
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
    
    %Set the load curves to apply to the compression and translational
    %force components. Note that the same load curve can be used for each
    %of these different forces as we want them to remain consistent across
    %the XYZ components.
    
    %We set the translational load curve reach similar force stability
    %ratios generated by Moroder et al. (2019). In essence, the design of
    %the translational step combined with the translational load curve here
    %results in the translational load increasing by 1N each step,
    %starting at 10N, in line with Moroder et al. (2019). We set the load
    %curve here to progressively reach the value of 100N over a 10 second
    %period (i.e. 10N.sec) starting from the 2 second point --- but limit
    %the number of steps later so that the simulation doesn't go forever
    %and unnecessarily high translational forces are applied. The time
    %frame for steps was based on Moroder et al. (2019) stability ratios
    %they derived, however some experimentation to get it right was also
    %required.
    %
    %This is slightly different to Moroder et al. (2019) which increased by
    %0.1N each step, however this is pretty cumbersome and results in super
    %long simulation times --- and only increases the resolution of the
    %stability ratio by one decimal place. It also results in the head not
    %reaching the glenoid edge smoothly (i.e. it progresses to the edge and
    %then rolls back) --- which is computationally expensive from the model
    %perspective.
    %
    %The translation driven approach causes 45mm displacement at 0.5mm
    %intervals, and hence takes 9 seconds to reach it's maximum of 30. This
    %is why we have a load curve of 2-11 seconds.
    
    %Compression load curve
    febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
    febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
    febio_spec.LoadData.loadcurve{1}.point.VAL = [0 0; 2 1];

    %Translation load curve
    %Driven by displacement or force
    if strcmp(simType,'translation')
        febio_spec.LoadData.loadcurve{2}.ATTR.id = 2;
        febio_spec.LoadData.loadcurve{2}.ATTR.type = 'linear';
        febio_spec.LoadData.loadcurve{2}.point.VAL = [2 0; 11 1];
    elseif strcmp(simType,'force')
        febio_spec.LoadData.loadcurve{2}.ATTR.id = 2;
        febio_spec.LoadData.loadcurve{2}.ATTR.type = 'linear';
        febio_spec.LoadData.loadcurve{2}.point.VAL = [2 0; 2.1 0.1; 11.1 1];
    end    

    %% Output section

    %Clear out original
    febio_spec = rmfield(febio_spec,'Output');

    %Repopulate
    febio_spec.Output.plotfile.ATTR.type = 'febio';
    febio_spec.Output.plotfile.var{1}.ATTR.type = 'displacement';

    %% Step section
    
    %Remove the control section to use step control
    febio_spec = rmfield(febio_spec,'Control');

    %Add first step

    %Control section
    %Compression step. This applies the compression force linearly
    %increasing from 0 to the desired 100N over 3 seconds at 0.1 step
    %intervals. This is necessary to softly position the head in the deep
    %glenoid point.
    febio_spec.Step{1}.ATTR.name = 'compress';
    febio_spec.Step{1}.Control.time_steps = 20;
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
    %Given the compression is occurring along the Z-axis, we can limit the
    %motion of the head in this step to only occur along this axis, so that
    %it gets seated nicely into the glenoid.
    febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 2;
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
    febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';

    %Load section
    %Add the compression load to the head elements, purely Z-axis force
    febio_spec.Step{1}.Loads.body_load{1}.ATTR.type = 'const';
    febio_spec.Step{1}.Loads.body_load{1}.ATTR.elem_set = 'headPart';
    febio_spec.Step{1}.Loads.body_load{1}.z.VAL = compressionForce;
    febio_spec.Step{1}.Loads.body_load{1}.z.ATTR.lc = 1;

    %Add second step

    %Control section
    %Translation step. This applies the translational force at 1N
    %intervals at each 0.1s step size. We set this to 50 steps in order to
    %reach a force of 60N, as the stability ratios provided in Moroder et
    %al. (2019) did not appear to reach this (not much above 50%).
    
    %Second step is either displacement vs. force control
    if strcmp(simType,'translation')
        
        %The specified steps here aim to move the head 45mm. It moves 0.5mm
        %at each step, hence we need 90 steps
        febio_spec.Step{2}.ATTR.name = 'translate';
        febio_spec.Step{2}.Control.time_steps = 90;
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
        %Given that compression is occurring along the Z-axis, and translation
        %is occurring purely along the X-axis --- we can just allow movements
        %along these directions.
        febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat = 2;
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'y';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'Rx';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Ry';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Rz';
        febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'x';
        febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 2;
        febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL = 45;
        
    elseif strcmp(simType,'force')    
        
        febio_spec.Step{2}.ATTR.name = 'translate';
        febio_spec.Step{2}.Control.time_steps = 50;
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
        %Given that compression is occurring along the Z-axis, and translation
        %is occurring purely along the X-axis --- we can just allow movements
        %along these directions.
        febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat = 2;
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'y';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'Rx';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Ry';
        febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Rz';

        %Load section
        %Add the translational load to the head elements, purely X-axis force 
        %but needs to be inverted to translate anteriorly along the glenoid
        febio_spec.Step{2}.Loads.body_load{1}.ATTR.type = 'const';
        febio_spec.Step{2}.Loads.body_load{1}.ATTR.elem_set = 'headPart';
        febio_spec.Step{2}.Loads.body_load{1}.x.VAL = translationForce*-1;
        febio_spec.Step{2}.Loads.body_load{1}.x.ATTR.lc = 2;
        
    end

    %% Output section 
    
    %Set filename variable
    fileName = [febioFebFileNamePart,'_',simType];
    
    %Defining file names
    febioFebFileName = [fileName,'.feb']; %FEB file name
    febioLogFileName = [fileName,'.txt']; %FEBio log file name
    febioLogFileName_dispHead = [fileName,'_dispHead_out.txt']; %Log file name for exporting displacement
    % febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

    %Log file
    %The displacement log file simply grabs the head nodes, as we can
    %assume that the glenoid remains fixed when animating later. The head
    %is also rigid, so we can reduce file size and later processing time by
    %simply grabbing a single node for head displacement.
    febio_spec.Output.logfile.ATTR.file = febioLogFileName;
    febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_dispHead;
    febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
    febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
    febio_spec.Output.logfile.node_data{1}.VAL = size(glenoidVolV,1)+1;%(1:size(headVolV,1));
    
    %% Save to file 
    
    %Save the .feb setup to .xml style file
    febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
    
    %% Outputs
    
    %Send the rotated meshes back out for animation after simulation
    glenoidMeshOutput.glenoidVolV = glenoidVolV;
    glenoidMeshOutput.glenoidVolE = glenoidVolE;
    glenoidMeshOutput.glenoidF = glenoidF;
    headMeshOutput.headVolV = headVolV;
    headMeshOutput.headVolE = headVolE;
    headMeshOutput.headVolFb = headVolFb;
    
%% %%%%% ----- End of createFEBioRunFile_congruentArcVstandardLatarjet.m ----- %%%%%
end