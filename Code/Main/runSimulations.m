
%%%%% enter notes...

%%%%% following a lot of the sphere sliding demo in GIBBON
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0007_sphere_sliding.html
%%%%% a fair bit for slicing and fixing up the hole also comes from the foot insole demo
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0050_foot_insole_01.html

%%%%% TODO: add generatePlots logical input
    
generatePlots = false;


%%%%%%%%%%% NEED TO FIX HUMERAL HEAD CARTILAGE - IT IS NOT CLOSED. THE
%%%%%%%%%%% REGIONTRIMESH3D APPROACH IS A BETTER OPTION FOR CREATING THE
%%%%%%%%%%% CARTILAGE THAN THE SWEEP LOFT (WHICH SEEMS TO WORK FOR THE
%%%%%%%%%%% GLENOID). AS PART OF THIS UPDATED PROCESS WE SHOULD ALSO INVERT
%%%%%%%%%%% THE NORMALS OF THE INTERNAL PIECE (FLIPLR APPLIED TO THE
%%%%%%%%%%% FACES). IT'S LIKELY THAT WE NEED TO TAKE A SIMILAR APPROACH TO
%%%%%%%%%%% WHAT'S DONE WITH THE SCAPULA CUTTING HERE AND FORCING THE
%%%%%%%%%%% VERTICES OF THE SEPARATE PIECES TO BE ON THE BOUNDARY (I.E. PIN
%%%%%%%%%%% THEM TO THE BOUNDARY)...

%%

%%%%% THE ABOVE POINT SEEMS TO HAVE BEEN SOLVED IN THE CREATE CARTILAGE
%%%%% FUNCTION - THE FUNCTION HAS BEEN EDITED WELL ENOUGH THAT IT COULD BE
%%%%% USED AS A SUPPLEMENTARY WITHIN THIS SCRIPT, AND RATHER THAN PASSING
%%%%% OUT TO .STL FILES, BRING THE CREATED FACES, VERTICES AND COLOURS BACK
%%%%% INTO THIS FUNCTION AS REMESHED & APPROPRIATE SURFACES...

%%%%% Would need to save the STL's in the right place if wanting to store
%%%%% them for latter use, but this wouldn't even be necessary

%%%%% Taking this apporach we could be more particular about using colours
%%%%% to separate surfaces, as these wouldn't get lost in the STL export

%% Current test code
%  Building up piece by piece to get to final simulation

%% Set-up

warning off

%Add supplementary code folder to path
addpath(genpath('..\Supplementary'));

%% Load and create relevant surface and volumetric meshes

% This step takes in relevant surface meshes exported from the 3matic
% processing to create a base system to run FEA simulations and adapt with
% bone defects to run subsequent simulations.
%
% As part of this the glenoid section is extracted from the scapula, the
% glenoid and humeral head cartilages are meshed, and the entire system is
% compiled and meshed as volumes. 
%
% Relevant files from the FEA\SurfaceMeshes directory are used in this
% process. The glenoid section is created from the 'Scapula.stl' file. The
% glenoid and humeral cartilages are created from the 'GlenoidFace.stl' and
% 'HumeralCartilageBase.stl'/'HumeralCartilageOuter.stl' files using the
% relevant supplementary functions. The idealised humeral head is meshged
% from the 'BaseHumeralHead.stl' file.

%Navigate to surface mesh directory
cd('..\..\FEA\SurfaceMeshes');

%Set-up the stl files structure to pass to the cartilage surfaces function
%This only contains the relevant files for the curvature method
stlFiles.GlenoidFace = fullfile(pwd,'GlenoidFace.stl');
stlFiles.HumeralCartilageBase = fullfile(pwd,'HumeralCartilageBase.stl');
stlFiles.HumeralCartilageOuter = fullfile(pwd,'HumeralCartilageOuter.stl');
stlFiles.HumeralCartilagePlane = fullfile(pwd,'HumeralCartilagePlane.xml');
% % % stlFiles.Scapula = fullfile(pwd,'Scapula.stl');
% % % stlFiles.Humerus = fullfile(pwd,'Humerus.stl');

%Create and obtain the cartilage surface meshes
outputSTL = createCartilageSurfaces(stlFiles,'curvature',false,pwd);

%Unpack surface mesh components

%Glenoid cartilage

%Access mesh parts
glenoidCartilageF = outputSTL.glenoidCartilage.FT;
glenoidCartilageV = outputSTL.glenoidCartilage.VT;
glenoidCartilageC = outputSTL.glenoidCartilage.CT;

%Uniformly remesh cartilage surface
[glenoidCartilageF,glenoidCartilageV] = triRemeshLabel(glenoidCartilageF,glenoidCartilageV,1.0);

%Humeral cartilage

%Access mesh parts
humeralCartilageF = outputSTL.humeralCartilage.FT;
humeralCartilageV = outputSTL.humeralCartilage.VT;
humeralCartilageC = outputSTL.humeralCartilage.CT;

%Uniformly remesh cartilage surface
[humeralCartilageF,humeralCartilageV] = triRemeshLabel(humeralCartilageF,humeralCartilageV,1.0);

%Visualise cartilages together
if generatePlots
    cFigure; hold on
    gpatch(glenoidCartilageF,glenoidCartilageV,'gw','k');
    gpatch(humeralCartilageF,humeralCartilageV,'bw','k');
    title('Cartilage Surfaces');
    axisGeom;
end

%Load in additional surface files
[headSTLstruct] = import_STL('BaseHumeralHead.stl');
[scapulaSTLstruct] = import_STL('Scapula.stl');

%Access the data from the STL structs

%Humeral head

%Get surface faces and vertices
headF = headSTLstruct.solidFaces{1}; %Faces
headV = headSTLstruct.solidVertices{1}; %Vertices
scapulaF = scapulaSTLstruct.solidFaces{1}; %Faces
scapulaV = scapulaSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[headF,headV] = mergeVertices(headF,headV);
[scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);

%Remesh humeral head sphere to reduce number of triangles
[headF,headV] = triRemeshLabel(headF,headV,3);

%Visualise all loaded in surfaces together
if generatePlots
    cFigure; hold on
    gpatch(scapulaF,scapulaV,'kw','none');
    gpatch(headF,headV,'kw','none');
    gpatch(glenoidCartilageF,glenoidCartilageV,'bw','none');
    gpatch(humeralCartilageF,humeralCartilageV,'bw','none');
    title('Complete Base System');
    camlight headlight;
    axisGeom;
end

%% Create the glenoid section from the scapula

%Slice the scapula 10mm back from the glenoid origin
%Generate settings for slicing
cutLevel = -10; %Set the cut level
snapTolerance = mean(patchEdgeLengths(scapulaF,scapulaV))/100;
n = vecnormalize([0 0 -1]); %Normal direction to plane
P = [0 0 cutLevel]; %Point on plane

%Slicing surface (note 3rd color data output is supressed)
[scapulaFc,scapulaVc,~,logicSide,scapulaEc] = triSurfSlice(scapulaF,scapulaV,[],P,n,snapTolerance);

%Visualise slice
if generatePlots
    %Plot split planes
    cFigure; subplot(1,2,1); hold on;
    hp1 = gpatch(scapulaFc(~logicSide,:),scapulaVc,'bw','k',1);
    hp2 = gpatch(scapulaFc(logicSide,:),scapulaVc,'rw','k',1);
    legend([hp1 hp2],{'Surface above plane','Surface below plane'})
    axisGeom; axis manual; camlight headligth;
    colormap gjet;
    set(gca,'FontSize',25);
    %Plot extracted surface and boundary
    subplot(1,2,2); hold on;
    gpatch(scapulaFc(logicSide,:),scapulaVc,'w','k',1);
    gpatch(scapulaFc(~logicSide,:),scapulaVc,'w','none',0.25);
    hp1=gpatch(scapulaEc,scapulaVc,'none','b',1,3);
    hp2=quiverVec(P,n,50,'k');
    legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
    axisGeom; axis manual; camlight headligth;
    set(gca,'FontSize',25);
end

%Extract the faces we want to keep
[scapulaKeepF,scapulaKeepV] = patchCleanUnused(scapulaFc(logicSide,:),scapulaVc);

%Use the grouping function to split the extra part of the scapula that
%comes through with the cut away
[indV,indF] = groupVertices(scapulaKeepF,scapulaKeepV,1);

%Visualise the grouping
if generatePlots
    %Plot ungrouped sections
    cFigure; subplot(1,2,1); hold on;
    title('Ungrouped')
    gpatch(scapulaKeepF,scapulaKeepV,'kw','none');
    axisGeom;
    camlight headlight;
    %Plot grouped sections
    subplot(1,2,2); hold on;
    title('Grouped')
    gpatch(scapulaKeepF,scapulaKeepV,indF,'none');
    axisGeom;
    camlight headlight;
    colormap gjet; icolorbar;
end

%Extract set 2, which is the glenoid
logicKeep = logical(indF == 2);
[extractedGlenoidF,extractedGlenoidV] = patchCleanUnused(scapulaKeepF(logicKeep,:),scapulaKeepV);

%Merge vertices
[extractedGlenoidF,extractedGlenoidV] = mergeVertices(extractedGlenoidF,extractedGlenoidV);

%Self triangulate the potentially jagged edge of the cut
extractedGlenoidEb = patchBoundary(extractedGlenoidF,extractedGlenoidV); %Get boundary edges
indBoundary = edgeListToCurve(extractedGlenoidEb); %Convert boundary edges to a curve list
indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
angleThreshold = pi*(120/180); %threshold for self triangulation
[extractedGlenoidF,extractedGlenoidV,indBoundaryBack] = ...
    triSurfSelfTriangulateBoundary(extractedGlenoidF,extractedGlenoidV,indBoundary,angleThreshold,1);

%Force boundary to have a Z level aligned with the cut
extractedGlenoidV(indBoundaryBack,3) = cutLevel;

%Visualise the boundary on the cut surface
if generatePlots
    cFigure; hold on;
    gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
    plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
    camlight('headlight');
    axisGeom;
end

%Create a surface that closes the back of the glenoid
%uses 1.5 point spacing
[backF,backV] = regionTriMesh2D({extractedGlenoidV(indBoundaryBack,[1 2])},1.5,0,0);
backV(:,3) = mean(extractedGlenoidV(indBoundaryBack,3)); %Add/set z-level converting to 3D mesh

%Visualise new meshes
if generatePlots
    cFigure; hold on;
    gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
    gpatch(backF,backV,'gw','k');
    plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
    camlight('headlight');
    axisGeom;
end

%Join the two element sets
[glenoidF,glenoidV,glenoidC] = joinElementSets({extractedGlenoidF,backF},{extractedGlenoidV,backV});

%Merge vertices
[glenoidF,glenoidV] = mergeVertices(glenoidF,glenoidV);

%Visualise the joined sets
if generatePlots
    cFigure;
    gpatch(glenoidF,glenoidV,glenoidC,'k');
    colormap gjet; icolorbar;
    axisGeom;
end

%Check boundaries to make sure there are no holes. If there is throw an
%error as the volumetric meshing won't work.
if ~isempty(patchBoundary(glenoidF,glenoidV))
    error('Holes detected in glenoid. Stopping here as volumetric meshing won''t work');
end

%% Mesh the surfaces using TetGen

%Humeral head

%Create tetgen input structure
headInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
headInputStruct.Faces = headF; %Boundary faces
headInputStruct.Nodes = headV; %Nodes of boundary
headInputStruct.regionPoints= getInnerPoint(headF,headV); %Interior points for regions
headInputStruct.holePoints = []; %Interior points for holes
headInputStruct.regionA = tetVolMeanEst(headF,headV); %Desired tetrahedral volume for each region

%Mesh model
[headMesh] = runTetGen(headInputStruct); %Run tetGen

%Get outputs of mesh structure
headVolE = headMesh.elements; %The elements
headVolV = headMesh.nodes; %The vertices or nodes
headVolCE = headMesh.elementMaterialID; %Element material or region id
headVolFb = headMesh.facesBoundary; %The boundary faces
headVolCb = headMesh.boundaryMarker; %The boundary markers

%Visualise humeral head mesh
if generatePlots
    meshView(headMesh);
    title('Tetrahedral Mesh: Humeral Head','FontSize',25);
end


%Glenoid

%Create tetgen input structure
glenoidInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
glenoidInputStruct.Faces = glenoidF; %Boundary faces
glenoidInputStruct.Nodes = glenoidV; %Nodes of boundary
glenoidInputStruct.regionPoints= getInnerPoint(glenoidF,glenoidV); %Interior points for regions
glenoidInputStruct.holePoints = []; %Interior points for holes
glenoidInputStruct.regionA = tetVolMeanEst(glenoidF,glenoidV); %Desired tetrahedral volume for each region

%Mesh model
[glenoidMesh] = runTetGen(glenoidInputStruct); %Run tetGen

%Get outputs of mesh structure
glenoidVolE = glenoidMesh.elements; %The elements
glenoidVolV = glenoidMesh.nodes; %The vertices or nodes
glenoidVolCE = glenoidMesh.elementMaterialID; %Element material or region id
glenoidVolFb = glenoidMesh.facesBoundary; %The boundary faces
glenoidVolCb = glenoidMesh.boundaryMarker; %The boundary markers

%Visualise glenoid mesh
if generatePlots
    meshView(glenoidMesh);
    title('Tetrahedral Mesh: Glenoid','FontSize',25);
end


%Glenoid cartilage

%Create tetgen input structure
glenoidCartilageInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
glenoidCartilageInputStruct.Faces = glenoidCartilageF; %Boundary faces
glenoidCartilageInputStruct.Nodes = glenoidCartilageV; %Nodes of boundary
glenoidCartilageInputStruct.regionPoints= getInnerPoint(glenoidCartilageF,glenoidCartilageV); %Interior points for regions
glenoidCartilageInputStruct.holePoints = []; %Interior points for holes
glenoidCartilageInputStruct.regionA = tetVolMeanEst(glenoidCartilageF,glenoidCartilageV); %Desired tetrahedral volume for each region

%Mesh model
[glenoidCartilageMesh] = runTetGen(glenoidCartilageInputStruct); %Run tetGen

%Get outputs of mesh structure
glenoidCartilageVolE = glenoidCartilageMesh.elements; %The elements
glenoidCartilageVolV = glenoidCartilageMesh.nodes; %The vertices or nodes
glenoidCartilageVolCE = glenoidCartilageMesh.elementMaterialID; %Element material or region id
glenoidCartilageVolFb = glenoidCartilageMesh.facesBoundary; %The boundary faces
glenoidCartilageVolCb = glenoidCartilageMesh.boundaryMarker; %The boundary markers

%Visualise glenoid mesh
if generatePlots
    meshView(glenoidCartilageMesh);
    title('Tetrahedral Mesh: Glenoid Cartilage','FontSize',25);
end

%% Define contact surfaces

%Glenoid cartilage to glenoid


%Humeral head to glenoid cartilage
%%%%% TODO: adapt this given that it should be cartilage to cartilage



%% FEBio details

%%%%% EVERYTHING SOMEWHAT WORKING AT THE MOMENT...

%Defining file names
%%%%% TODO: ensure these work without full file paths, or adapt appropriately
febioFebFileNamePart = 'baseModel';
febioFebFileName = [febioFebFileNamePart,'.feb']; %FEB file name
febioLogFileName = [febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force = [febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Material parameter set
%%%%% TODO: adapt appropriately, currently don't match what is needed
c = 1e-3; %Shear-modulus-like parameter
m = 8; %Material parameter setting degree of non-linearity
k_factor = 1e2; %Bulk modulus factor
k = c*k_factor; %Bulk modulus

%FEA control settings
%%%%% TODO: adapt appropriately
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 10; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size

%Contact parameters
contactAlg = 2;
switch contactAlg
    case 1
        contactType='sticky';
    case 2
        contactType='facet-to-facet sliding';
    case 3
        contactType='sliding_with_gaps';
    case 4
        contactType='sliding2';
end

%Get a template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type = 'static';
stepStruct.Control.time_steps = numTimeSteps;
stepStruct.Control.step_size = 1/numTimeSteps; %0.1 right now
stepStruct.Control.time_stepper.dtmin = dtmin;
stepStruct.Control.time_stepper.dtmax = dtmax;
stepStruct.Control.time_stepper.max_retries = max_retries;
stepStruct.Control.time_stepper.opt_iter = opt_iter;
stepStruct.Control.max_refs = max_refs;
stepStruct.Control.max_ups = max_ups;

%Add template based default settings to proposed control section
[stepStruct.Control] = structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec = rmfield(febio_spec,'Control');

%Step specific control section
%%%%% TODO: only one step for now
febio_spec.Step{1}.Control = stepStruct.Control;
febio_spec.Step{1}.ATTR.id = 1;
% % % febio_spec.Step{2}.Control = stepStruct.Control;
% % % febio_spec.Step{2}.ATTR.id = 2;

%Material section
%Rigid body for humeral head
febio_spec.Material.material{1}.ATTR.type = 'rigid body';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.density = 1;
febio_spec.Material.material{1}.center_of_mass = mean(headV,1);
%Rigid body for glenoid
febio_spec.Material.material{2}.ATTR.type = 'rigid body';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.density = 1;
febio_spec.Material.material{2}.center_of_mass = mean(glenoidV,1);

% % % %A default like Ogden material
% % % %%%%% TODO: applying to glenoid for now - adapt
% % % febio_spec.Material.material{2}.ATTR.type = 'Ogden';
% % % febio_spec.Material.material{2}.ATTR.id = 2;
% % % febio_spec.Material.material{2}.c1 = c;
% % % febio_spec.Material.material{2}.m1 = m;
% % % febio_spec.Material.material{2}.c2 = c;
% % % febio_spec.Material.material{2}.m2 = -m;
% % % febio_spec.Material.material{2}.k = k;

%Geometry section
% -> Nodes
%%%%% TODO: adapt with additional bodies
%Combine node sets
V = [headVolV;glenoidVolV];
%Fixed element indices
glenoidVolE = glenoidVolE+size(headVolV,1);
%Set
%%%% COULD STILL TRY SPLITTING TO MULTIPLE NODES, but this works...I think...
febio_spec.Geometry.Nodes{1}.ATTR.name = 'allNodes'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = V; %The nodel coordinates

% -> Elements
%Humeral head
febio_spec.Geometry.Elements{1}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name = 'HumeralHead'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(headVolE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = headVolE;
%Glenoid
febio_spec.Geometry.Elements{2}.ATTR.type = 'tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name = 'Glenoid'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id = (1:1:size(glenoidVolE,1))' + length(headVolE); %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL = glenoidVolE;

% % % % -> NodeSets
%%%%% TODO: this might be the place to set something that makes collecting
%%%%% displacement data easier...
% % % febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
% % % febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%%%%% TODO: add surfaces when contact added
%%%%% Seems like this is needed, as the faces aren't anywhere in the data
%%%%% at the moment...
% % % % -> Surfaces
% % % febio_spec.Geometry.Surface{1}.ATTR.name='contact_master';
% % % febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(F_contact_master,1))';
% % % febio_spec.Geometry.Surface{1}.tri3.VAL=F_contact_master;
% % % 
% % % febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave';
% % % febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_slave,1))';
% % % febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_slave;

% % % % -> Surface pairs
% % % febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1';
% % % febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
% % % febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

% % % %Boundary condition section
% % % % -> Fix boundary conditions
% % % febio_spec.Boundary.fix{1}.ATTR.bc='x';
% % % febio_spec.Boundary.fix{1}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% % % febio_spec.Boundary.fix{2}.ATTR.bc='y';
% % % febio_spec.Boundary.fix{2}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;
% % % febio_spec.Boundary.fix{3}.ATTR.bc='z';
% % % febio_spec.Boundary.fix{3}.ATTR.node_set=febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body
%Humeral head
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL = 45; %prescribed displacement of 45mm
%Glenoid
febio_spec.Step{1}.Boundary.rigid_body{2}.ATTR.mat = 2;
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{3}.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{4}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{5}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{6}.ATTR.bc = 'Rz';

% % % %Contact section
% % % switch contactType
% % %     case 'sticky'
% % %         febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
% % %         febio_spec.Contact.contact{1}.ATTR.type='sticky';
% % %         febio_spec.Contact.contact{1}.penalty=100;
% % %         febio_spec.Contact.contact{1}.laugon=0;
% % %         febio_spec.Contact.contact{1}.tolerance=0.1;
% % %         febio_spec.Contact.contact{1}.minaug=0;
% % %         febio_spec.Contact.contact{1}.maxaug=10;
% % %         febio_spec.Contact.contact{1}.snap_tol=0;
% % %         febio_spec.Contact.contact{1}.max_traction=0;
% % %         febio_spec.Contact.contact{1}.search_tolerance=0.1;
% % %     case 'facet-to-facet sliding'
% % %         febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
% % %         febio_spec.Contact.contact{1}.ATTR.type='facet-to-facet sliding';
% % %         febio_spec.Contact.contact{1}.penalty=100;
% % %         febio_spec.Contact.contact{1}.auto_penalty=1;
% % %         febio_spec.Contact.contact{1}.two_pass=0;
% % %         febio_spec.Contact.contact{1}.laugon=0;
% % %         febio_spec.Contact.contact{1}.tolerance=0.1;
% % %         febio_spec.Contact.contact{1}.gaptol=0;
% % %         febio_spec.Contact.contact{1}.minaug=0;
% % %         febio_spec.Contact.contact{1}.maxaug=10;
% % %         febio_spec.Contact.contact{1}.search_tol=0.01;
% % %         febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
% % %     case 'sliding_with_gaps'
% % %         febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
% % %         febio_spec.Contact.contact{1}.ATTR.type='sliding_with_gaps';
% % %         febio_spec.Contact.contact{1}.penalty=100;
% % %         febio_spec.Contact.contact{1}.auto_penalty=1;
% % %         febio_spec.Contact.contact{1}.two_pass=0;
% % %         febio_spec.Contact.contact{1}.laugon=0;
% % %         febio_spec.Contact.contact{1}.tolerance=0.1;
% % %         febio_spec.Contact.contact{1}.gaptol=0;
% % %         febio_spec.Contact.contact{1}.minaug=0;
% % %         febio_spec.Contact.contact{1}.maxaug=10;
% % %         febio_spec.Contact.contact{1}.fric_coeff=0;
% % %         febio_spec.Contact.contact{1}.fric_penalty=0;
% % %         febio_spec.Contact.contact{1}.ktmult=1;
% % %         febio_spec.Contact.contact{1}.seg_up=0;
% % %         febio_spec.Contact.contact{1}.search_tol=0.01;
% % %     case 'sliding2'
% % %         febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
% % %         febio_spec.Contact.contact{1}.ATTR.type='sliding2';
% % %         febio_spec.Contact.contact{1}.penalty=30;
% % %         febio_spec.Contact.contact{1}.auto_penalty=1;
% % %         febio_spec.Contact.contact{1}.two_pass=0;
% % %         febio_spec.Contact.contact{1}.laugon=0;
% % %         febio_spec.Contact.contact{1}.tolerance=0.1;
% % %         febio_spec.Contact.contact{1}.gaptol=0;
% % %         febio_spec.Contact.contact{1}.symmetric_stiffness=0;
% % %         febio_spec.Contact.contact{1}.search_tol=0.01;
% % %         febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
% % % end

% LoadData section
%%%% TODO: current basic load curve for displacement testing
febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1];

% % % febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
% % % febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
% % % febio_spec.LoadData.loadcurve{2}.point.VAL=[0 0; 1 0; 2 1];

%Output section
% -> log file
%Collect humeral head displacement
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.node_data{1}.VAL = 1:size(headVolV,1); %%%only collect the first head set of nodes

%%%%% TODO: CURRENTLY MATCHING ID'S GIVEN THAT THE SETS HAVEN'T BEEN JOINED
%%%%% TOGETHER, DON'T KNOW WHAT THIS MEANS FOR DATA THOUGH?
    %%%%% SEEMS LIKE IT CAUSES AN ERROR
    %%%%% LOOKS LIKE NODESET NEEDS TO BE JOINED, BUT ELEMENT SETS DON'T FOR
    %%%%% THE ID NUMBERING (I.E. ELEMENT ID NUMBERING CAN START AT WHATEVER
    %%%%% IT WANTS TO...

% % % febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
% % % febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
% % % febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
% % % febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

%View febio spec file
febView(febio_spec);

%Export feb file
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode

%% Run the FEBio analysis

%Setup analysis details
febioAnalysis.run_filename = febioFebFileName; %The input file name
febioAnalysis.run_logname = febioLogFileName; %The name for the log file
febioAnalysis.disp_on = 1; %Display information on the command window
febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
febioAnalysis.runMode = 'internal';%'external';
febioAnalysis.t_check = 0.25; %Time for checking log file (dont set too small)
febioAnalysis.maxtpi = 1e99; %Max analysis time
febioAnalysis.maxLogCheckTime = 10; %Max log file checking time

%Run FEBio
clc
[runFlag] = runMonitorFEBio(febioAnalysis);


%% Import FEBio results

%Check for successful run
if runFlag == 1
    
    %Importing nodal displacements from a log file
    %%%%% TODO: comment this better and adjust variable names
    %%%%% TODO: might be a better way to extract results too...
    
    %%%%% TODO: anim8 is not working right, there is some disconnect
    %%%%% between the elements it's animating and the data...
    
    %Import the output logfile for displacements
    [time_mat, N_disp_mat,~] = importFEBio_logfile(febioLogFileName_disp); %Nodal displacements
    
    %Extract the time data
    time_mat = [0; time_mat(:)];
    
    %Extract the displacement data (raw)
    N_disp_mat = N_disp_mat(:,2:end,:);
    
    %Get data from each time step in the matrix
    %Get size of matrix and time steps
    sizImport = size(N_disp_mat);
    sizImport(3) = sizImport(3)+1;
    %Extract displacement at each step
    N_disp_mat_n = zeros(sizImport);
    N_disp_mat_n(:,:,2:end) = N_disp_mat;
    %Replace original matrix
    N_disp_mat = N_disp_mat_n;
    %Get end displacement    
    DN = N_disp_mat(:,:,end);
    DN_magnitude = sqrt(sum(DN(:,3).^2,2));
    
    %Define positions of objects across steps by taking node positions and
    %adding displacement to these 
    
    %%%%% CURRENTLY ONLY COLLECTING HEAD DISPLACEMENT (VOLUME)
    
    V_def = headVolV+DN;
    V_DEF = N_disp_mat+repmat(headVolV,[1 1 size(N_disp_mat,3)]);
    X_DEF = V_DEF(:,1,:);
    Y_DEF = V_DEF(:,2,:);
    Z_DEF = V_DEF(:,3,:);
% % %     [CF] = vertexToFaceMeasure(Fb1,DN_magnitude);

    %Plot simulated results using anim8
    
    %%%%% TODO: using elements to plot faces here means no faces
    %%%%% visualised...
    
    %Create basic view and store graphics handle to initiate animation
    hf = cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1 = gpatch(headVolFb,V,'gw'); %Add graphics object to animate
    %%%%% IF WANT TO USE glenoidVolFb NEED TO ADD INDEXING TO THIS LIKE
    %%%%% EARLIER WITH THE GLENOID VOLUME ELEMENT ID'S
    hp2 = gpatch(glenoidVolE,V,'bw'); %Add graphics object to animate
% % %     hp2=gpatch(E2,V_def,'kw','none',faceAlpha2); %Add graphics object to animate
% % %     gpatch(Fb1,V,0.5*ones(1,3),'none',0.25); %A static graphics object

    axisGeom(gca,25);
% % %     colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    
    %Set axes to min and max displacement values
% % %     axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    %The above isn't appropriate when the glenoid is considered, as it
    %goes outside of the humeral head displacement we recorded. We'l only
    %adjust the X-axis values for now as this is where the head goes...
    ax = gca;
    axis([min(X_DEF(:)) max(X_DEF(:)) ax.YLim(1) ax.YLim(2) ax.ZLim(1) ax.ZLim(2)]);
    camlight headlight;
    
    %Set up animation features
    animStruct.Time = time_mat; %The time vector
    for qt = 1:1:size(N_disp_mat,3) %Loop over time increments
        
        %Get the current displacement position from matrix
        DN = N_disp_mat(:,:,qt);
        
        %Get the current displacement magnitude
        DN_magnitude = sqrt(sum(DN.^2,2));
        
        %Get the current nodal coordinates, factoring in original position
        %and current displacement
        V_def = headVolV+DN;
        
% % %         [CF]=vertexToFaceMeasure(Fb1,DN_magnitude); %Current color data to use

        %Set entries in animation structure
        animStruct.Handles{qt} = [hp1]; %[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt} = {'Vertices'}; %{'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt} = {V_def};%{V_def,CF,V_def}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature

end
