
%%%%% enter notes...references to relevant papers such as Walia et al.
%%%%% studies...Klemt et al. study

%%%%% following a lot of the sphere sliding demo in GIBBON
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0007_sphere_sliding.html
%%%%% a fair bit for slicing and fixing up the hole also comes from the foot insole demo
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0050_foot_insole_01.html

%%%%% TODO: add generatePlots logical input
    
generatePlots = false;

%%%%% TODO: mesh convergence...

%%%%% TODO: identify displacement relative to 50N force to convert to
%%%%% static problem...

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

%Humeral cartilage

%Create tetgen input structure
humeralCartilageInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
humeralCartilageInputStruct.Faces = humeralCartilageF; %Boundary faces
humeralCartilageInputStruct.Nodes = humeralCartilageV; %Nodes of boundary
humeralCartilageInputStruct.regionPoints= getInnerPoint(humeralCartilageF,humeralCartilageV); %Interior points for regions
humeralCartilageInputStruct.holePoints = []; %Interior points for holes
humeralCartilageInputStruct.regionA = tetVolMeanEst(humeralCartilageF,humeralCartilageV); %Desired tetrahedral volume for each region

%Mesh model
[humeralCartilageMesh] = runTetGen(humeralCartilageInputStruct); %Run tetGen

%Get outputs of mesh structure
humeralCartilageVolE = humeralCartilageMesh.elements; %The elements
humeralCartilageVolV = humeralCartilageMesh.nodes; %The vertices or nodes
humeralCartilageVolCE = humeralCartilageMesh.elementMaterialID; %Element material or region id
humeralCartilageVolFb = humeralCartilageMesh.facesBoundary; %The boundary faces
humeralCartilageVolCb = humeralCartilageMesh.boundaryMarker; %The boundary markers

%Visualise glenoid mesh
if generatePlots
    meshView(humeralCartilageMesh);
    title('Tetrahedral Mesh: Humeral Cartilage','FontSize',25);
end

%Visualise all volumes together
if generatePlots
    cFigure; hold on
    gpatch(glenoidVolFb,glenoidVolV,'kw','none');
    gpatch(headVolFb,headVolV,'kw','none');
    gpatch(glenoidCartilageVolFb,glenoidCartilageVolV,'bw','none');
    gpatch(humeralCartilageVolFb,humeralCartilageVolV,'bw','none');
    title('Complete Volumetric Base System');
    camlight headlight;
    axisGeom;
end

%% Join the node sets together and specify new variables for indexing

%Combine the node sets into one space
allV = [glenoidVolV; glenoidCartilageVolV; headVolV; humeralCartilageVolV];
    
%New indexing for faces against nodes
all_glenoidVolFb = glenoidVolFb;
all_glenoidCartilageVolFb = glenoidCartilageVolFb + size(glenoidVolV,1);
all_headVolFb = headVolFb + (size(glenoidVolV,1) + size(glenoidCartilageVolV,1));
all_humeralCartilageVolFb = humeralCartilageVolFb + (size(glenoidVolV,1) + size(glenoidCartilageVolV,1) + size(headVolV,1));

%New indexing for elements against nodes
all_glenoidVolE = glenoidVolE;
all_glenoidCartilageVolE = glenoidCartilageVolE + size(glenoidVolV,1);
all_headVolE = headVolE + (size(glenoidVolV,1) + size(glenoidCartilageVolV,1));
all_humeralCartilageVolE = humeralCartilageVolE + (size(glenoidVolV,1) + size(glenoidCartilageVolV,1) + size(headVolV,1));

%New indexing for node ID
all_glenoidVolV_ind = (1:1:size(glenoidVolV,1))';
all_glenoidCartilageVolV_ind = (length(glenoidVolV)+1:1:size(glenoidCartilageVolV,1)+length(glenoidVolV))';
all_headVolV_ind = (length(glenoidVolV)+length(glenoidCartilageVolV)+1:1:size(headVolV,1)+length(glenoidVolV)+length(glenoidCartilageVolV))';
all_humeralCartilageVolV_ind = (length(glenoidVolV)+length(glenoidCartilageVolV)+length(headVolV)+1:1:size(humeralCartilageVolV,1)+length(glenoidVolV)+length(glenoidCartilageVolV)+length(headVolV))';

%New indexing for element ID
all_glenoidVolE_ind = (1:1:size(glenoidVolE,1))';
all_glenoidCartilageVolE_ind = (length(glenoidVolE)+1:1:size(glenoidCartilageVolE,1)+length(glenoidVolE))';
all_headVolE_ind = (length(glenoidVolE)+length(glenoidCartilageVolE)+1:1:size(headVolE,1)+length(glenoidVolE)+length(glenoidCartilageVolE))';
all_humeralCartilageVolE_ind = (length(glenoidVolE)+length(glenoidCartilageVolE)+length(headVolE)+1:1:size(humeralCartilageVolE,1)+length(glenoidVolE)+length(glenoidCartilageVolE)+length(headVolE))';

%New indexing for face ID
all_glenoidVolFb_ind = (1:1:size(glenoidVolFb,1))';
all_glenoidCartilageVolFb_ind = (length(glenoidVolFb)+1:1:size(glenoidCartilageVolFb,1)+length(glenoidVolFb))';
all_headVolFb_ind = (length(glenoidVolFb)+length(glenoidCartilageVolFb)+1:1:size(headVolFb,1)+length(glenoidVolFb)+length(glenoidCartilageVolFb))';
all_humeralCartilageVolFb_ind = (length(glenoidVolFb)+length(glenoidCartilageVolFb)+length(headVolFb)+1:1:size(humeralCartilageVolFb,1)+length(glenoidVolFb)+length(glenoidCartilageVolFb)+length(headVolFb))';

%Visualise joined volumes
if generatePlots
    cFigure; hold on
    gpatch(all_glenoidVolFb,allV,'kw','none');
    gpatch(all_headVolFb,allV,'kw','none');
    gpatch(all_glenoidCartilageVolFb,allV,'bw','none');
    gpatch(all_humeralCartilageVolFb,allV,'bw','none');
    title('Joined Base System');
    camlight headlight;
    axisGeom;
end

%%

%% Define contact surfaces

%Glenoid cartilage to glenoid


%Humeral head to glenoid cartilage
%%%%% TODO: adapt this given that it should be cartilage to cartilage



%% Build FEBio model file up

%%%%% unsure if below is working well --- see febio compiled for example,
%%%%% parts seemingly should work...? id's for everything do need to
%%%%% incrementally increase...tet4 instead of tri3 for elements --- tri3
%%%%% for surfaces...

%% Create structure and defaults

%%%%% TODO: Febio studio uses the febio3 executable, which isn't what the
%%%%% GIBBON toolbox is using -- so this could cause some of the issues
%%%%% between the FEBio studio FE model definitions and GIBBON coding...?

%Create FEBio template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%% Model controls

%FEA control settings
%%%%% TODO: adapt appropriately
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 10; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size

%Add initial compression step
febio_spec.Step{1}.ATTR.name = 'initialCompression';

%Create control structure for first step
stepStruct.Control.time_steps = numTimeSteps;
stepStruct.Control.step_size = 1/numTimeSteps; %0.1 right now
stepStruct.Control.max_refs = max_refs;
stepStruct.Control.max_ups = max_ups;
stepStruct.Control.max_ups = max_ups;
stepStruct.Control.time_stepper.dtmin = dtmin;
stepStruct.Control.time_stepper.dtmax = dtmax;
stepStruct.Control.time_stepper.max_retries = max_retries;
stepStruct.Control.time_stepper.opt_iter = opt_iter;
stepStruct.Control.analysis.ATTR.type = 'dynamic'; %Dynamic appears necessary when running body loads

%Add template based default settings to proposed control section
[stepStruct.Control] = structComplete(stepStruct.Control,febio_spec.Control,1); %Complement provided with default if missing

%Remove control field (part of template) since step specific control sections are used
febio_spec = rmfield(febio_spec,'Control');

%%%%% TODO: remove the elements of the stepStruct defaults that FEBio
%%%%% studio seems to reject -- see errors on import...
    %%%%% This includes the tags plot_range, plot_zero_state, cmax,
    %%%%% print_level, output_level

%Add control structure to initial compression step
febio_spec.Step{1}.Control = stepStruct.Control;

% % % febio_spec.Step{2}.Control = stepStruct.Control;
% % % febio_spec.Step{2}.ATTR.id = 2;

%% Model materials

%Set constants for incompressible neo-Hookean materials based on Youngs
%modulus and Poisson's ratio from Walia et al. (2015)
%Set Youngs and Poisson
E = 10;
v = 0.4;
%Calculate constants using Walia et al. (2015) equations
%Multiplication by factor of 1000 appears necessary to get appropriate
%behaviour and likely has to do with variable mesh units
%%%%% TODO: check relation of these to the shear and bulk factors in FEBio
c10 = (E / (4*(1+v)))*1000;
d10 = (E / (6*(1-2*v)))*1000;

%Glenoid

%Rigid body parameters
febio_spec.Material.material{1}.ATTR.type = 'rigid body';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.ATTR.name = 'rigidGlenoid';
febio_spec.Material.material{1}.density = 1;
febio_spec.Material.material{1}.center_of_mass = mean(glenoidV,1);

%Glenoid cartilage

%Elastic material parameters
febio_spec.Material.material{2}.ATTR.type = 'incomp neo-Hookean';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.ATTR.name = 'elasticGlenoidCartilage';
febio_spec.Material.material{2}.density = 1;
febio_spec.Material.material{2}.G = c10;
febio_spec.Material.material{2}.k = d10;

%Humeral Head

%Rigid body parameters
febio_spec.Material.material{3}.ATTR.type = 'rigid body';
febio_spec.Material.material{3}.ATTR.id = 3;
febio_spec.Material.material{3}.ATTR.name = 'rigidHead';
febio_spec.Material.material{3}.density = 1;
febio_spec.Material.material{3}.center_of_mass = mean(headV,1);

%Humeral cartilage

%Elastic material parameters
febio_spec.Material.material{4}.ATTR.type = 'incomp neo-Hookean';
febio_spec.Material.material{4}.ATTR.id = 4;
febio_spec.Material.material{4}.ATTR.name = 'elasticHumeralCartilage';
febio_spec.Material.material{4}.density = 1;
febio_spec.Material.material{4}.G = c10;
febio_spec.Material.material{4}.k = d10;

%% Model geometry

%Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name = 'allNodes'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(allV,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = allV; %The nodel coordinates

%Elements

%Glenoid
febio_spec.Geometry.Elements{1}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material name for this set; can also use index if this doesn't work
febio_spec.Geometry.Elements{1}.ATTR.name = 'glenoidElements'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = all_glenoidVolE_ind; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = all_glenoidVolE;

%Glenoid cartilage
febio_spec.Geometry.Elements{2}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material name for this set; can also use index if this doesn't work
febio_spec.Geometry.Elements{2}.ATTR.name = 'glenoidCartilageElements'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id = all_glenoidCartilageVolE_ind; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL = all_glenoidCartilageVolE;

%Humeral head
febio_spec.Geometry.Elements{3}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat = 3; %material name for this set; can also use index if this doesn't work
febio_spec.Geometry.Elements{3}.ATTR.name = 'headElements'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id = all_headVolE_ind; %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL = all_headVolE;

%Humeral cartilage
febio_spec.Geometry.Elements{4}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{4}.ATTR.mat = 4; %material name for this set; can also use index if this doesn't work
febio_spec.Geometry.Elements{4}.ATTR.name = 'humeralCartilageElements'; %Name of the element set
febio_spec.Geometry.Elements{4}.elem.ATTR.id = all_humeralCartilageVolE_ind; %Element id's
febio_spec.Geometry.Elements{4}.elem.VAL = all_humeralCartilageVolE;

%Node sets

% % % %Glenoid nodes
% % % febio_spec.Geometry.NodeSet{1}.ATTR.name = 'glenoidNodes';
% % % febio_spec.Geometry.NodeSet{1}.node.ATTR.id = all_glenoidVolV_ind;
% % % 
% % % %Glenoid cartilage nodes
% % % febio_spec.Geometry.NodeSet{2}.ATTR.name = 'glenoidCartilageNodes';
% % % febio_spec.Geometry.NodeSet{2}.node.ATTR.id = all_glenoidCartilageVolV_ind;
% % % 
% % % %Humeral head nodes
% % % febio_spec.Geometry.NodeSet{3}.ATTR.name = 'headNodes';
% % % febio_spec.Geometry.NodeSet{3}.node.ATTR.id = all_headVolV_ind;
% % % 
% % % %Humeral cartilage nodes
% % % febio_spec.Geometry.NodeSet{4}.ATTR.name = 'humeralCartilageNodes';
% % % febio_spec.Geometry.NodeSet{4}.node.ATTR.id = all_humeralCartilageVolV_ind;

%Surfaces

%%%%%% TODO: surfaces don't seem to come across to FEBio studio, other
%%%%%% surfaces are automatically identified...

%Glenoid
febio_spec.Geometry.Surface{1}.ATTR.name = 'glenoidSurface';
% % % febio_spec.Geometry.Surface{1}.tri3.ATTR.id = all_glenoidVolFb_ind;
febio_spec.Geometry.Surface{1}.tri3.ATTR.id = (1:1:size(all_glenoidVolFb_ind,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL = all_glenoidVolFb;

%Glenoid cartilage
febio_spec.Geometry.Surface{2}.ATTR.name = 'glenoidCartilageSurface';
% % % febio_spec.Geometry.Surface{2}.tri3.ATTR.id = all_glenoidCartilageVolFb_ind;
febio_spec.Geometry.Surface{2}.tri3.ATTR.id = (1:1:size(all_glenoidCartilageVolFb_ind,1))';
febio_spec.Geometry.Surface{2}.tri3.VAL = all_glenoidCartilageVolFb;

%Humeral head
febio_spec.Geometry.Surface{3}.ATTR.name = 'headSurface';
% % % febio_spec.Geometry.Surface{3}.tri3.ATTR.id = all_headVolFb_ind;
febio_spec.Geometry.Surface{3}.tri3.ATTR.id = (1:1:size(all_headVolFb_ind,1))';
febio_spec.Geometry.Surface{3}.tri3.VAL = all_headVolFb;

%Glenoid cartilage
febio_spec.Geometry.Surface{4}.ATTR.name = 'humeralCartilageSurface';
% % % febio_spec.Geometry.Surface{4}.tri3.ATTR.id = all_humeralCartilageVolFb_ind;
febio_spec.Geometry.Surface{4}.tri3.ATTR.id = (1:1:size(all_humeralCartilageVolFb_ind,1))';
febio_spec.Geometry.Surface{4}.tri3.VAL = all_humeralCartilageVolFb;

%% Model contact interfaces

%The rigid body to elastic cartilage contact is described by a tied facet to
%facet interface as described by Klemt et al. (2019). The contact between
%the cartilage surfaces is described by a sliding facet to facet interface
%as per Walia et al. (2013).

%Humeral head to cartilage

% % % %Create surface pair
% % % febio_spec.Geometry.SurfacePair{1}.ATTR.name = 'surfacePair_headToCartilage';
% % % ebio_spec.Geometry.SurfacePair{1}.master.ATTR.surface = 'humeralCartilageSurface';
% % % febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface = 'headSurface';
% % % 
% % % %Contact parameters
% % % febio_spec.Contact.contact{1}.ATTR.name = 'headToCartilageContact';
% % % febio_spec.Contact.contact{1}.ATTR.type = 'tied-facet-on-facet';
% % % febio_spec.Contact.contact{1}.ATTR.surface_pair = 'surfacePair_headToCartilage';
% % % febio_spec.Contact.contact{1}.laugon = 0;
% % % febio_spec.Contact.contact{1}.tolerance = 0.1;
% % % febio_spec.Contact.contact{1}.penalty = 1000;
% % % febio_spec.Contact.contact{1}.minaug = 0;
% % % febio_spec.Contact.contact{1}.maxaug = 10;
% % % 
% % % %Glenoid to cartilage
% % % 
% % % %Create surface pair
% % % febio_spec.Geometry.SurfacePair{2}.ATTR.name = 'surfacePair_glenoidToCartilage';
% % % ebio_spec.Geometry.SurfacePair{2}.master.ATTR.surface = 'glenoidCartilageSurface';
% % % febio_spec.Geometry.SurfacePair{2}.slave.ATTR.surface = 'glenoidSurface';
% % % 
% % % %Contact parameters
% % % febio_spec.Contact.contact{2}.ATTR.name = 'glenoidToCartilageContact';
% % % febio_spec.Contact.contact{2}.ATTR.type = 'tied-facet-on-facet';
% % % febio_spec.Contact.contact{2}.ATTR.surface_pair = 'surfacePair_glenoidToCartilage';
% % % febio_spec.Contact.contact{2}.laugon = 0;
% % % febio_spec.Contact.contact{2}.tolerance = 0.1;
% % % febio_spec.Contact.contact{2}.penalty = 1000;
% % % febio_spec.Contact.contact{2}.minaug = 0;
% % % febio_spec.Contact.contact{2}.maxaug = 10;

%% Model loads

%%%%% TODO: ensure that this basic 50 application relates to a 50N force
%%%%% being applied to the body - a basic sphere with density and mass
%%%%% against displacement would be good test

%Add compressive load to humeral head
febio_spec.Loads.body_load{1}.ATTR.name = 'compressiveForce';
febio_spec.Loads.body_load{1}.ATTR.type = 'const';

%%%%%%% TODO: elem_set doesn't seem to get recognised in body load when
%%%%%%% importing to FEBio studio

%%%%% Maybe need to consider adding all geometry aspects under a Part
%%%%% within the geometry section?

febio_spec.Loads.body_load{1}.ATTR.elem_set = 'headElements';
febio_spec.Loads.body_load{1}.z.VAL = 50; %50N compressive force
febio_spec.Loads.body_load{1}.z.ATTR.lc = 1;

%Add load curve for compressive force. Builds up to maximum over first second
febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1];

%% Model boundary conditions

%Boundary conditions on rigid bodies from step 1

%Glenoid is completely fixed
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{6}.ATTR.bc = 'Rz';

%Humeral head is fixed with prescribed displacement
%%%%% TODO: this needs to be edited after adding first force control step
%Fixed conditions
febio_spec.Step{1}.Boundary.rigid_body{2}.ATTR.mat = 3;
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{2}.fixed{5}.ATTR.bc = 'Rz';
% % % %Prescribed conditions
% % % febio_spec.Step{1}.Boundary.rigid_body{3}.ATTR.mat = 3;
% % % febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.ATTR.bc = 'x';
% % % febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.ATTR.lc = 1;
% % % febio_spec.Step{1}.Boundary.rigid_body{3}.prescribed.VAL = 45; %prescribed displacement of 45mm

%%

febioStruct2xml(febio_spec,'test.feb'); %Exporting to file and domNode

%%%%%% NOTES: simulations work without any contact in current form -- head
%%%%%% still translates along appropriately

%%%%% Adding sticky contact causes errors in both static and dynamic mode,
%%%%% negative jacobian and problem diverges...

%%%%% Displacement works well with fairly basic (see sphere sliding)
%%%%% sliding facet on facet contact parameters between the humeral head
%%%%% and humeral cartilage -- in short, with just displacement of the
%%%%% rigid body connected with facet to facet contact, things go well...

%%%%% could get force control sim to terminate normally by changing the
%%%%% step to dynamic, but nothing happens...?

%%%%% displacement control seems to work fine - contact interaction between
%%%%% cartilage surfaces starts causing issues though -- solver starts
%%%%% struggling and reduces time steps to manage this, don't think it will
%%%%% converge though -- reduced step size to 0.01 and it still can't seem
%%%%% to get through the initial time step...

%%%%% Adding force as a 'Load' and changing the step to dynamic mode seemed
%%%%% to progress a little better (without cartilage to cartilage contact).
%%%%% This is similar to the sphere cone slide body force example in
%%%%% GIBBON. Seemingly a rigid body force would work, but this might be
%%%%% the better option. Force seems to be opposite to what you'd expect
%%%%% though (+ve Z vs. -ve Z) -- this approach works fine to move the
%%%%% head, and the contact interface works well, in fact, too well -- the
%%%%% cartilage gets stretched and doesn't move rather than moving with the
%%%%% head! Body load is applied to 'headElements'

%%%%% Body load without any contact works fine as you would expect...

%%%%% Tied facet to facet contact seems to solve the problem, dragging the
%%%%% cartilage with the head. Penalty of 100 resulted in surface
%%%%% penetration though before dragging it along...1e+07 works better, but
%%%%% still not perfect...continually increasing this might end up working
%%%%% out...tried a tied-elastic contact though as well, with some basic
%%%%% parameters (e.g. auto penalty, penalty factor = 100) the solver finds
%%%%% it difficult to finish (progresses but keeps reducing time step),
%%%%% tied elastic interface seems to stop the head from moving rather than
%%%%% moving with it -- 1e+10 doesn't seem to have penetration but does
%%%%% 'squish' out the cartilage, this could be to do with material
%%%%% properties though. Simply increasing Poisson to 0.4 seemed to help
%%%%% the time stepper follow the 0.1 specification, although this was in
%%%%% conjunction with adding cartilage to cartilage facet to facet
%%%%% interface -- the time stepper did start struggling later in the
%%%%% problem though, likely when cartilage contact was occurring...

%%%%% cartilage to cartilage contact failed when using facet to facet auto
%%%%% penalty and 100 penalty - turned off auto penalty and cranked up
%%%%% penalty to 1e+10 --- solver still fails at cartilage contact

%%%%% Refining humeral head mesh could improve contact interfaces, it is
%%%%% somewhat coarse at the moment

%%%%% No matter what's contacting it, facet on facet contact seems
%%%%% problematic when the glenoid cartilage is considered

%%%%% Note that rigid surface is recommended as master for contacts, I've
%%%%% been doing the opposite for this in all cases so far -- could improve
%%%%% the tied head to humeral cartilage contact, doesn't seem to improve
%%%%% contact with glenoid cartilage

%%%%% The laugon option is also recommended for compression dominated
%%%%% problems, given that the penalty factor needs to be so high to avoid
%%%%% penetration in these instances -- perhaps an option to combine with
%%%%% auto penalty? -- laugon and auto penalty on with penalty of 1000
%%%%% still seems to struggle with glenoid cartilage contact (in sliding
%%%%% elastic format)

%%

%%%%% EVERYTHING SOMEWHAT WORKING AT THE MOMENT...below here...above is
%%%%% still testing with doing geometry as parts...

%%

%%%% theoretically the above should work better, it seems to go into FEBio
%%%% studio but with some inconsistencies -- the below structure may work
%%%% with the older FEBio executable?

%%%%% Consider use of parts in geometry??? Node ID still needs to
%%%%% progressively increase, but may separate geometries better???

%Create FEBio template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';


%Defining file names
%%%%% TODO: ensure these work without full file paths, or adapt appropriately
febioFebFileNamePart = 'baseModel';
febioFebFileName = [febioFebFileNamePart,'.feb']; %FEB file name
febioLogFileName = [febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp = [febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
febioLogFileName_force = [febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%FEA control settings
%%%%% TODO: adapt appropriately
numTimeSteps = 10; %Number of time steps desired
max_refs = 25; %Max reforms
max_ups = 0; %Set to zero to use full-Newton iterations
opt_iter = 10; %Optimum number of iterations
max_retries = 5; %Maximum number of retires
dtmin = (1/numTimeSteps)/100; %Minimum time step size
dtmax = 1/numTimeSteps; %Maximum time step size

% % % %Contact parameters
% % % contactAlg = 2;
% % % switch contactAlg
% % %     case 1
% % %         contactType='sticky';
% % %     case 2
% % %         contactType='facet-to-facet sliding';
% % %     case 3
% % %         contactType='sliding_with_gaps';
% % %     case 4
% % %         contactType='sliding2';
% % % end



%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type = 'dynamic';
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
febio_spec.Material.material{1}.ATTR.name = 'rigidHead';
febio_spec.Material.material{1}.ATTR.id = 1;
febio_spec.Material.material{1}.density = 1;
febio_spec.Material.material{1}.center_of_mass = mean(headV,1);
%Rigid body for glenoid
febio_spec.Material.material{2}.ATTR.type = 'rigid body';
febio_spec.Material.material{2}.ATTR.name = 'rigidGlenoid';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.density = 1;
febio_spec.Material.material{2}.center_of_mass = mean(glenoidV,1);


%Set constants for incompressible neo-Hookean materials based on Youngs
%modulus and Poisson's ratio from Walia et al. (2015)
%Set Youngs and Poisson
E = 10;
v = 0.4;
%Calculate constants using Walia et al. (2015) equations
%Multiplication by factor of 1000 appears necessary to get appropriate
%behaviour and likely has to do with variable mesh units
%%%%% TODO: check relation of these to the shear and bulk factors in FEBio
c10 = (E / (4*(1+v)))*1000;
d10 = (E / (6*(1-2*v)))*1000;

%Glenoid cartilage

%Elastic material parameters
febio_spec.Material.material{3}.ATTR.type = 'incomp neo-Hookean';
febio_spec.Material.material{3}.ATTR.id = 3;
febio_spec.Material.material{3}.ATTR.name = 'elasticGlenoidCartilage';
febio_spec.Material.material{3}.density = 1;
febio_spec.Material.material{3}.G = c10;
febio_spec.Material.material{3}.k = d10;


%Humeral cartilage

%Elastic material parameters
febio_spec.Material.material{4}.ATTR.type = 'incomp neo-Hookean';
febio_spec.Material.material{4}.ATTR.id = 4;
febio_spec.Material.material{4}.ATTR.name = 'elasticHumeralCartilage';
febio_spec.Material.material{4}.density = 1;
febio_spec.Material.material{4}.G = c10;
febio_spec.Material.material{4}.k = d10;


%Geometry section
% -> Nodes
%%%%% TODO: adapt with additional bodies
%Combine node sets
V = [headVolV;glenoidVolV;glenoidCartilageVolV];
% % % V = [headVolV;glenoidVolV;glenoidCartilageVolV;humeralCartilageVolV];
%Fixed element indices
glenoidVolE = glenoidVolE+size(headVolV,1);
glenoidCartilageVolE = glenoidCartilageVolE+(size(headVolV,1)+size(glenoidVolV,1));
% % % humeralCartilageVolE = humeralCartilageVolE+(size(headVolV,1)+size(glenoidVolV,1)+size(glenoidCartilageVolV,1));
%Set
%%%% COULD STILL TRY SPLITTING TO MULTIPLE NODES, but this works...I think...
febio_spec.Geometry.Nodes{1}.ATTR.name = 'allNodes'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = V; %The nodel coordinates

% -> Elements
%Humeral head
febio_spec.Geometry.Elements{1}.ATTR.type='tet4'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name = 'HumeralHead'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(headVolE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = headVolE;
%Glenoid
febio_spec.Geometry.Elements{2}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat = 2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name = 'Glenoid'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id = (1:1:size(glenoidVolE,1))' + length(headVolE); %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL = glenoidVolE;
%Glenoid cartilage
febio_spec.Geometry.Elements{3}.ATTR.type = 'tet4'; %Element type of this set
febio_spec.Geometry.Elements{3}.ATTR.mat = 3; %material index for this set
febio_spec.Geometry.Elements{3}.ATTR.name = 'GlenoidCartilage'; %Name of the element set
febio_spec.Geometry.Elements{3}.elem.ATTR.id = (1:1:size(glenoidCartilageVolE,1))' + (length(headVolE) + length(glenoidVolE)); %Element id's
febio_spec.Geometry.Elements{3}.elem.VAL = glenoidCartilageVolE;
% % % %Humeral cartilage
% % % febio_spec.Geometry.Elements{4}.ATTR.type = 'tet4'; %Element type of this set
% % % febio_spec.Geometry.Elements{4}.ATTR.mat = 4; %material index for this set
% % % febio_spec.Geometry.Elements{4}.ATTR.name = 'HumeralCartilage'; %Name of the element set
% % % febio_spec.Geometry.Elements{4}.elem.ATTR.id = (1:1:size(humeralCartilageVolE,1))' + (length(headVolE) + length(glenoidVolE) + length(glenoidCartilageVolE)); %Element id's
% % % febio_spec.Geometry.Elements{4}.elem.VAL = humeralCartilageVolE;

% -> NodeSets
febio_spec.Geometry.NodeSet{1}.ATTR.name = 'headNodes';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id = (1:1:size(headVolV,1))';

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

%Humeral head
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
% % % febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'x';
% % % febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 1;
% % % febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL = 45; %prescribed displacement of 45mm
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

%Body loads
% % % %%%%% TODO: this doesn't seem to specify what to apply the load to?
% % % febio_spec.Loads.body_load{1}.ATTR.type = 'const';
% % % febio_spec.Loads.body_load{1}.ATTR.elem_set = 'HumeralHead';
% % % febio_spec.Loads.body_load{1}.z.VAL = 50;
% % % febio_spec.Loads.body_load{1}.z.ATTR.lc = 1;

%Nodal force
febio_spec.Loads.nodal_load{1}.ATTR.bc = 'z';
febio_spec.Loads.nodal_load{1}.ATTR.node_set = 'headNodes';
febio_spec.Loads.nodal_load{1}.scale.ATTR.lc = 1;
febio_spec.Loads.nodal_load{1}.scale.VAL = 1.0;
febio_spec.Loads.nodal_load{1}.value = 50;

%%%%% Nodal loads maybe work, but move the sphere far too much. Is it
%%%%% summing and applying 50N on each node anbd needs to be divided?
%%%%% Using a really low number errored out in FEBio though? Anything under
%%%%% 1 for either value or scale seems to cause an issue

%%%%% Even applying to one node generates massive displacement - something
%%%%% to do with the mass of the head? Mass is super low, so it makes sense
%%%%% if you apply a 50N force to nodes that you would get crazy
%%%%% displacement

%%%%% Rigid bodies don't seem to need mass for prescribed displacement, but
%%%%% probably do for nodal forces. Also rigid bodies potentially don't
%%%%% need to be volume meshed.

%%%%% Changed head to neo-Hookean material and nodal loads seemingly work
%%%%% as you'd expect. Applying desired force to a signle node on the
%%%%% outside of the humeral head would probably achieve the desired
%%%%% outcome.


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
%%%%% this may be wrong now given different nodal organisation?
febio_spec.Output.logfile.node_data{1}.VAL = (1:1:size(headVolV,1)); %%%only collect the first head set of nodes

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
% % % febView(febio_spec);

%Export feb file
febioStruct2xml(febio_spec,febioFebFileName); %Exporting to file and domNode
% % % febioStruct2xml(febio_spec,'test.feb'); %Exporting to file and domNode

%%%% Works with the glenoid cartilage elements in there, but doesn't
%%%% progress, seems like that might be something to do with the force
%%%% application

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
