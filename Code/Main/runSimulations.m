
%%%%% enter notes...

%%%%% following a lot of the sphere sliding demo in GIBBON
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0007_sphere_sliding.html

%% Current test code
%  Building up piece by piece to get to final simulation

warning off

%Navigate to surface mesh directory
cd('..\..\FEA\SurfaceMeshes');


%Load surface files
[headSTLstruct] = import_STL('BaseHumeralHead.stl');
    [scapulaSTLstruct] = import_STL('Scapula_r.stl');
[glenoidSTLstruct] = import_STL('BaseGlenoid.stl');

%Access the data from the STL structs

%Humeral head

%Get surface faces and vertices
headF = headSTLstruct.solidFaces{1}; %Faces
headV = headSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[headF,headV] = mergeVertices(headF,headV);

%Remesh sphere to reduce number of triangles
[headF,headV] = triRemeshLabel(headF,headV,3);

%% Scapula testing
%  TODO: fix variable naming throughout

%Get single surface faces and vertices
scapulaF = scapulaSTLstruct.solidFaces{1}; %Faces
scapulaV = scapulaSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);

% % % cFigure;
% % % gpatch(scapulaF,scapulaV);
% % % axisGeom;

%TEST FULL SCAPULA VOLUMETRIC MESHING

%Create tetgen input structure
inputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
inputStruct.Faces = scapulaF; %Boundary faces
inputStruct.Nodes = scapulaV; %Nodes of boundary
inputStruct.regionPoints= getInnerPoint(scapulaF,scapulaV); %Interior points for regions
inputStruct.holePoints = []; %Interior points for holes
inputStruct.regionA = tetVolMeanEst(scapulaF,scapulaV); %Desired tetrahedral volume for each region

%Mesh model
[scapulaMesh] = runTetGen(inputStruct); %Run tetGen

%Visualise scapula mesh
meshView(scapulaMesh);

%%%%% SCAPULA MESHES WITH VOLUME WITH NO ERRORS

%% Test extracting glenoid from scapula with GIBBON

%Slice the scapula 10mm back from the glenoid origin
%Generate settings for slicing
snapTolerance = mean(patchEdgeLengths(scapulaF,scapulaV))/100;
n = vecnormalize([0 0 -1]); %Normal direction to plane
P = [0 0 -10]; %Point on plane

%Slicing surface (note 3rd color data output is supressed)
[Fc,Vc,~,logicSide,Eb] = triSurfSlice(scapulaF,scapulaV,[],P,n,snapTolerance);

%Visualise slice
cFigure; subplot(1,2,1); hold on;
hp1 = gpatch(Fc(~logicSide,:),Vc,'bw','k',1);
hp2 = gpatch(Fc(logicSide,:),Vc,'rw','k',1);
legend([hp1 hp2],{'Surface above plane','Surface below plane'})
axisGeom; axis manual; camlight headligth;
colormap gjet;
set(gca,'FontSize',25);
subplot(1,2,2); hold on;
gpatch(Fc(logicSide,:),Vc,'w','k',1);
gpatch(Fc(~logicSide,:),Vc,'w','none',0.25);
hp1=gpatch(Eb,Vc,'none','b',1,3);
hp2=quiverVec(P,n,50,'k');
legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
axisGeom; axis manual; camlight headligth;
set(gca,'FontSize',25);

%Extract the faces we want to keep
[glenoidF,glenoidV] = patchCleanUnused(Fc(logicSide,:),Vc);

%Use the grouping function to split the extra part of the scapula that
%comes through with the cut away
[groupIndexVertices,groupIndexFaces] = groupVertices(glenoidF,glenoidV,1);

%Visualise the grouping
cFigure; subplot(1,2,1); hold on;
title('Ungrouped')
gpatch(glenoidF,glenoidV,'kw','none');
axisGeom;
camlight headlight;
subplot(1,2,2); hold on;
title('Grouped')
gpatch(glenoidF,glenoidV,'kw','none');
scatterV(glenoidV,15,groupIndexVertices,'filled');
axisGeom;
camlight headlight;
colormap gjet; icolorbar;

%Extract set 2, which is the glenoid
logicKeep = logical(groupIndexFaces == 2);
[extractedGlenoidF,extractedGlenoidV] = patchCleanUnused(glenoidF(logicKeep,:),glenoidV);

%Merge vertices
[extractedGlenoidF,extractedGlenoidV] = mergeVertices(extractedGlenoidF,extractedGlenoidV);

%Attempt to self triangulate potentially jagged edge
Eb = patchBoundary(extractedGlenoidF,extractedGlenoidV); %Get boundary edges
indBoundary = edgeListToCurve(Eb); %Convert boundary edges to a curve list
indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
angleThreshold = pi*(120/180); %threshold for self triangulation
[extractedGlenoidF,extractedGlenoidV,indBoundaryBack] = ...
    triSurfSelfTriangulateBoundary(extractedGlenoidF,extractedGlenoidV,indBoundary,angleThreshold,1);

%Force boundary to have the Z level chosen
%i.e. where the cut was made at -10
extractedGlenoidV(indBoundaryBack,3) = -10;

%Visualise
cFigure; hold on;
gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
camlight('headlight');
axisGeom;

%Close the back of the glenoid
[backF,backV] = regionTriMesh2D({extractedGlenoidV(indBoundaryBack,[1 2])},1.5,0,0);
backV(:,3) = mean(extractedGlenoidV(indBoundaryBack,3)); %Add/set z-level

%Visualise
cFigure; hold on;
gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
gpatch(backF,backV,'gw','k');
plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
camlight('headlight');
axisGeom;

clear glenoidF glenoidV

%Joint element sets
[glenoidF,glenoidV,glenoidC] = ...
    joinElementSets({extractedGlenoidF,backF},...
    {extractedGlenoidV,backV});

%Merge vertices
[glenoidF,glenoidV] = mergeVertices(glenoidF,glenoidV);

%Visualise
cFigure;
gpatch(glenoidF,glenoidV,glenoidC,'k');
colormap gjet; icolorbar;
axisGeom;

%Check boundaries
if isempty(patchBoundary(glenoidF,glenoidV))
    disp('No boundaries identified');
else
    disp('Boundaries identified. Oh no...!');
end

%Test volumetric meshing of glenoid

%Create tetgen input structure
inputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
inputStruct.Faces = glenoidF; %Boundary faces
inputStruct.Nodes = glenoidV; %Nodes of boundary
inputStruct.regionPoints= getInnerPoint(glenoidF,glenoidV); %Interior points for regions
inputStruct.holePoints = []; %Interior points for holes
inputStruct.regionA = tetVolMeanEst(glenoidF,glenoidV); %Desired tetrahedral volume for each region

%Mesh model
[glenoidMesh] = runTetGen(inputStruct); %Run tetGen

%Visualise scapula mesh
meshView(glenoidMesh);

%%%%% THIS APPROACH WORKS TO MESH THE GLENOID! IT IS POSSIBLE THAT THIS
%%%%% APPROACH MIGHT BE ABLE TO BE APPLIED TO THE GLENOID? NOPE...THERE ARE
%%%%% STILL SELF-INTERSECTING FACES -- IT SEEMS LIKE ANY SLICING APPLIED TO
%%%%% SURFACE MESHES MIGHT BE BEST DONE WITH GIBBON. STRAIGHT LINE CUTS
%%%%% ALONG AXES WILL BE EASY ENOUGH TO APPLY, ANGLED MIGHT BE TRICKY --
%%%%% SOLUTION MIGHT BE TO ALIGN THE ANGLE CUT TO AN AXES, CUT AND FIX, AND
%%%%% THEN REALIGN BACK TO ORIGINAL AXES...


%% This section isn't needed but may have some useful sections...


%Glenoid

%TEST different approach of filling the back with our own mesh

%%%%% TODO: fix variable labelling through here

%Start by extracting the two separate front and back surfaces
F1 = glenoidSTLstruct.solidFaces{1}; %Faces
V1 = glenoidSTLstruct.solidVertices{1}; %Vertices
% % % F2 = glenoidSTLstruct.solidFaces{2}; %Faces
% % % V2 = glenoidSTLstruct.solidVertices{2}; %Vertices

%Get a clean slice on the back of the glenoid. The slice in 3matic was
%created 10mm back from the origin point (-Z axis), so we can just sneak in
%a little bit from this
%Generate settings for slicing
snapTolerance = mean(patchEdgeLengths(F1,V1))/100;
n = vecnormalize([0 0 1]); %Normal direction to plane
P = [0 0 -9]; %Point on plane

%Slicing surface (note 3rd color data output is supressed)
[Fc,Vc,~,logicSide,~] = triSurfSlice(F1,V1,[],P,n,snapTolerance);

%Visualise slice
% % % cFigure; subplot(1,2,1); hold on;
% % % hp1 = gpatch(Fc(~logicSide,:),Vc,'bw','k',1);
% % % hp2 = gpatch(Fc(logicSide,:),Vc,'rw','k',1);
% % % legend([hp1 hp2],{'Surface above plane','Surface below plane'})
% % % axisGeom; axis manual; camlight headligth;
% % % colormap gjet;
% % % set(gca,'FontSize',25);
% % % subplot(1,2,2); hold on;
% % % gpatch(Fc(logicSide,:),Vc,'w','k',1);
% % % gpatch(Fc(~logicSide,:),Vc,'w','none',0.25);
% % % hp1=gpatch(Eb,Vc,'none','b',1,3);
% % % hp2=quiverVec(P,n,50,'k');
% % % legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
% % % axisGeom; axis manual; camlight headligth;
% % % set(gca,'FontSize',25);

%Extract faces from upper side
[F2,V2] = patchCleanUnused(Fc(~logicSide,:),Vc);

%Merge vertices
[F2,V2] = mergeVertices(F2,V2);

%Get boundary of glenoid en face surface
Eb = patchBoundary(F2,V2);
indBoundary = edgeListToCurve(Eb); %Convert boundary edges to a curve list
indBoundary = indBoundary(1:end-1); %Remove the final point from the boundary to not have closed

%Fill area
%Defining a region and control parameters
regionCell = {V2(indBoundary,:)};
pointSpacing = 1.5; %Desired point spacing
resampleCurveOpt = 1;
interpMethod = 'linear'; %or 'natural'
[backF,backV] = regionTriMesh3D(regionCell,pointSpacing,resampleCurveOpt,interpMethod);

%Plot meshes together
cFigure; hold on;
gpatch(F2,V2,'gw','k');
gpatch(backF,backV,'kw','k');
plotV(V2(indBoundary,:),'r-','LineWidth',2);
axisGeom;

%Join element sets
F_cell{1} = F2; F_cell{2} = backF;
V_cell{1} = V2; V_cell{2} = backV;
[glenoidF,glenoidV] = joinElementSets(F_cell,V_cell);

%Include a colour label for different surfaces
glenoidC = [ones(length(F2),1); ones(length(backF),1)*2];

%Close holes in glenoid surface
[glenodoidF,glenoidV] = triSurfCloseHoles(glenoidF,glenoidV);


%Cleanup
clear F1 F2 V1 V2 F_cell V_cell

[glenoidF,glenoidV] = mergeVertices(glenoidF,glenoidV);

%Visualise the model
cFigure; hold on
title('Imported patch data from STL','fontSize',25);
gpatch(headF,headV,'kw','k');
gpatch(glenoidF,glenoidV,glenoidC,'k');
axisGeom;
%%%%% TODO: identify the appropriate view(x,y) settings for a good view


%% Mesh using tetgen

%Humeral head

%Create tetgen input structure
inputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
inputStruct.Faces = headF; %Boundary faces
inputStruct.Nodes = headV; %Nodes of boundary
inputStruct.regionPoints= getInnerPoint(headF,headV); %Interior points for regions
inputStruct.holePoints = []; %Interior points for holes
inputStruct.regionA = tetVolMeanEst(headF,headV); %Desired tetrahedral volume for each region

%Mesh model
[headMesh] = runTetGen(inputStruct); %Run tetGen

%Get outputs of mesh structure
headVolE = headMesh.elements; %The elements
headVolV = headMesh.nodes; %The vertices or nodes
headVolCE = headMesh.elementMaterialID; %Element material or region id
headVolFb = headMesh.facesBoundary; %The boundary faces
headVolCb = headMesh.boundaryMarker; %The boundary markers

%Visualise humeral head mesh
meshView(headMesh);
title('Tetrahedral Mesh: Humeral Head','FontSize',25);



%Glenoid

%%%%% SELF INTER-SECTING FACES ARE A PROBLEM...!!!!!
stlStruct.solidNames={'intersectingGlenoid'}; %names of parts
stlStruct.solidVertices={glenoidV}; %Vertices
stlStruct.solidFaces={glenoidF}; %Faces
stlStruct.solidNormals={[]}; %Face normals (optional)
export_STL_txt('intersectingGlenoid.stl',stlStruct);

[newGlenoidSTLstruct] = import_STL('intersectingGlenoid.stl');

%Humeral head

%Get single surface faces and vertices
newGlenoidF = newGlenoidSTLstruct.solidFaces{1}; %Faces
newGlenoidV = newGlenoidSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[newGlenoidF,newGlenoidV] = mergeVertices(newGlenoidF,newGlenoidV);


%Create tetgen input structure
inputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options; add -d for self intersecting faces
inputStruct.Faces = newGlenoidF; %Boundary faces
inputStruct.Nodes = newGlenoidV; %Nodes of boundary
inputStruct.regionPoints= getInnerPoint(newGlenoidF,newGlenoidV); %Interior points for regions
inputStruct.holePoints = []; %Interior points for holes
inputStruct.regionA = tetVolMeanEst(newGlenoidF,newGlenoidV); %Desired tetrahedral volume for each region

%Mesh model
[glenoidMesh] = runTetGen(inputStruct); %Run tetGen

%Get outputs of mesh structure
glenoidVolE = glenoidMesh.elements; %The elements
glenoidVolV = glenoidMesh.nodes; %The vertices or nodes
glenoidVolCE = glenoidMesh.elementMaterialID; %Element material or region id
glenoidVolFb = glenoidMesh.facesBoundary; %The boundary faces
glenoidVolCb = glenoidMesh.boundaryMarker; %The boundary markers

%Visualise humeral head mesh
meshView(glenoidMesh);
title('Tetrahedral Mesh: Glenoid','FontSize',25);



%% FEBio details

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

% % % %Contact parameters
% % % contactInitialOffset=0.1;
% % % contactAlg=2;
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

%Get a template with default settings
[febio_spec] = febioStructTemplate;

%febio_spec version
febio_spec.ATTR.version = '2.5';

%Module section
febio_spec.Module.ATTR.type = 'solid';

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type = 'static';
stepStruct.Control.time_steps = numTimeSteps;
stepStruct.Control.step_size = 1/numTimeSteps;
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

% % % febio_spec.Material.material{2}.ATTR.type='Ogden';
% % % febio_spec.Material.material{2}.ATTR.id=1;
% % % febio_spec.Material.material{2}.c1=c1;
% % % febio_spec.Material.material{2}.m1=m1;
% % % febio_spec.Material.material{2}.c2=c1;
% % % febio_spec.Material.material{2}.m2=-m1;
% % % febio_spec.Material.material{2}.k=k;

%Geometry section
% -> Nodes
%%%%% TODO: adapt with additional bodies
febio_spec.Geometry.Nodes{1}.ATTR.name = 'nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id = (1:size(headVolV,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL = headVolV; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat = 1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name = 'HumeralHead'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id = (1:1:size(headVolE,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL = headVolE;

% % % febio_spec.Geometry.Elements{2}.ATTR.type='hex8'; %Element type of this set
% % % febio_spec.Geometry.Elements{2}.ATTR.mat=1; %material index for this set
% % % febio_spec.Geometry.Elements{2}.ATTR.name='Slab'; %Name of the element set
% % % febio_spec.Geometry.Elements{2}.elem.ATTR.id=(1:1:size(E1,1))'; %Element id's
% % % febio_spec.Geometry.Elements{2}.elem.VAL=E1;

% % % % -> NodeSets
% % % febio_spec.Geometry.NodeSet{1}.ATTR.name='bcSupportList';
% % % febio_spec.Geometry.NodeSet{1}.node.ATTR.id=bcSupportList(:);

%%%%% TODO: add surfaces when contact added
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
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL = 45; %prescribed displacement of 45mm

% % % febio_spec.Step{2}.Boundary.rigid_body{1}.ATTR.mat=2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='z';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='Rx';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Ry';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Rz';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.bc='x';
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.ATTR.lc=2;
% % % febio_spec.Step{2}.Boundary.rigid_body{1}.prescribed.VAL=bcPrescribeMagnitudes(1);

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
febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; 1 1];

% % % febio_spec.LoadData.loadcurve{2}.ATTR.id=2;
% % % febio_spec.LoadData.loadcurve{2}.ATTR.type='linear';
% % % febio_spec.LoadData.loadcurve{2}.point.VAL=[0 0; 1 0; 2 1];

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file = febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file = febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data = 'ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim = ',';
febio_spec.Output.logfile.node_data{1}.VAL = 1:size(headVolV,1);

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
[runFlag] = runMonitorFEBio(febioAnalysis);


%% Import FEBio results

%Check for successful run
if runFlag == 1
    
    %Importing nodal displacements from a log file
    %%%%% TODO: comment this better and adjust variable names
    %%%%% TODO: might be a better way to extract results too...
    [time_mat, N_disp_mat,~] = importFEBio_logfile(febioLogFileName_disp); %Nodal displacements
    time_mat = [0; time_mat(:)]; %Time
    N_disp_mat = N_disp_mat(:,2:end,:);
    sizImport = size(N_disp_mat);
    sizImport(3) = sizImport(3)+1;
    N_disp_mat_n = zeros(sizImport);
    N_disp_mat_n(:,:,2:end) = N_disp_mat;
    N_disp_mat = N_disp_mat_n;
    DN = N_disp_mat(:,:,end);
    DN_magnitude = sqrt(sum(DN(:,3).^2,2));
    V_def = headVolV+DN;
    V_DEF = N_disp_mat+repmat(headVolV,[1 1 size(N_disp_mat,3)]);
    X_DEF = V_DEF(:,1,:);
    Y_DEF = V_DEF(:,2,:);
    Z_DEF = V_DEF(:,3,:);
% % %     [CF] = vertexToFaceMeasure(Fb1,DN_magnitude);

    %Plot simulated results using anim8
    
    %Create basic view and store graphics handle to initiate animation
    hf = cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1 = gpatch(headVolFb,V_def,'gw'); %Add graphics object to animate
% % %     hp2=gpatch(E2,V_def,'kw','none',faceAlpha2); %Add graphics object to animate
% % %     gpatch(Fb1,V,0.5*ones(1,3),'none',0.25); %A static graphics object

    axisGeom(gca,25);
% % %     colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;
    
    %Set up animation features
    animStruct.Time = time_mat; %The time vector
    for qt = 1:1:size(N_disp_mat,3) %Loop over time increments
        
        DN = N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude = sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def = headVolV+DN; %Current nodal coordinates
% % %         [CF]=vertexToFaceMeasure(Fb1,DN_magnitude); %Current color data to use

        %Set entries in animation structure
        animStruct.Handles{qt} = [hp1]; %[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt} = {'Vertices'}; %{'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt} = {V_def};%{V_def,CF,V_def}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;

end
