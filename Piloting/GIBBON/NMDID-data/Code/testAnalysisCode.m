%% This code tests out some protocols for processing the 3d surface data
%  in MATLAB from a few landmarks rather than using 3matic features.
%
%  This approach requires GIBBON and geom3d; also the xml_read file

%%%%% TODO: add note and reference to Wu et al. study for coordinate system
%%%%% notation...

%%%%% TODO: function is quite long, could be split up or broken into
%%%%% separate functions to call within

%%%%% TODO: whole scapula mesh doesn't seem to be connected right, get
%%%%% error message about triangles not being connected -- doesn't progress
%%%%% to the glenoid, so there may be some floating or disconnected
%%%%% triangles somewhere else in the mesh

%%%%% TODO: glenoid mesh still a little bumpy in parts from segmentation,
%%%%% should be cleaned up

addpath(genpath(pwd));

warning off

generatePlots = false; %%%set as function input

%% Set-up

%Navigate up to where the data is stored
cd('..');

%Load in scapula and humerus surfaces
[scapulaSTLstruct] = import_STL('Scapula_r.stl');
[humerusSTLstruct] = import_STL('Humerus_r.stl');

%Access the data from the STL structs

%Humeral head

%Get surface faces and vertices
scapulaF = scapulaSTLstruct.solidFaces{1}; %Faces
scapulaV = scapulaSTLstruct.solidVertices{1}; %Vertices
humerusF = humerusSTLstruct.solidFaces{1}; %Faces
humerusV = humerusSTLstruct.solidVertices{1}; %Vertices

%Merge vertices
[scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);
[humerusF,humerusV] = mergeVertices(humerusF,humerusV);

% % % %Convert mm nodes to m
% % % scapulaV = scapulaV / 1000;
% % % humerusV = humerusV / 1000;

%% Load in the landmark data. Convert to m while importing

%Scapula landmarks
scapPoints = [{'AA'},{'AI'},{'TS'},{'DeepGlenoid'},{'SGT'},{'IGT'}];
for pp = 1:length(scapPoints)
    tree = xml_read([scapPoints{pp},'.txt']);
    landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
    clear tree
end
clear pp

%Humeral landmarks
humPoints = [{'GHJC'},{'EL'},{'EM'}];
for pp = 1:length(humPoints)
    tree = xml_read([humPoints{pp},'.txt']);
    landmarks.(char(tree.Point.Name)) = tree.Point.Coordinate;% / 1000;
    clear tree
end
clear pp

%Humeral head sphere
tree = xml_read('HeadSphere.txt');
shapes.headSphere.centre = tree.Sphere.CenterPoint;% / 1000;
shapes.headSphere.radius = tree.Sphere.Radius;% / 1000;
clear tree

%Glenoid plane
tree = xml_read('GlenoidPlane.txt');
planes.glenoid.origin = tree.Plane.Origin;% / 1000;
planes.glenoid.normal = tree.Plane.Normal;
clear tree

%Visualise
%Landmarks
currLandmarks = fieldnames(landmarks);
if generatePlots
    cFigure; hold on
    %Surfaces
    gpatch(scapulaF,scapulaV,'kw','none');
    gpatch(humerusF,humerusV,'kw','none');
    for ff = 1:length(currLandmarks)
        plotV(landmarks.(currLandmarks{ff}),'r.','MarkerSize',15);
    end
    clear ff
    camlight headlight;
    axisGeom;
end

%% Align objects with world coordinate systems
%  Also make deep glenoid point the origin of the WCS

%Convert the exported glenoid plane data to a geom3d plane
glenoidPlane = createPlane(planes.glenoid.origin,planes.glenoid.normal);

%Create the world XY plane
xyPlane = createPlane([0,0,0],[1,0,0],[0,1,0]);

%Create the basis transform between the glenoid plane and world XY plane
worldTransform = createBasisTransform3d(xyPlane,glenoidPlane);

%Transform surfaces
for pp = 1:length(scapulaV)
    scapulaV(pp,:) = transformPoint3d(scapulaV(pp,:),worldTransform);    
end
clear pp
for pp = 1:length(humerusV)
    humerusV(pp,:) = transformPoint3d(humerusV(pp,:),worldTransform);    
end
clear pp

%Transform landmarks
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),worldTransform);
end
clear ff

%Transform shapes
shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,worldTransform);

%Transform planes
planes.glenoid.origin = transformPoint3d(planes.glenoid.origin,worldTransform);
planes.glenoid.normal = transformPoint3d(planes.glenoid.normal,worldTransform);

%Do the secondary rotation to align the Y-axis with the SGT to IGT plane
%Identify the angle need to rotate around the Z-axis to align the SGT and
%IGT points.

%Identify distance between IGT and SGT
pt1 = [landmarks.IGT(1),landmarks.IGT(2)];
pt2 = [landmarks.SGT(1),landmarks.SGT(2)];
dIGT_SGT = sqrt((pt2(2) - pt1(2))^2 + (pt2(1) - pt1(1))^2);

%Identify distance of SGT from IGT along the x-axis
pt3 = [landmarks.IGT(1),landmarks.SGT(2)];
dIGT_x = sqrt((pt3(2) - pt2(2))^2 + (pt3(1) - pt2(1))^2);

%Calculate angle made between pt1 to pt2 and y-axis
rotAng = deg2rad(180) - asin(dIGT_x/dIGT_SGT);

%%%%% TODO: application/calculation of rotation may differ from person to
%%%%% person, relative to whether the SGT is left/right or above/below the
%%%%% IGT -- might expect it to be relatively consistent if scan parameters
%%%%% are all the same, but would change with right vs. left scapula

%Create rotation matrix around Z-axis by specified angle (in radians)
rotMatrix = createRotationOz(rotAng*-1); %-ve for anti-clockwise

%Rotate surfaces
for pp = 1:length(scapulaV)
    scapulaV(pp,:) = transformPoint3d(scapulaV(pp,:),rotMatrix);    
end
clear pp
for pp = 1:length(humerusV)
    humerusV(pp,:) = transformPoint3d(humerusV(pp,:),rotMatrix);    
end
clear pp

%Rotate landmarks
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),rotMatrix);
end
clear ff

%Rotate shapes
shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,rotMatrix);

%Rotate planes
planes.glenoid.origin = transformPoint3d(planes.glenoid.origin,rotMatrix);
planes.glenoid.normal = transformPoint3d(planes.glenoid.normal,rotMatrix);

%Create translation matrix to make deep glenoid point the origin
transMatrix = createTranslation3d([0,0,0] - landmarks.DeepGlenoid);

%Translate surfaces
for pp = 1:length(scapulaV)
    scapulaV(pp,:) = transformPoint3d(scapulaV(pp,:),transMatrix);    
end
clear pp
for pp = 1:length(humerusV)
    humerusV(pp,:) = transformPoint3d(humerusV(pp,:),transMatrix);    
end
clear pp

%Translate landmarks
for ff = 1:length(currLandmarks)
    landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),transMatrix);
end
clear ff

%Translate shapes
shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,transMatrix);

%Translate planes
planes.glenoid.origin = transformPoint3d(planes.glenoid.origin,transMatrix);
planes.glenoid.normal = transformPoint3d(planes.glenoid.normal,transMatrix);

%Visualise rotated
if generatePlots
    cFigure; hold on
    %Surfaces
    gpatch(scapulaF,scapulaV,'kw','none');
    gpatch(humerusF,humerusV,'kw','none',0.2);
    %Landmarks
    currLandmarks = fieldnames(landmarks);
    for ff = 1:length(currLandmarks)
        plotV(landmarks.(currLandmarks{ff}),'r.','MarkerSize',15);
    end
    clear ff
    camlight headlight;
    axisGeom;
end

%% Create scapula coordinate system

%Create line that extends from AA, along the same direction as TS to AA
%This represents the Z-axis of the scapula
%First get direction of TS to AA
TS_AA = createLine3d(landmarks.TS,landmarks.AA);
%Keep the direction but start the line from the AA origin
scapulaCS.Zc = createLine3d(landmarks.AA,TS_AA(4),TS_AA(5),TS_AA(6));

%Define the plane passing through AA, TS and AI
scapBodyPlane = createPlane(landmarks.TS,landmarks.AA,landmarks.AI);
%Get normal of the plane
scapBodyNormal = planeNormal(scapBodyPlane);
%Create line at AA origin with normal of this plane
scapulaCS.Xc = createLine3d(landmarks.AA,scapBodyNormal(1),scapBodyNormal(2),scapBodyNormal(3));

%Take cross product of other two axes
Ytemp = crossProduct3d(scapulaCS.Xc(4:end),scapulaCS.Zc(4:end))*-1;
scapulaCS.Yc = [landmarks.AA,Ytemp];

%Visualise axes
csLabels = [{'Xc'}; {'Yc'}; {'Zc'}];
csColours = [{'r'}; {'g'}; {'b'}];
if generatePlots
    cFigure; hold on
    %Plot scapula
    gpatch(scapulaF,scapulaV,'kw','none',0.2);
    axisGeom; camlight headlight;
    %Plot scapula origin point (use x-axes coordinates --- same across other axes)
    scatter3(scapulaCS.Xc(1),scapulaCS.Xc(2),scapulaCS.Xc(3),1e3,'y.');
    %Plot axes
    for cc = 1:length(csLabels)
        %Get points and vectors
        x = scapulaCS.(csLabels{cc})(1); y = scapulaCS.(csLabels{cc})(2); z = scapulaCS.(csLabels{cc})(3);
        u = scapulaCS.(csLabels{cc})(4); v = scapulaCS.(csLabels{cc})(5); w = scapulaCS.(csLabels{cc})(6);
        %Normalise lengths
        Ln = sqrt(u.^2 + v.^2 + w.^2);
        u = u./Ln; v = v./Ln; w = w./Ln;  %normalize vectors
        MaxLen = 1e-1*1000; %max length preference
        %Set vector length as max length
        u = u*MaxLen;
        v = v*MaxLen;  
        w = w*MaxLen;
        %Plot axes
        quiver3(x,y,z,u,v,w,csColours{cc},'LineWidth',2)
    end
    clear cc
end

%% Create humeral coordinate system

%Create landmark at the midpoint of EL and EM
landmarks.EJC(1) = (landmarks.EL(1) + landmarks.EM(1)) / 2;
landmarks.EJC(2) = (landmarks.EL(2) + landmarks.EM(2)) / 2;
landmarks.EJC(3) = (landmarks.EL(3) + landmarks.EM(3)) / 2;

%Create line connecting GHJC and EJC
Ytemp = createLine3d(landmarks.EJC,landmarks.GHJC);
humerusCS.Yc = [landmarks.GHJC,Ytemp(4:end)];

%Create plane that goes through EL, EM and GHJC
humBodyPlane = createPlane(landmarks.EL,landmarks.EM,landmarks.GHJC);
%Get normal of plane
humBodyNormal = planeNormal(humBodyPlane);
%Create line at GHJC origin with normal of this plane
humerusCS.Xc = createLine3d(landmarks.GHJC,humBodyNormal(1),humBodyNormal(2),humBodyNormal(3));

%Take cross product of other two axes
Ztemp = crossProduct3d(humerusCS.Yc(4:end),humerusCS.Xc(4:end))*-1;
humerusCS.Zc = [landmarks.GHJC,Ztemp];

%Visualise axes
if generatePlots
    cFigure; hold on
    %Plot scapula
    gpatch(humerusF,humerusV,'kw','none',0.2);
    axisGeom; camlight headlight;
    %Plot scapula origin point (use x-axes coordinates --- same across other axes)
    scatter3(humerusCS.Xc(1),humerusCS.Xc(2),humerusCS.Xc(3),1e3,'y.');
    %Plot axes
    for cc = 1:length(csLabels)
        %Get points and vectors
        x = humerusCS.(csLabels{cc})(1); y = humerusCS.(csLabels{cc})(2); z = humerusCS.(csLabels{cc})(3);
        u = humerusCS.(csLabels{cc})(4); v = humerusCS.(csLabels{cc})(5); w = humerusCS.(csLabels{cc})(6);
        %Normalise lengths
        Ln = sqrt(u.^2 + v.^2 + w.^2);
        u = u./Ln; v = v./Ln; w = w./Ln;  %normalize vectors
        MaxLen = 1e-1*1000; %max length preference
        %Set vector length as max length
        u = u*MaxLen;
        v = v*MaxLen;  
        w = w*MaxLen;
        %Plot axes
        quiver3(x,y,z,u,v,w,csColours{cc},'LineWidth',2)
    end
    clear cc
end

%% Align and translate the humerus to a neutral position
%  NOTE: a specified distance of 10mm for the GHJC from the deep glenoid
%  point is specified

%Rotate the humerus and associated landmarks to align with the scapula CS

%Do X-axes alignment first

%Create rotation matrix to align Xc vectors
humerusRotX = createRotationVector3d(humerusCS.Xc(4:end),scapulaCS.Xc(4:end));

%Rotate humerus points
for pp = 1:length(humerusV)
    humerusV(pp,:) = transformPoint3d(humerusV(pp,:),humerusRotX);    
end
clear pp

%Rotate coordinates system points
for cc = 1:length(csLabels)
    humerusCS.(csLabels{cc})(1:3) = transformPoint3d(humerusCS.(csLabels{cc})(1:3),humerusRotX);
    humerusCS.(csLabels{cc})(4:6) = transformPoint3d(humerusCS.(csLabels{cc})(4:6),humerusRotX);
end
clear cc

%Rotate landmarks
landmarks.GHJC = transformPoint3d(landmarks.GHJC,humerusRotX);
landmarks.EL = transformPoint3d(landmarks.EL,humerusRotX);
landmarks.EM = transformPoint3d(landmarks.EM,humerusRotX);
landmarks.EJC = transformPoint3d(landmarks.EJC,humerusRotX);

%Rotate shapes
shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,humerusRotX);

%Do Z-axes alignment second

%Create rotation matrix to align Zc vectors
humerusRotZ = createRotationVector3d(humerusCS.Zc(4:end),scapulaCS.Zc(4:end));

%Rotate humerus points
for pp = 1:length(humerusV)
    humerusV(pp,:) = transformPoint3d(humerusV(pp,:),humerusRotZ);    
end
clear pp

%Rotate coordinates system points
for cc = 1:length(csLabels)
    humerusCS.(csLabels{cc})(1:3) = transformPoint3d(humerusCS.(csLabels{cc})(1:3),humerusRotZ);
    humerusCS.(csLabels{cc})(4:6) = transformPoint3d(humerusCS.(csLabels{cc})(4:6),humerusRotZ);
end
clear cc

%Rotate landmarks
landmarks.GHJC = transformPoint3d(landmarks.GHJC,humerusRotZ);
landmarks.EL = transformPoint3d(landmarks.EL,humerusRotZ);
landmarks.EM = transformPoint3d(landmarks.EM,humerusRotZ);
landmarks.EJC = transformPoint3d(landmarks.EJC,humerusRotZ);

%Rotate shapes
shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,humerusRotZ);

%Translate the GHJC to align with the deep glenoid. The joint centre will
%be placed at a distance of 10mm plus the humeral head sphere radius from
%the deep glenoid point along the world Z-axes. Theoretically this should
%position the edge of the humeral head 10mm from the deep glenoid point.
%Note that the data is in m, not mm

%Identify the point that the centre of the humeral head should sit at
% % % desiredGHJCpoint = [0,0,-0.01 - shapes.headSphere.radius];
desiredGHJCpoint = [0,0,-10 - shapes.headSphere.radius];

%Calculate required shift to get this based on the current GHJC position
humerusTranslate = landmarks.GHJC - desiredGHJCpoint;

%Shift humerus points
for pp = 1:length(humerusV)
    humerusV(pp,:) = humerusV(pp,:) - humerusTranslate;
end
clear pp

%Shift coordinate system points
for cc = 1:length(csLabels)
    humerusCS.(csLabels{cc})(1:3) = humerusCS.(csLabels{cc})(1:3) - humerusTranslate;
end
clear cc

%Shift landmarks
landmarks.GHJC = landmarks.GHJC - humerusTranslate;
landmarks.EL = landmarks.EL - humerusTranslate;
landmarks.EM = landmarks.EM - humerusTranslate;
landmarks.EJC = landmarks.EJC - humerusTranslate;

%Shift shapes
shapes.headSphere.centre = shapes.headSphere.centre - humerusTranslate;

%Visualise
if generatePlots
    cFigure; hold on
    %Bones
    gpatch(scapulaF,scapulaV,'kw','none',0.2);
    gpatch(humerusF,humerusV,'kw','none',0.2);
    axisGeom; camlight headlight;
    %Scapula coordinate system
    for cc = 1:length(csLabels)
        %Get points and vectors
        x = scapulaCS.(csLabels{cc})(1); y = scapulaCS.(csLabels{cc})(2); z = scapulaCS.(csLabels{cc})(3);
        u = scapulaCS.(csLabels{cc})(4); v = scapulaCS.(csLabels{cc})(5); w = scapulaCS.(csLabels{cc})(6);
        %Normalise lengths
        Ln = sqrt(u.^2 + v.^2 + w.^2);
        u = u./Ln; v = v./Ln; w = w./Ln;  %normalize vectors
        MaxLen = 1e-1*1000; %max length preference
        %Set vector length as max length
        u = u*MaxLen;
        v = v*MaxLen;  
        w = w*MaxLen;
        %Plot axes
        quiver3(x,y,z,u,v,w,csColours{cc},'LineWidth',2)
    end
    clear cc
    %Humerus coordinate system
    for cc = 1:length(csLabels)
        %Get points and vectors
        x = humerusCS.(csLabels{cc})(1); y = humerusCS.(csLabels{cc})(2); z = humerusCS.(csLabels{cc})(3);
        u = humerusCS.(csLabels{cc})(4); v = humerusCS.(csLabels{cc})(5); w = humerusCS.(csLabels{cc})(6);
        %Normalise lengths
        Ln = sqrt(u.^2 + v.^2 + w.^2);
        u = u./Ln; v = v./Ln; w = w./Ln;  %normalize vectors
        MaxLen = 1e-1*1000; %max length preference
        %Set vector length as max length
        u = u*MaxLen;
        v = v*MaxLen;  
        w = w*MaxLen;
        %Plot axes
        quiver3(x,y,z,u,v,w,csColours{cc},'LineWidth',2)
    end
    clear cc
end

%% Extract glenoid section from whole scapula

%Slice the scapula 10mm back from the glenoid origin (+ve Z-direction)
%Generate settings for slicing
cutLevel = 0.01*1000; %Set the cut level
snapTolerance = mean(patchEdgeLengths(scapulaF,scapulaV))/100;
n = vecnormalize([0 0 1]); %Normal direction to plane
P = [0 0 cutLevel]; %Point on plane

%Slicing surface (note 3rd color data output is supressed)
[scapulaFc,scapulaVc,~,logicSide,scapulaEc] = triSurfSlice(scapulaF,scapulaV,[],P,n,snapTolerance);

%Visualise slice
if generatePlots
    %Plot split planes
    cFigure; subplot(1,2,1); hold on;
    hp1 = gpatch(scapulaFc(~logicSide,:),scapulaVc,'bw','none',1);
    hp2 = gpatch(scapulaFc(logicSide,:),scapulaVc,'rw','none',1);
    legend([hp1 hp2],{'Surface above plane','Surface below plane'})
    axisGeom; axis manual; camlight headligth;
    colormap gjet;
    set(gca,'FontSize',25);
    %Plot extracted surface and boundary
    subplot(1,2,2); hold on;
    gpatch(scapulaFc(logicSide,:),scapulaVc,'w','none',1);
    gpatch(scapulaFc(~logicSide,:),scapulaVc,'w','none',0.25);
    hp1=gpatch(scapulaEc,scapulaVc,'none','b',1,3);
    hp2=quiverVec(P,n,0.05,'k');
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

%% Create a humeral head mesh from the sphere

%Create sphere surface
[headF,headV] = geoSphere(3,shapes.headSphere.radius);

%Shift the head to the GHJC position
%Can just subtract the GHJC position as the sphere gets centred at 0,0,0
for hh = 1:length(headV)
    headV(hh,:) = shapes.headSphere.centre - headV(hh,:);
end
clear hh

%Translate humeral head surface to sit just off the glenoid without making
%direct contact
headV(:,3) = headV(:,3) + 9.5; %%%this is relative to the humeral position

%Visualise extracted surfaces together
if generatePlots
    cFigure; hold on
    gpatch(glenoidF,glenoidV,'kw','none',0.5);
    gpatch(headF,headV,'kw','none',0.5);
    axisGeom;
end

%% TODO: Mesh refinement...

%Remeshing too...

%do a basic remesh for now of the glenoid
%in it's current format there seems to be too many nodes for the febio spec
%text file to be created in computer memory
[Fb,Vb] = triRemeshLabel(glenoidF,glenoidV,1.0);

% % % cFigure;
% % % gpatch(Fb,Vb);
% % % axisGeom;

glenoidF = Fb;
glenoidV = Vb;

%% Create volumetric mesh of glenoid

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

%% Test adaption of the sphere sliding demo

%% Set-up for simulation

%Face normals need to be flipped on both surfaces for some reason
headF = fliplr(headF);
glenoidVolFb = fliplr(glenoidVolFb);

%Joining node sets
V = [glenoidVolV; headV; ]; %Combined node sets
headF = headF + size(glenoidVolV,1); %Fixed element indices for head surface

%Define boundary conditions

%Identify back facing surface on glenoid to support

%Extract normals of faces
[glenoidN] = patchNormal(glenoidVolFb,glenoidVolV);

%Find face normals that are pointing exactly along the Z-axis
%Also check that the X and Y direction of these normals is quite small. It
%doesn't seem to ever be EXACTLY zero
logicRigid = logical(glenoidN(:,1) < 1e-10) & ...
    logical(glenoidN(:,2) < 1e-10) &logical(glenoidN(:,3) == 1);
Fr = glenoidVolFb(logicRigid,:);
bcSupportList = unique(Fr(:));

%Visualize BC's
if generatePlots
    cFigure; hold on;
    title('Boundary Conditions: Support Nodes');
    gpatch(glenoidVolFb,V,'kw','none',0.3); 
    plotV(V(bcSupportList,:),'k.','MarkerSize',10);
    axisGeom;
    camlight headlight;
end

%Define contact surfaces

%Rigid master surface of the sphere
F_contact_master = headF;

%Elastic slave surface of the glenoid
%Use faces that aren't on the rear of the glenoid
%%%%% TODO: investigate capacity of isolating faces on the glenoid surface
%%%%% in a similar manner to the connected vertebrae GIBBON demp
F_contact_slave = glenoidVolFb(~logicRigid,:);

% Plotting surface models and contact faces
if generatePlots
    cFigure; hold on;
    title('Contact sets and normal directions');
    gpatch(glenoidVolFb,V,'kw','none',0.3); 
    hl(1) = gpatch(F_contact_master,V,'g','k',1); 
    patchNormPlot(F_contact_master,V);
    hl(2) = gpatch(F_contact_slave,V,'b','k',1);
    patchNormPlot(F_contact_slave,V);
    legend(hl,{'Master','Slave'});
    axisGeom;
    camlight headlight;
    drawnow;
end

%% Define the FEBio input structure

%%%%% TODO: adapt from basic parameters of the sphere sliding demo

%Get a template with default settings 
[febio_spec] = febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version = '2.5'; 

%Module section
febio_spec.Module.ATTR.type = 'solid'; 

%Create control structure for use by all steps
stepStruct.Control.analysis.ATTR.type = 'static';
stepStruct.Control.time_steps = 10;
stepStruct.Control.step_size = 1/10;
stepStruct.Control.time_stepper.dtmin = (1/10)/100;
stepStruct.Control.time_stepper.dtmax = 1/10; 
stepStruct.Control.time_stepper.max_retries = 5;
stepStruct.Control.time_stepper.opt_iter = 10;
stepStruct.Control.max_refs = 25;
stepStruct.Control.max_ups = 0;

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

% % % %Sphere sliding material
febio_spec.Material.material{1}.ATTR.type='Ogden';
febio_spec.Material.material{1}.ATTR.id=1;
febio_spec.Material.material{1}.c1=1e-3;
febio_spec.Material.material{1}.m1=8;
febio_spec.Material.material{1}.c2=1e-3;
febio_spec.Material.material{1}.m2=-8;
febio_spec.Material.material{1}.k=1e-3*1e2;

%Head material
febio_spec.Material.material{2}.ATTR.type = 'rigid body';
febio_spec.Material.material{2}.ATTR.name = 'rigidHead';
febio_spec.Material.material{2}.ATTR.id = 2;
febio_spec.Material.material{2}.density = 1;
% febio_spec.Material.material{2}.center_of_mass = center_of_mass;

%Geometry section

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

% -> NodeSets
%Nodes to support back of glenoid
febio_spec.Geometry.NodeSet{1}.ATTR.name = 'bcSupportList';
febio_spec.Geometry.NodeSet{1}.node.ATTR.id = bcSupportList(:);

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

%Fix the glenoid support nodes
febio_spec.Boundary.fix{1}.ATTR.bc = 'x';
febio_spec.Boundary.fix{1}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{2}.ATTR.bc = 'y';
febio_spec.Boundary.fix{2}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;
febio_spec.Boundary.fix{3}.ATTR.bc = 'z';
febio_spec.Boundary.fix{3}.ATTR.node_set = febio_spec.Geometry.NodeSet{1}.ATTR.name;

% -> Prescribed boundary conditions on the rigid body

%Initial displacement
%Note that distance to contact surface is 0.01m from the glenoid, so we add
%a tiny little bit to ensure contact
febio_spec.Step{1}.Boundary.rigid_body{1}.ATTR.mat = 2;
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{1}.ATTR.bc = 'x';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{2}.ATTR.bc = 'y';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{3}.ATTR.bc = 'Rx';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{4}.ATTR.bc = 'Ry';
febio_spec.Step{1}.Boundary.rigid_body{1}.fixed{5}.ATTR.bc = 'Rz';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.bc = 'z';
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.ATTR.lc = 1;
febio_spec.Step{1}.Boundary.rigid_body{1}.prescribed.VAL = 2.5; %10-0.5;

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

%Contact section
contactType = 'facet-to-facet-sliding';
pointSpacings = mean(patchEdgeLengths(glenoidVolFb,V)); %value used earlier for basic remesh
switch contactType
    case 'sticky'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sticky';
        febio_spec.Contact.contact{1}.penalty=100;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.snap_tol=0;
        febio_spec.Contact.contact{1}.max_traction=0;
        febio_spec.Contact.contact{1}.search_tolerance=0.1;
    case 'facet-to-facet-sliding'

        %sphere sliding parameters
        
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='facet-to-facet sliding';
        
        %%% 1e+07 crashes at first time step; same with 1e+04 (penalty)
        
        %%% Augmented Lagragian method is likely the better solution given
        %%% the large penalties seem to cause errors in process.
        
        %%% Seems like there is two approaches with laugon = 1, the first
        %%% being to use a convergence tolerance with the 'tolerance'
        %%% parameter - a tolerance of 0.01 will converge when there is less
        %%% than a 1% change in the L2 norm of the augmented Lagrangian
        %%% multiplier vector between successive augmentations. A tolerance
        %%% of 0.01 seems to start struggling when contact occurs. Same
        %%% thing seemingly happens with 0.1 tolerance. This can be set to
        %%% zero and the gaptol parameter can be used. In this case, the
        %%% iterations will terminate (i.e. progress) when the gap norm is
        %%% less than the user specified gap tolerance. The gaptol
        %%% parameter is an absolute value, so this parameter depends on
        %%% the dimensions of the model. There are cases when a gap
        %%% tolerance can't be reached simply due to model geometry, so the
        %%% problem will never converge. A gaptol of 1 struggles even
        %%% earlier than the tolerance approach. Lower gaptol of 0.001 goes
        %%% even worse. gaptol of 10 gets to the same 0.4 point as the
        %%% tolerance approach before it starts to diverge and reformat the
        %%% stiffness matrix - but did eventually push through to converge
        %%% at 0.4 and onwards --- it takes a bit of time though due to
        %%% all the reformations and augmentations. The time steps it
        %%% reduces to are still similar to the sphere sliding example
        %%% though when using the same material parameter set. This method
        %%% does take a long time, but it's possibly due to the amount of
        %%% deformation that the squishy ogden material produces changing
        %%% the interface between the objects with the amount of
        %%% displacement being applied (this was for 2.5mm into glenoid).
        %%% Cancelled this but it was basically about to finish (~0.999
        %%% through). There is a little bit of penetration through the
        %%% glenoid surface - which could probably be reduced by lowering
        %%% the gaptol? Nonethless, a penalty of 100, laugon=1,
        %%% tolerance=0, gaptol = 10, search_tol=0.01 works with the Ogden
        %%% material, but with some penetration...
        %%%
        %%% All of the above was done with a penalty of 100, not sure
        %%% whether this is considered in the laugon method? It is, but
        %%% scales the Lagrange multiplier increment, rather than the gap.
        %%% Lower values are suggested when using laugon, so this might not
        %%% need to be so high? Eqaully, this could be tried in conjunction
        %%% with the auto-penalty flag too to set it a little better.
        %%%
        %%% The sliding-elastic method is also recommended in the manual
        %%% for compression focused problems - converges better but is
        %%% non-symemetric so takes longer. There are different parameters
        %%% to this, but the sphere tube slide body force demo uses
        %%% sliding-elastic contact so would serve as a similarly relevant
        %%% example.
        
        febio_spec.Contact.contact{1}.penalty=100; %upped from 100
        febio_spec.Contact.contact{1}.auto_penalty=0; %changed
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=1;
        febio_spec.Contact.contact{1}.tolerance=0; %convergence tolerance for laug method
        febio_spec.Contact.contact{1}.gaptol=10;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
    
    case 'sliding-elastic'
        
        %Start by using laugon = 1, gaptol = 10, penalty = 100 as per facet
        %to facet that worked OK, but took a long time...
        
        %%% Current implementation will use the default 1.0 for ktmult
        %%% (i.e. the tangential stiffness multiplier) --- this is a
        %%% scaling factor seemingly for friction, which at 1.0 wouldn't
        %%% scale, and given the contact is frictionless probably doesn't
        %%% matter
        
        %%% Same with the segmentation updates defaulting to off --- which
        %%% means segment updates will occur across the whole simulation
        
        %%% Same with the smooth Lagrangian augmentation defaulting to off
        
        %%% THIS DID NOT GO WELL...a lot of negative jacobian
        %%% errors...potentially need to try with non symmetric stiffness
        %%% formulation in the control section
        
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
        febio_spec.Contact.contact{1}.two_pass=1; %different to f2f
        febio_spec.Contact.contact{1}.laugon=1;
        febio_spec.Contact.contact{1}.tolerance=0;
        febio_spec.Contact.contact{1}.gaptol=10;
        febio_spec.Contact.contact{1}.minaug=1;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2; %pointSpacings/10;
        febio_spec.Contact.contact{1}.symmetric_stiffness=0;
        febio_spec.Contact.contact{1}.auto_penalty=0;
        febio_spec.Contact.contact{1}.penalty=100;
        febio_spec.Contact.contact{1}.fric_coeff=0; %no friction
        
    case 'sliding_with_gaps'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding_with_gaps';
        febio_spec.Contact.contact{1}.penalty=100;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.minaug=0;
        febio_spec.Contact.contact{1}.maxaug=10;
        febio_spec.Contact.contact{1}.fric_coeff=0;
        febio_spec.Contact.contact{1}.fric_penalty=0;
        febio_spec.Contact.contact{1}.ktmult=1;
        febio_spec.Contact.contact{1}.seg_up=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
    case 'sliding2'
        febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
        febio_spec.Contact.contact{1}.ATTR.type='sliding2';
        febio_spec.Contact.contact{1}.penalty=30;
        febio_spec.Contact.contact{1}.auto_penalty=1;
        febio_spec.Contact.contact{1}.two_pass=0;
        febio_spec.Contact.contact{1}.laugon=0;
        febio_spec.Contact.contact{1}.tolerance=0.1;
        febio_spec.Contact.contact{1}.gaptol=0;
        febio_spec.Contact.contact{1}.symmetric_stiffness=0;
        febio_spec.Contact.contact{1}.search_tol=0.01;
        febio_spec.Contact.contact{1}.search_radius=mean(pointSpacings)/2;
end

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
febioFebFileNamePart='sphereSliding_adaptation';
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

%% test import of data and display

if runFlag==1 %i.e. a succesful run
    
    %%%%% TODO: still think this import process can be tidied up...
    
    % Importing nodal displacements from a log file
    [time_mat, N_disp_mat,~]=importFEBio_logfile(febioLogFileName_disp); %Nodal displacements
    time_mat=[0; time_mat(:)]; %Time

    N_disp_mat=N_disp_mat(:,2:end,:);
    sizImport=size(N_disp_mat);
    sizImport(3)=sizImport(3)+1;
    N_disp_mat_n=zeros(sizImport);
    N_disp_mat_n(:,:,2:end)=N_disp_mat;
    N_disp_mat=N_disp_mat_n;
    DN=N_disp_mat(:,:,end);
    DN_magnitude=sqrt(sum(DN(:,3).^2,2));
    V_def=V+DN;
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
    X_DEF=V_DEF(:,1,:);
    Y_DEF=V_DEF(:,2,:);
    Z_DEF=V_DEF(:,3,:);
    [CF]=vertexToFaceMeasure(glenoidVolFb,DN_magnitude);
    
    % Create basic view and store graphics handle to initiate animation
    hf=cFigure; %Open figure
    gtitle([febioFebFileNamePart,': Press play to animate']);
    hp1=gpatch(glenoidVolFb,V_def,CF,'k',1); %Add graphics object to animate
    hp2=gpatch(headF,V_def,'kw','none',0.3); %Add graphics object to animate
    gpatch(glenoidVolFb,V,0.5*ones(1,3),'none',0.25); %A static graphics object

    axisGeom;
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);
    axis([min(X_DEF(:)) max(X_DEF(:)) min(Y_DEF(:)) max(Y_DEF(:)) min(Z_DEF(:)) max(Z_DEF(:))]);
    camlight headlight;

    % Set up animation features
    animStruct.Time=time_mat; %The time vector
    for qt=1:1:size(N_disp_mat,3) %Loop over time increments
        DN=N_disp_mat(:,:,qt); %Current displacement
        DN_magnitude=sqrt(sum(DN.^2,2)); %Current displacement magnitude
        V_def=V+DN; %Current nodal coordinates
        [CF]=vertexToFaceMeasure(glenoidVolFb,DN_magnitude); %Current color data to use
        %%%%% plots colour for displacement magnitude

        %Set entries in animation structure
        animStruct.Handles{qt}=[hp1 hp1 hp2]; %Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData','Vertices'}; %Properties of objects to animate
        animStruct.Set{qt}={V_def,CF,V_def}; %Property values for to set in order to animate
    end
    anim8(hf,animStruct); %Initiate animation feature
    drawnow;
    
    
    
end















%%
