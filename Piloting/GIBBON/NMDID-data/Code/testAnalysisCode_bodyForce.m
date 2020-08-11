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

% % % %Create sphere surface
% % % [headF,headV] = geoSphere(3,shapes.headSphere.radius);
% % % 
% % % %Shift the head to the GHJC position
% % % %Can just subtract the GHJC position as the sphere gets centred at 0,0,0
% % % for hh = 1:length(headV)
% % %     headV(hh,:) = shapes.headSphere.centre - headV(hh,:);
% % % end
% % % clear hh
% % % 
% % % %Translate humeral head surface to sit just off the glenoid without making
% % % %direct contact
% % % headV(:,3) = headV(:,3) + 9.5; %%%this is relative to the humeral position
% % % 
% % % %Visualise extracted surfaces together
% % % if generatePlots
% % %     cFigure; hold on
% % %     gpatch(glenoidF,glenoidV,'kw','none',0.5);
% % %     gpatch(headF,headV,'kw','none',0.5);
% % %     axisGeom;
% % % end

%% Create solid hexahedral sphere mesh
% Currently using spher tube sliding body force example parameters

%Control settings
pointSpacing = 1.0;
optionStruct.sphereRadius = shapes.headSphere.radius;
optionStruct.coreRadius = shapes.headSphere.radius.*0.75; %won't need this...
optionStruct.numElementsCore=ceil((shapes.headSphere.radius/2)/pointSpacing);
optionStruct.numElementsMantel=ceil((shapes.headSphere.radius-optionStruct.coreRadius)/(2*pointSpacing));
optionStruct.makeHollow=0;
optionStruct.outputStructType=2;

%Creating sphere
[meshOutput] = hexMeshSphere(optionStruct);

% Access model element and patch data
Fb_blob = meshOutput.facesBoundary;
Cb_blob = meshOutput.boundaryMarker;
V_blob = meshOutput.nodes;
E_blob = meshOutput.elements;

if generatePlots
    hFig=cFigure;
    subplot(1,2,1); hold on;
    title('Boundary surfaces');
    gpatch(Fb_blob,V_blob,Cb_blob,'k',0.5);
    axisGeom;
    colormap(gjet); icolorbar;
    camlight headlight;

    hs=subplot(1,2,2); hold on;
    title('Cut view of solid mesh');
    optionStruct.hFig=[hFig hs];
    gpatch(Fb_blob,V_blob,'kw','none',0.25);
    meshView(meshOutput,optionStruct);
    axisGeom(gca,fontSize);
    drawnow;
end

%Shift the head to the GHJC position
%Can just subtract the GHJC position as the sphere gets centred at 0,0,0
for hh = 1:length(V_blob)
    V_blob(hh,:) = shapes.headSphere.centre - V_blob(hh,:);
end
clear hh

%Translate humeral head surface to sit just off the glenoid without making
%direct contact
V_blob(:,3) = V_blob(:,3) + 9.5; %%%this is relative to the humeral position


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

%Visualise together
if generatePlots
    cFigure; hold on;
    gpatch(Fb_blob,V_blob,Cb_blob,'k',0.5);
    gpatch(glenoidF,glenoidV,'kw','k',0.5);
    axisGeom;
    camlight headlight;
end

%% Create volumetric mesh of glenoid

% % % %Create tetgen input structure
% % % glenoidInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
% % % glenoidInputStruct.Faces = glenoidF; %Boundary faces
% % % glenoidInputStruct.Nodes = glenoidV; %Nodes of boundary
% % % glenoidInputStruct.regionPoints= getInnerPoint(glenoidF,glenoidV); %Interior points for regions
% % % glenoidInputStruct.holePoints = []; %Interior points for holes
% % % glenoidInputStruct.regionA = tetVolMeanEst(glenoidF,glenoidV); %Desired tetrahedral volume for each region
% % % 
% % % %Mesh model
% % % [glenoidMesh] = runTetGen(glenoidInputStruct); %Run tetGen
% % % 
% % % %Get outputs of mesh structure
% % % glenoidVolE = glenoidMesh.elements; %The elements
% % % glenoidVolV = glenoidMesh.nodes; %The vertices or nodes
% % % glenoidVolCE = glenoidMesh.elementMaterialID; %Element material or region id
% % % glenoidVolFb = glenoidMesh.facesBoundary; %The boundary faces
% % % glenoidVolCb = glenoidMesh.boundaryMarker; %The boundary markers
% % % 
% % % %Visualise glenoid mesh
% % % if generatePlots
% % %     meshView(glenoidMesh);
% % %     title('Tetrahedral Mesh: Glenoid','FontSize',25);
% % % end

%% Test adaption of the sphere sliding demo

%% Set-up for simulation

%Face normals need to be flipped on both surfaces for some reason
Fb_blob = fliplr(Fb_blob);
% % glenoidVolFb = fliplr(glenoidVolFb);

%Joining node sets
V = [V_blob; glenoidV; ]; %Combined node sets
glenoidF = glenoidF + size(V_blob,1); %Fixed element indices for head surface

%Visualise joined node sets
if generatePlots
    cFigure; hold on;
    gtitle('Model components');
    hl(1)=gpatch(Fb_blob,V,'rw','k',0.8);
    hl(2)=gpatch(glenoidF,V,'bw','k',0.8);
    legend(hl,{'Blob','Glenoid'}); clear hl;
    axisGeom;
    camlight headlight;
    drawnow;
end

%Define contact surfaces
F_contact_blob = Fb_blob;

%Visualise contact surfaces
if generatePlots
    cFigure; hold on;
    title('The contact pair');
    hl(1)=gpatch(glenoidF,V,'rw','k',1);
    patchNormPlot(glenoidF,V);
    hl(2)=gpatch(F_contact_blob,V,'gw','k',1);
    patchNormPlot(F_contact_blob,V);
    legend(hl,{'Master','Slave'}); clear hl;
    axisGeom;
    camlight headlight;
    drawnow;
end

%% Define the FEBio input structure

% Material parameter set
materialType = 1;

c1=1e-3; %Shear-modulus-like parameter
m1=2; %Material parameter setting degree of non-linearity
k_factor=10; %Bulk modulus factor
k=c1*k_factor; %Bulk modulus
g1=0.5; %Viscoelastic QLV proportional coefficient
t1=12; %Viscoelastic QLV time coefficient
materialDensity=1e-9*100000; %Density

% FEA control settings
timeRamp=1; %Time for ramping up body force from 0-100%
timeTotal=5; %Total simulation time
stepSizeDesired=0.25; %Desired time increment size
numTimeSteps=round(timeTotal/stepSizeDesired); %Number of time steps desired
max_refs=25; %Max reforms
max_ups=0; %Set to zero to use full-Newton iterations
opt_iter=10; %Optimum number of iterations
max_retries=10; %Maximum number of retires
step_size=timeTotal/numTimeSteps;
dtmin=(timeTotal/numTimeSteps)/100; %Minimum time step size
dtmax=timeTotal/(numTimeSteps); %Maximum time step size
symmetric_stiffness=0;
min_residual=1e-20;
analysisType='dynamic';

%Contact parameters
two_pass=1;
contactPenalty=25;
fric_coeff=0.1;

bodyLoadMagnitude=5;

%%%%% TODO: adapt from basic parameters of the sphere tube demo

%Get a template with default settings 
[febio_spec] = febioStructTemplate;

%febio_spec version 
febio_spec.ATTR.version = '2.5'; 

%Module section
febio_spec.Module.ATTR.type = 'solid'; 

%Control section
febio_spec.Control.analysis.ATTR.type=analysisType;
febio_spec.Control.time_steps=numTimeSteps;
febio_spec.Control.step_size=step_size;
febio_spec.Control.time_stepper.dtmin=dtmin;
febio_spec.Control.time_stepper.dtmax=dtmax;
febio_spec.Control.time_stepper.max_retries=max_retries;
febio_spec.Control.time_stepper.opt_iter=opt_iter;
febio_spec.Control.max_refs=max_refs;
febio_spec.Control.max_ups=max_ups;
febio_spec.Control.symmetric_stiffness=symmetric_stiffness;
febio_spec.Control.min_residual=min_residual;

%Material section
switch materialType
    case 0
        febio_spec.Material.material{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.c1=c1;
        febio_spec.Material.material{1}.m1=m1;
        febio_spec.Material.material{1}.c2=c1;
        febio_spec.Material.material{1}.m2=-m1;
        febio_spec.Material.material{1}.cp=k;
        febio_spec.Material.material{1}.density=materialDensity;
    case 1
        %Viscoelastic part
        febio_spec.Material.material{1}.ATTR.type='viscoelastic';
        febio_spec.Material.material{1}.ATTR.Name='Block_material';
        febio_spec.Material.material{1}.ATTR.id=1;
        febio_spec.Material.material{1}.g1=g1;
        febio_spec.Material.material{1}.t1=t1;
        febio_spec.Material.material{1}.density=materialDensity;

        %Elastic part
        febio_spec.Material.material{1}.elastic{1}.ATTR.type='Ogden unconstrained';
        febio_spec.Material.material{1}.elastic{1}.c1=c1;
        febio_spec.Material.material{1}.elastic{1}.m1=m1;
        febio_spec.Material.material{1}.elastic{1}.c2=c1;
        febio_spec.Material.material{1}.elastic{1}.m2=-m1;
        febio_spec.Material.material{1}.elastic{1}.cp=k;
        febio_spec.Material.material{1}.elastic{1}.density=materialDensity;
end
febio_spec.Material.material{2}.ATTR.type='rigid body';
febio_spec.Material.material{2}.ATTR.id=2;
febio_spec.Material.material{2}.density=1;
febio_spec.Material.material{2}.center_of_mass=mean(glenoidV,1);

%Geometry section
% -> Nodes
febio_spec.Geometry.Nodes{1}.ATTR.name='nodeSet_all'; %The node set name
febio_spec.Geometry.Nodes{1}.node.ATTR.id=(1:size(V,1))'; %The node id's
febio_spec.Geometry.Nodes{1}.node.VAL=V; %The nodel coordinates

% -> Elements
febio_spec.Geometry.Elements{1}.ATTR.type='hex8'; %Element type of this set
febio_spec.Geometry.Elements{1}.ATTR.mat=1; %material index for this set
febio_spec.Geometry.Elements{1}.ATTR.name='Blob'; %Name of the element set
febio_spec.Geometry.Elements{1}.elem.ATTR.id=(1:1:size(E_blob,1))'; %Element id's
febio_spec.Geometry.Elements{1}.elem.VAL=E_blob;

febio_spec.Geometry.Elements{2}.ATTR.type='tri3'; %Element type of this set
febio_spec.Geometry.Elements{2}.ATTR.mat=2; %material index for this set
febio_spec.Geometry.Elements{2}.ATTR.name='glenoid'; %Name of the element set
febio_spec.Geometry.Elements{2}.elem.ATTR.id=size(E_blob,1)+(1:1:size(glenoidF,1))'; %Element id's
febio_spec.Geometry.Elements{2}.elem.VAL=glenoidF;

% % -> Surfaces
febio_spec.Geometry.Surface{1}.ATTR.name='contact_master1';
febio_spec.Geometry.Surface{1}.tri3.ATTR.lid=(1:1:size(glenoidF,1))';
febio_spec.Geometry.Surface{1}.tri3.VAL=glenoidF;

febio_spec.Geometry.Surface{2}.ATTR.name='contact_slave1';
febio_spec.Geometry.Surface{2}.quad4.ATTR.lid=(1:1:size(F_contact_blob,1))';
febio_spec.Geometry.Surface{2}.quad4.VAL=F_contact_blob;

% -> Surface pairs
febio_spec.Geometry.SurfacePair{1}.ATTR.name='Contact1_tube_blob';
febio_spec.Geometry.SurfacePair{1}.master.ATTR.surface=febio_spec.Geometry.Surface{1}.ATTR.name;
febio_spec.Geometry.SurfacePair{1}.slave.ATTR.surface=febio_spec.Geometry.Surface{2}.ATTR.name;

%Loads
febio_spec.Loads.body_load{1}.ATTR.type='const';
febio_spec.Loads.body_load{1}.z.VAL=bodyLoadMagnitude;
febio_spec.Loads.body_load{1}.z.ATTR.lc=1;

%LoadData section
febio_spec.LoadData.loadcurve{1}.ATTR.id=1;
febio_spec.LoadData.loadcurve{1}.ATTR.type='linear';
febio_spec.LoadData.loadcurve{1}.point.VAL=[0 0; timeRamp 1];

% -> Prescribed boundary conditions on the rigid body
febio_spec.Boundary.rigid_body{1}.ATTR.mat=2;
febio_spec.Boundary.rigid_body{1}.fixed{1}.ATTR.bc='x';
febio_spec.Boundary.rigid_body{1}.fixed{2}.ATTR.bc='y';
febio_spec.Boundary.rigid_body{1}.fixed{3}.ATTR.bc='z';
febio_spec.Boundary.rigid_body{1}.fixed{4}.ATTR.bc='Rx';
febio_spec.Boundary.rigid_body{1}.fixed{5}.ATTR.bc='Ry';
febio_spec.Boundary.rigid_body{1}.fixed{6}.ATTR.bc='Rz';

%Contact section
febio_spec.Contact.contact{1}.ATTR.surface_pair=febio_spec.Geometry.SurfacePair{1}.ATTR.name;
febio_spec.Contact.contact{1}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{1}.two_pass=two_pass;
febio_spec.Contact.contact{1}.laugon=0;
febio_spec.Contact.contact{1}.tolerance=0.2;
febio_spec.Contact.contact{1}.gaptol=0;
febio_spec.Contact.contact{1}.minaug=1;
febio_spec.Contact.contact{1}.maxaug=10;
febio_spec.Contact.contact{1}.search_tol=0.01;
febio_spec.Contact.contact{1}.search_radius=pointSpacing/10;
febio_spec.Contact.contact{1}.symmetric_stiffness=0;
febio_spec.Contact.contact{1}.auto_penalty=1;
febio_spec.Contact.contact{1}.penalty=contactPenalty;
febio_spec.Contact.contact{1}.fric_coeff=fric_coeff;

% Defining file names
febioFebFileNamePart='tubeSlide_adaptation';
febioFebFileName=[febioFebFileNamePart,'.feb']; %FEB file name
febioLogFileName=[febioFebFileNamePart,'.txt']; %FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; %Log file name for exporting displacement
% febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; %Log file name for exporting force

%Output section
% -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

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
febioAnalysis.runMode='internal';%'internal';
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
