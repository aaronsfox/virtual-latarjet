function [scapulaOutput,humerusOutput,glenoidOutput,headOutput,...
    scapulaCS,humerusCS,landmarks,shapes] = importAndCreateSurfaces(participantDir,generatePlots)

%% This function serves to import in and create the base surface system for
%  running the FEA analysis of the humeral head against the glenoid.
%
%  Inputs:
%
%   participantDir      directory for the current participant containing
%   generatePlots       flag whether to generate figures from the
%                       processing throughout the function (default = false)
%
%  Outputs
%
%   scapulaOutput       surface mesh data for scapula
%   humerusOutput       surface mesh data for humerus
%   glenoidOutput       volumetric mesh data for glenoid
%   headOutput          volumetric mesh data for head hemisphere
%   scapulaCS           coordinate system axes for scapula
%   humerusCS           coordinate system axes for humerus
%   landmarks           imported landmark points on scapula/humerus
%   shapes              imported shapes data for scapula/humerus (i.e. humeral head radius)
%
%  This function, like parts of the main code, uses elements of the GIBBON
%  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
%  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).

    %% Check inputs

    %Participant directory
    if nargin < 1
        error('The participant directory where the surfaces are stored is required as a function input')
    end

    %Generate plots flag
    if nargin < 2
        generatePlots = false;
    else
        %Check whether it is a logical operator
        if ~islogical(generatePlots)
            error('generatePlots function inpur must be a logical of true or false')
        end
    end

    %% Set-up

    %Get starting directory to return to
    home_dir = cd;

    %Navigate to participant directory
    cd(participantDir);

    %Identify whether right or left limb data is present. Search for the
    %humerus .stl file and check it's trailing limb label
    checkFile = dir('Humerus_*.stl');
    if size(checkFile,1) > 1
        error('More than one humerus surface file identified. Check file labels.');
    end
    %Identify the limb based on file name
    limbCheck = strsplit(checkFile.name,'_');
    if strcmp(limbCheck{2}(1),'r')
        limb = 'r';
    elseif strcmp(limbCheck{2}(1),'l')
        limb = 'l';
    else
        error('Limb can''t be identified from humerus surface file name. Check file labels.');
    end

    %Start the waitbar
    wbar = waitbar(0/9,'Importing surfaces...');
    
    %% Load and clean up surface files
    
    %Update the waitbar
    waitbar(0/9,wbar,'Importing scapula...');
    
    %Scapula
    [scapulaSTLstruct] = import_STL(['Scapula_',limb,'.stl']);
    scapulaF = scapulaSTLstruct.solidFaces{1}; %Faces
    scapulaV = scapulaSTLstruct.solidVertices{1}; %Vertices
    [scapulaF,scapulaV] = mergeVertices(scapulaF,scapulaV);
    
    %Update the waitbar
    waitbar(0.33/9,wbar,'Importing humerus...');
    
    %Humerus
    [humerusSTLstruct] = import_STL(['Humerus_',limb,'.stl']);
    humerusF = humerusSTLstruct.solidFaces{1}; %Faces
    humerusV = humerusSTLstruct.solidVertices{1}; %Vertices
    [humerusF,humerusV] = mergeVertices(humerusF,humerusV);
    
    %Update the waitbar
    waitbar(0.66/9,wbar,'Importing articulating surface...');
    
    %Humeral articulating surface
    [headSTLstruct] = import_STL(['HumerusArticulating_',limb,'.stl']);
    headF = headSTLstruct.solidFaces{1}; %Faces
    headV = headSTLstruct.solidVertices{1}; %Vertices
    [headF,headV] = mergeVertices(headF,headV);

    % % % %Convert mm nodes to m
    % % % scapulaV = scapulaV / 1000;
    % % % humerusV = humerusV / 1000;
    % % % headV = headV / 1000;

    %% Load in the landmark data.
    %  NOTE: uses xml_read supplementary function
    %  NOTE: commented out method for converting from mm to m

    %Update waitbar
    waitbar(1/9,wbar,'Importing landmarks...');
    
    %Scapula landmarks
    scapPoints = [{'AA'},{'AI'},{'TS'},{'DeepGlenoid'},{'SGT'},{'IGT'},...
        {'AntGlenoid'},{'PostGlenoid'},{'SupGlenoid'},{'InfGlenoid'}];
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
    %Get current landmarks
    currLandmarks = fieldnames(landmarks);
    if generatePlots
        cFigure; hold on
        %Surfaces
        gpatch(scapulaF,scapulaV,'kw','none',0.8);
        gpatch(humerusF,humerusV,'kw','none',0.8);
        gpatch(headF,headV,'rw','none',1);
        for ff = 1:length(currLandmarks)
            plotV(landmarks.(currLandmarks{ff}),'g.','MarkerSize',15);
        end
        clear ff
        axisGeom;
        title('Raw Imported Surfaces and Landmarks');
    end

    %% Align objects with world coordinate systems
    %  As part of this process the deep glenoid also becomes the origin
    %  point (i.e. [0,0,0])
    
    %Update waitbar
    waitbar(2/9,wbar,'Aligning surfaces...');
    
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
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),worldTransform);    
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
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),rotMatrix);    
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
    
    %Rotate around the y-axis so that X points forward and Z points out.
    %This requires a 180 degree rotation
    rotMatrixY = createRotationOy(deg2rad(180));
    
    %Rotate surfaces
    for pp = 1:length(scapulaV)
        scapulaV(pp,:) = transformPoint3d(scapulaV(pp,:),rotMatrixY);    
    end
    clear pp
    for pp = 1:length(humerusV)
        humerusV(pp,:) = transformPoint3d(humerusV(pp,:),rotMatrixY);    
    end
    clear pp
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),rotMatrixY);    
    end
    clear pp

    %Rotate landmarks
    for ff = 1:length(currLandmarks)
        landmarks.(currLandmarks{ff}) = transformPoint3d(landmarks.(currLandmarks{ff}),rotMatrixY);
    end
    clear ff

    %Rotate shapes
    shapes.headSphere.centre = transformPoint3d(shapes.headSphere.centre,rotMatrixY);

    %Rotate planes
    planes.glenoid.origin = transformPoint3d(planes.glenoid.origin,rotMatrixY);
    planes.glenoid.normal = transformPoint3d(planes.glenoid.normal,rotMatrixY);
    
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
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),transMatrix);    
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
        gpatch(scapulaF,scapulaV,'kw','none',0.8);
        gpatch(humerusF,humerusV,'kw','none',0.8);
        gpatch(headF,headV,'rw','none',1);
        %Landmarks
        for ff = 1:length(currLandmarks)
            plotV(landmarks.(currLandmarks{ff}),'g.','MarkerSize',15);
        end
        clear ff
        axisGeom;
        title('Coordinate System Aligned Surfaces');
    end
    
    %% Create object coordinate systems to align with one another
    
    %% Create scapula coordinate system
    
    %Update waitbar
    waitbar(3/9,wbar,'Creating object coordinate systems...');
    
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
        gpatch(scapulaF,scapulaV,'kw','none',0.5);
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
        %Title
        title('Scapula Coordinate System');
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
        gpatch(humerusF,humerusV,'kw','none',0.5);
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
        %Title
        title('Humerus Coordinate System')
    end

    %% Align humerus with scapula and translate to a 'neutral' position
    %  Here the neutral position corresponds to a distance of 10mm for the
    %  GHJC from the deep glenoid point. This assists with the basic force
    %  application later in FEA simulations.
    
    %Update waitbar
    waitbar(4/9,wbar,'Aligning coordinate systems...');
    
    %Rotate the humerus and associated landmarks to align with the scapula CS

    %Do X-axes alignment first

    %Create rotation matrix to align Xc vectors
    humerusRotX = createRotationVector3d(humerusCS.Xc(4:end),scapulaCS.Xc(4:end));

    %Rotate humerus points
    for pp = 1:length(humerusV)
        humerusV(pp,:) = transformPoint3d(humerusV(pp,:),humerusRotX);    
    end
    clear pp
    
    %Rotate head points
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),humerusRotX);    
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
    
    %Rotate head points
    for pp = 1:length(headV)
        headV(pp,:) = transformPoint3d(headV(pp,:),humerusRotZ);    
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
    desiredGHJCpoint = [0,0,10 + shapes.headSphere.radius];

    %Calculate required shift to get this based on the current GHJC position
    humerusTranslate = landmarks.GHJC - desiredGHJCpoint;

    %Shift humerus points
    for pp = 1:length(humerusV)
        humerusV(pp,:) = humerusV(pp,:) - humerusTranslate;
    end
    clear pp
    
    %Shift head points
    for pp = 1:length(headV)
        headV(pp,:) = headV(pp,:) - humerusTranslate;
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

    %Visualise overall aligned system
    if generatePlots
        cFigure; hold on
        %Bones
        gpatch(scapulaF,scapulaV,'kw','none',0.5);
        gpatch(humerusF,humerusV,'kw','none',0.5);
        gpatch(headF,headV,'rw','none',0.8);
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
        %Title
        title('Rotated and Aligned Surfaces');
    end
    
    %% Estimate GH rotation centre from regression...
    
    %PC - 
    %AI - 
    %AA - 
    
    xC = 18.9743 + (PCx * 0.2434) + (AIx * 0.2341) + (L_AI_AA * 0.1590) + (PCy * 0.0558)

    %% Extract glenoid section from whole scapula

    %Update waitbar
    waitbar(5/9,wbar,'Extracting glenoid...');
    
    %Slice the scapula 10mm back from the glenoid origin (+ve Z-direction)
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
        hp1 = gpatch(scapulaFc(~logicSide,:),scapulaVc,'bw','none',1);
        hp2 = gpatch(scapulaFc(logicSide,:),scapulaVc,'rw','none',1);
        legend([hp1 hp2],{'Surface above plane','Surface below plane'})
        axisGeom; axis manual; camlight headligth;
        colormap gjet;
        %Plot extracted surface and boundary
        subplot(1,2,2); hold on;
        gpatch(scapulaFc(logicSide,:),scapulaVc,'w','none',1);
        gpatch(scapulaFc(~logicSide,:),scapulaVc,'w','none',0.25);
        hp1=gpatch(scapulaEc,scapulaVc,'none','b',1,3);
        hp2=quiverVec(P,n,0.05,'k');
        legend([hp1 hp2],{'Intersection curve','Plane normal vector'})
        axisGeom; axis manual; camlight headligth;
    end

    %Extract the faces we want to keep
    [scapulaKeepF,scapulaKeepV] = patchCleanUnused(scapulaFc(logicSide,:),scapulaVc);

    %Use the grouping function to split the extra part of the scapula that
    %comes through with the cut away
    [~,indF] = groupVertices(scapulaKeepF,scapulaKeepV,1);

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
        title('Boundary of Cut Glenoid Surface')
    end

    %Create a surface that closes the back of the glenoid
    %Currently uses the default 0.5 point spacing of the scapula
    [backF,backV] = regionTriMesh2D({extractedGlenoidV(indBoundaryBack,[1 2])},0.5,0,0);
    backV(:,3) = mean(extractedGlenoidV(indBoundaryBack,3)); %Add/set z-level converting to 3D mesh

    %Visualise new meshes
    if generatePlots
        cFigure; hold on;
        gpatch(extractedGlenoidF,extractedGlenoidV,'bw','k');
        gpatch(backF,backV,'gw','k');
        plotV(extractedGlenoidV(indBoundaryBack,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        title('Filled Back of Cut Glenoid Surface');
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
        title('Joined Surfaces for Cut and Filled Glenoid');
    end

    %Check boundaries to make sure there are no holes. If there is throw an
    %error as the volumetric meshing won't work.
    if ~isempty(patchBoundary(glenoidF,glenoidV))
        error('Holes detected in glenoid. Stopping here as volumetric meshing won''t work');
    end
    
    %% Create a hemispherical mesh for the articulating humeral head surface
    %  NOTE: this creates a solid hexahedral hemisphere mesh, unlike the
    %  existing aspects which have, so far, only created surface meshes.
    
    %Update waitbar
    waitbar(6/9,wbar,'Creating hemispherical mesh...');
    
    %Identify the hemisphere radius required based on the articulating
    %surface extracted from the humeral head. This will identify the
    %boundary points of the articulating surface mesh and fit an idealised
    %circle to these 3D points. We'll take this as the radius for the
    %hemisphere, and align it so that it's centre sits on an edge point of
    %the approximate humeral head. We project a point from the centre of
    %the fitted 3D circle to the idealised humeral head sphere, and align
    %the edge of the hemisphere to this along the circles normal.
    
    %Create a boundary at the articulating surface mesh edges
    headBoundaryEb = patchBoundary(headF,headV); %Get boundary edges
    indBoundary = edgeListToCurve(headBoundaryEb); %Convert boundary edges to a curve list
    indBoundary = indBoundary(1:end-1); %Trim off last point since it is equal to first on a closed loop
    
    %Extract the points corresponding to the boundary
    headBoundaryPts = headV(indBoundary,:);
    
    %Fit a 3D circle to these boundary points
    [headBoundaryCircle,headBoundaryNormal] = fitCircle3d(headBoundaryPts);

    %Visualise the boundary on the cut surface and the fitted circle
    if generatePlots
        cFigure; hold on;
        gpatch(headF,headV,'kw','none',0.5);
        plotV(headV(indBoundary,:),'r-','LineWidth',2);
        camlight('headlight');
        axisGeom;
        drawCircle3d(headBoundaryCircle,'g','LineWidth',2);
        title('Boundary of Humeral Articulating Surface with Fitted Circle')
    end
    
    %Set hemisphere radius (4th value of 3d circle)
    hemiRadius = headBoundaryCircle(4);
    
    %Identify the intersecting point of a line along the 3d circles normal
    %to that of an idealised sphere fit to the humeral head articulating
    %circle.
    
    %Create the sphere (i.e. xc yc zc r)
    headSphere = [shapes.headSphere.centre(1),...
        shapes.headSphere.centre(2),...
        shapes.headSphere.centre(3),...
        shapes.headSphere.radius];
    
    %Create the line (i.e. x0 y0 z0 dx dy dz).
    %These points correspond to the fitted circle centre and it's normal
    headLine = [headBoundaryCircle(1), headBoundaryCircle(2), headBoundaryCircle(3),...
        headBoundaryNormal(1), headBoundaryNormal(2), headBoundaryNormal(3)];
    
    %Identify the intersection points of the line with the sphere
    headIntersectPts = intersectLineSphere(headLine,headSphere);
    
    %The point we're interested in will be the one with a Z-value closer to
    %zero, as this will be the intersection closer to the glenoid.
    if headIntersectPts(1,3) < headIntersectPts(2,3)
        keepIntersectPt = headIntersectPts(1,:);
    else
        keepIntersectPt = headIntersectPts(2,:);
    end
    
    %Visualise the circle, line and intersection along with the head surface
    if generatePlots
        cFigure; hold on;
        gpatch(headF,headV,'kw','none',0.5);
        drawCircle3d(headBoundaryCircle,'g','LineWidth',2);
        drawLine3d(headLine,'b')
        drawPoint3d(keepIntersectPt,'rx')
        camlight('headlight');
        axisGeom;        
        title('Intersection of Fitted Circle Normal to Idealised Humeral Head Sphere')
    end
    
    %Create control structure for creating hemisphere mesh
    %NOTE: this process creates a sphere with a core and mantel, but in
    %essence we just consider it one solid structure. We also use some
    %default parameters for mesh size here --- TODO: adapt with mesh
    %refinement?
    optionStruct.sphereRadius = hemiRadius;
    optionStruct.coreRadius = optionStruct.sphereRadius/5;
    optionStruct.numElementsMantel = 6;
    optionStruct.numElementsCore = optionStruct.numElementsMantel*2;
    optionStruct.outputStructType = 2;
    optionStruct.makeHollow = 0;
    optionStruct.cParSmooth.n = 25;

    %Create head mesh
    [headMeshOutput] = hexMeshHemiSphere(optionStruct);

    %Access the head mesh data
    hemiFb = headMeshOutput.facesBoundary;
    hemiCb = headMeshOutput.boundaryMarker;
    hemiV = headMeshOutput.nodes;
    hemiE = headMeshOutput.elements;
    
    %The hemispheres back surface firstly needs to be aligned to the same
    %plane as the fitted 3d circle.
    
    %Create a plane representing the 3d circle, using it's centre point and
    %normal direction
    headBoundaryPlane = createPlane(headBoundaryCircle(1:3),headBoundaryNormal);
    
    %The hemisphere sits on the global plane, so we can create the
    %transform between these two planes
    hemiTransform = createBasisTransform3d(headBoundaryPlane,'global');
    
    %Rotate hemisphere points
    for pp = 1:length(hemiV)
        hemiV(pp,:) = transformPoint3d(hemiV(pp,:),hemiTransform);    
    end
    clear pp
    
    %Calculate the difference between the hemisphere radius and the
    %distance from the 3d circle centre to the intersection point
    hemiSizeDiff = hemiRadius - (distancePoints3d(headBoundaryCircle(1:3),keepIntersectPt));
    
    %Create the translation matrix for moving hemisphere
    hemiTransMat = createTranslation3d((headBoundaryNormal*-1)*hemiSizeDiff);
    
    %Translate the hemisphere points
    for pp = 1:length(hemiV)
        hemiV(pp,:) = transformPoint3d(hemiV(pp,:),hemiTransMat);    
    end
    clear pp
        
    %Visualise aligned hemisphere
    if generatePlots
        cFigure; hold on
        gpatch(hemiFb,hemiV,'bw','none',0.8)
        gpatch(headF,headV,'kw','none',0.5);
        drawCircle3d(headBoundaryCircle,'g','LineWidth',2);
        camlight('headlight');        
        axisGeom;
        title('Hemisphere Aligned to Humeral Articulating Surface')
    end
    
    %% Create volumetric mesh of glenoid
    
    %Update waitbar
    waitbar(7/9,wbar,'Creating glenoid mesh...');
    
    %Remesh glenoid to start at a 1.0mm edge length value
    [glenoidF,glenoidV] = triRemeshLabel(glenoidF,glenoidV,1.0);
    
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
    glenoidVolFb = glenoidMesh.facesBoundary; %The boundary faces
    glenoidVolCb = glenoidMesh.boundaryMarker; %The boundary markers

    %Visualise glenoid mesh
    if generatePlots
        meshView(glenoidMesh);
        title('Tetrahedral Mesh: Glenoid');
    end
    
    %% Visualise entire system together
    
    if generatePlots
        cFigure; hold on;
        gpatch(scapulaF,scapulaV,'kw','none',0.8);
        gpatch(humerusF,humerusV,'kw','none',0.8);
        gpatch(hemiFb,hemiV,'bw','none',1.0)
        gpatch(glenoidF,glenoidV,'rw','none',1.0)
        camlight headlight;
        axisGeom;
        title('Total Meshed Glenohumeral Base System');
    end
    
    %% Package mesh data for output
    
    %Update waitbar
    waitbar(8/9,wbar,'Setting up outputs...');
    
    %Glenoid volumetric mesh
    glenoidOutput.glenoidVolFb = glenoidVolFb;
    glenoidOutput.glenoidVolV = glenoidVolV;
    glenoidOutput.glenoidVolE = glenoidVolE;
    glenoidOutput.glenoidV = glenoidV;
    glenoidOutput.glenoidF = glenoidF;
    
    %Hemisphere volumetric mesh
    headOutput.headVolFb = hemiFb;
    headOutput.headVolV = hemiV;
    headOutput.headVolE = hemiE;
    headOutput.headVolCb = hemiCb;
    
    %Scapula surface mesh
    scapulaOutput.scapulaF = scapulaF;
    scapulaOutput.scapulaV = scapulaV;
    
    %Humerus surface mesh
    humerusOutput.humerusF = humerusF;
    humerusOutput.humerusV = humerusV;
    
    %% Return to home directory
    cd(home_dir)
    
    %Update waitbar
    waitbar(9/9,wbar,'Import and creating surfaces done!');
    
    %Close waitbar
    close(wbar)

%% %%%%% ----- End of importAndCreateSurfaces.m ----- %%%%%
end