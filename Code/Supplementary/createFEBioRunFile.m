function [glenoidMeshOutput,headMeshOutput] = createFEBioRunFile(glenoidMesh,headMesh,febioFebFileNamePart,...
    elevationAngle,rotationAngle,translateDirection,scapulaCS,humerusCS,landmarks,pointSpacing,generatePlots)

    %% This function serves to import in and create the base surface system for
    %  running the FEA analysis of the humeral head against the glenoid.
    %
    %  Inputs:
    %
    %   glenoidMesh         structure containing surface and volumetric aspects
    %                       of the current glenoid mesh
    %   headMesh            structure containing surface and volumetric aspects
    %                       of the current head mesh
    %   febioFileNamePart   name to store FEBio file outputs under
    %   elevationAngle      humeral elevation angle to run simulations at (in
    %                       degrees)
    %   rotationAngle       humeral axial rotation angle to run simulations at
    %                       (in degrees). Internal rotation is +ve while
    %                       external rotation is -ve
    %   translateDirection  a value containing the 'clock-face' direction
    %                       to translate the humeral head in. Currently
    %                       supports 3, 4 and 5
    %   scapulaCS           structure containing details of scapula coordinate
    %                       system so that the glenoid can be rotated
    %   humerusCS           structure containing details of humerus coordinate
    %                       system so that the glenoid can be rotated
    %   landmarks           structure containing imported landmark points
    %                       on scapula/humerus
    %   pointSpacing        point spacing of glenoid mesh
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
    if nargin < 9
        error('Need all function inputs except for pointSpacing and generatePlots')
    end
    if nargin < 11
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
% % %     if translateDirection == 3
% % %         translationVec = [-1 0 0];
% % %     elseif translateDirection == 4
% % %         translationVec = [-1 0 0];
% % %     elseif translateDirection == 5
% % %         translationVec = [-1 0 0];
% % %     else
% % %         error('Translation direction must be 3, 4 or 5 - representing this direction on a clock face.')
% % %     end
    
    %Create a string for the coordinate system labels
    csLabels = [{'Xc'}, {'Yc'}, {'Zc'}];
    csColours = [{'r'}; {'g'}; {'b'}];
    
    %Extract the nodes, elements and surfaces from the mesh structures
    glenoidVolV = glenoidMesh.glenoidVolV;
    headVolV = headMesh.headVolV;
    glenoidVolE = glenoidMesh.glenoidVolE;
    headVolE = headMesh.headVolE;
    glenoidF = glenoidMesh.glenoidF;
    headF = headMesh.headVolFb;
    
    %Calculate point spacing if not input
    if nargin < 9
        pointSpacing = mean(patchEdgeLengths(glenoidF,glenoidVolV));        
    end
    
    %Create a list humeral landmarks
    humerusLandmarks = [{'GHJC'},{'EL'},{'EM'},{'EJC'}];
    
    %Create a list of scapula landmarks
    scapulaLandmarks = [{'AA'},{'AI'},{'TS'},{'DeepGlenoid'},{'SGT'},{'IGT'},{'AntGlenoid'},{'PostGlenoid'},{'SupGlenoid'},{'InfGlenoid'}];
    
    %% Humeral elevation
    
    %Rotate the head and glenoid mesh according to the elevation angle
    %provided in the function. A 3:2 elevation ratio is adopted for moving
    %the glenoid relative to the humerus as per Walia et al. (2013, 2015).
    
    if elevationAngle > 0
    
        %Rotate relevant aspects of the humeral head

        %Create the transform about humeral coordinate system X axis
        humerusElvTrans = createRotation3dLineAngle(humerusCS.Xc,deg2rad(elevationAngle*-1)); %inverted to rotate the right way

        %Rotate nodes
        for pp = 1:length(headVolV)
            headVolV(pp,:) = transformPoint3d(headVolV(pp,:),humerusElvTrans);
        end
        clear pp

        %Rotate humeral coordinate system
        for cc = 1:length(csLabels)
            humerusCS.(csLabels{cc}) = transformLine3d(humerusCS.(csLabels{cc}),humerusElvTrans);
        end
        clear cc
        
        %Rotate humeral landmarks
        for mm = 1:length(humerusLandmarks)
            landmarks.(humerusLandmarks{mm}) = transformPoint3d(landmarks.(humerusLandmarks{mm}),humerusElvTrans);            
        end
        clear mm
           

        %Rotate the relevant aspects of the glenoid

        %Create the transform about scapula coordinate system X axis
        glenoidElvTrans = createRotation3dLineAngle(scapulaCS.Xc,deg2rad((elevationAngle)*(2/3)*-1)); %inverted to rotate the right way

        %Rotate nodes
        for pp = 1:length(glenoidVolV)
            glenoidVolV(pp,:) = transformPoint3d(glenoidVolV(pp,:),glenoidElvTrans); %needs to be inverted for appropriate rotation
        end
        clear pp
        
        %Rotate scapula coordinate system
        for cc = 1:length(csLabels)
            scapulaCS.(csLabels{cc}) = transformLine3d(scapulaCS.(csLabels{cc}),glenoidElvTrans);
        end
        clear cc
        
        %Rotate scapula landmarks
        for mm = 1:length(scapulaLandmarks)
            landmarks.(scapulaLandmarks{mm}) = transformPoint3d(landmarks.(scapulaLandmarks{mm}),glenoidElvTrans);            
        end
        clear mm
        
        %Rotate the glenoid origin point by the glenoid translation
        originPt = transformPoint3d(originPt,glenoidElvTrans);
        
        %Rotate the compression and translation vectors by the glenoid rotation
        compressionVec = transformVector3d(compressionVec,glenoidElvTrans);
% % %         translationVec = transformVector3d(translationVec,glenoidElvTrans);
        
    end
    
    %% Humeral axial rotation
    
    %Rotate the head by the axial rotation specified
    
    if rotationAngle ~= 0
    
        %Rotate relevant aspects of the humeral head

        %Create the transform about humeral coordinate system X axis
        humerusRotTrans = createRotation3dLineAngle(humerusCS.Yc,deg2rad(rotationAngle));

        %Rotate nodes
        for pp = 1:length(headVolV)
            headVolV(pp,:) = transformPoint3d(headVolV(pp,:),humerusRotTrans); %needs to be inverted for appropriate rotation
        end
        clear pp

        %Rotate humeral coordinate system
        for cc = 1:length(csLabels)
            humerusCS.(csLabels{cc}) = transformLine3d(humerusCS.(csLabels{cc}),humerusRotTrans);
        end
        clear cc
        
        %Rotate humeral landmarks
        for mm = 1:length(humerusLandmarks)
            landmarks.(humerusLandmarks{mm}) = transformPoint3d(landmarks.(humerusLandmarks{mm}),humerusRotTrans);            
        end
        clear mm
        
    end    
    
    %% Position humeral head against glenoid
    
    % In the starting position the humeral head centre was positioned 10mm 
    % plus its radius from the deep glenoid origin. With the glenoid likely
    % being rotated, the origin point has changed, and hence the humeral
    % head will be a different distance from this point now. Given we have
    % the updated origin point and the humeral head has simply rotated
    % about it's central point, we have both values to calculate the
    % updated distance. The humeral head was aligned with the original
    % Z-axis or compression vector, so given we have the updated
    % compression vector we also have the direction the humeral head sits
    % from the deep glenoid point.
    
    %Create the translation line from the humeral head to the current deep
    %glenoid origin point
    headToGlenoid = createLine3d(headPt,originPt);
    
    %Calculate the original and new distance from the head to glenoid
    origHeadDist = distancePoints3d(headPt,[0 0 0]);
    newHeadDist = distancePoints3d(headPt,originPt);
    
    %Calculate difference between original and new head distance
    diffHeadDist = newHeadDist - origHeadDist;
    
    %The original distance the head needs to move to be the desired 0.5mm
    %is 9.5mm. Set this as a variable
    origShift = 9.5;
    
    %The summation between the original shift and the difference between
    %the original and new head position is how far the head now needs to
    %move.
    newShift = diffHeadDist + origShift;
    
    %We need the relative magnitude of the new shift against the current
    %distance of the head from the glenoid point as to multiply the line we
    %created earlier by this relative factor in order to generate the
    %appropriate translation matrix
    relHeadShift = newShift / newHeadDist;
    
    %Create translation matrix to shift the humeral head along the head to
    %glenoid line, taking into account the relative difference in the new
    %shift distance by the current head to glenoid distance
    headTransMat = createTranslation3d(headToGlenoid(4:end)*relHeadShift);
    
    %Translate the head nodes
    for pp = 1:length(headVolV)
        headVolV(pp,:) = transformPoint3d(headVolV(pp,:),headTransMat);
    end
    clear pp
    
    %Translate humeral landmarks
    for mm = 1:length(humerusLandmarks)
        landmarks.(humerusLandmarks{mm}) = transformPoint3d(landmarks.(humerusLandmarks{mm}),headTransMat);            
    end
    clear mm
    
    %% Realign the glenoid and head to the world XY plane
    
    % To make things easier for the FE simulations with respect to the
    % application of forces and boundary conditions - this step realigns 
    % the bodies so that the glenoid origin point returns to the world
    % coordinate system origin, and the glenoid plane is coincident with
    % the world XY plane. This therefore means that the compression and
    % translation vectors can be applied along the Z and X axes,
    % respectively.
    
    %Create a plane at the new origin and glenoid Z-axis (i.e. compression)
    newGlenoidPlane = createPlane(originPt,compressionVec);
    
    %Create the world XY plane
    xyPlane = createPlane([0,0,0],[1,0,0],[0,1,0]);
    
    %Create the basis transform between the new glenoid plane and world XY plane
    worldTransform = createBasisTransform3d(xyPlane,newGlenoidPlane);
    
    %Transform surfaces
    for pp = 1:length(glenoidVolV)
        glenoidVolV(pp,:) = transformPoint3d(glenoidVolV(pp,:),worldTransform);    
    end
    clear pp
    for pp = 1:length(headVolV)
        headVolV(pp,:) = transformPoint3d(headVolV(pp,:),worldTransform);    
    end
    clear pp
    
    %Transform landmarks
    for mm = 1:length(scapulaLandmarks)
        landmarks.(scapulaLandmarks{mm}) = transformPoint3d(landmarks.(scapulaLandmarks{mm}),worldTransform);            
    end
    clear mm
    
    for mm = 1:length(humerusLandmarks)
        landmarks.(humerusLandmarks{mm}) = transformPoint3d(landmarks.(humerusLandmarks{mm}),worldTransform);            
    end
    clear mm
    
    %Do the rotation that aligns the Y-axis with the SGT to IGT plane.
    %Identify the angle need to rotate around the Z-axis to align the SGT
    %and IGT points.

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
    rotMatrix = createRotationOz(rotAng); %-ve for anti-clockwise

    %Rotate surfaces
    for pp = 1:length(glenoidVolV)
        glenoidVolV(pp,:) = transformPoint3d(glenoidVolV(pp,:),rotMatrix);    
    end
    clear pp
    for pp = 1:length(headVolV)
        headVolV(pp,:) = transformPoint3d(headVolV(pp,:),rotMatrix);    
    end
    clear pp

    %Rotate landmarks
    for mm = 1:length(scapulaLandmarks)
        landmarks.(scapulaLandmarks{mm}) = transformPoint3d(landmarks.(scapulaLandmarks{mm}),rotMatrix);            
    end
    clear mm    
    for mm = 1:length(humerusLandmarks)
        landmarks.(humerusLandmarks{mm}) = transformPoint3d(landmarks.(humerusLandmarks{mm}),rotMatrix);            
    end
    clear mm
    
    %%%%% TODO: check the robustness of these rotations with different
    %%%%% glenoid/head configurations --- theoretically it should be OK
    %%%%% here as the glenoid and humerus have already been configured to
    %%%%% be consistent, unlike at the import and create surfaces stage.
    
    %The glenoid is now aligned that the 3 o clock direction represents the
    %positive X-axis. If this is the translation direction input, no
    %further translation about the Z-axis is required. However, if we want
    %to do the 4 or 5 o clock translation, the objects need to be rotated
    %by a specific amount so that the x-axis lines up to these directions.
    %Given a clock face has 12 intervals, the amount of rotation for one
    %'tick' is 360/12 degrees --- we therefore use this information to
    %apply the appropriate amount of rotation around the Z-axis.
    
    %Check for error
    if ~ismember(translateDirection,[3 4 5])
        error('Translation direction must be 3, 4 or 5 - representing this direction on a clock face.')
    end
    
    %Check for the amount of rotation needed
    if translateDirection == 3
        
        %No rotation needs to be done
        
    else
    
        if translateDirection == 4
        
            %Specify the rotation magnitude (in degrees)
            rotMagDeg = 360/12 * 1;

            %Create the rotation transform about the pure Z-axis
            zRot = createRotationOz(deg2rad(rotMagDeg));
        
        elseif translateDirection == 5
            
            %Specify the rotation magnitude (in degrees)
            rotMagDeg = 360/12 * 2;

            %Create the rotation transform about the pure Z-axis
            zRot = createRotationOz(deg2rad(rotMagDeg));
            
        end
        
        %Apply the rotation to the surfaces
        for pp = 1:length(glenoidVolV)
            glenoidVolV(pp,:) = transformPoint3d(glenoidVolV(pp,:),zRot);    
        end
        clear pp
        for pp = 1:length(headVolV)
            headVolV(pp,:) = transformPoint3d(headVolV(pp,:),zRot);    
        end
        clear pp

    end
        
    
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

    %%%%% TODO: hemisphere does seem kind of 'un-centred' but is based on
    %%%%% alignment of humerus and scapula axes????? May be best to check
    %%%%% this on another participant at some stage
    
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
    
    %Compression load curve
    febio_spec.LoadData.loadcurve{1}.ATTR.id = 1;
    febio_spec.LoadData.loadcurve{1}.ATTR.type = 'linear';
    febio_spec.LoadData.loadcurve{1}.point.VAL = [0 0; 2 1];

    %Translation load curve
    febio_spec.LoadData.loadcurve{2}.ATTR.id = 2;
    febio_spec.LoadData.loadcurve{2}.ATTR.type = 'linear';
    febio_spec.LoadData.loadcurve{2}.point.VAL = [2 0; 2.1 0.1; 11.1 1];

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
    
    %%%%% TODO: is this step size // number of steps going to take too
    %%%%% long?
    
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

    %% Output section 
    
    %Defining file names
    febioFebFileName = [febioFebFileNamePart,'.feb']; %FEB file name
    febioLogFileName = [febioFebFileNamePart,'.txt']; %FEBio log file name
    febioLogFileName_dispHead = [febioFebFileNamePart,'_dispHead_out.txt']; %Log file name for exporting displacement
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
    headMeshOutput.headVolFb = headF;
    
%% %%%%% ----- End of createFEBioRunFile.m ----- %%%%%
end