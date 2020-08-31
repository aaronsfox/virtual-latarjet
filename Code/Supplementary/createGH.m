function GH = createGH(glenoidSurfacePoints,shapes,generatePlots)

%% This function calculates and creates the glenohumeral centre landmark as
%  per the fitting procedures outlined in Meskers et al. (1998) for use in
%  creating the glenhumeral joint coordinate system
%
%  Inputs:
%
%   glenoidSurfaceV     vertice points from the glenoid surface STL object
%   shapes              structure containing imported shapes 
%   generatePlots       flag whether to generate figures from the
%                       processing throughout the function (default = false)
%
%  Outputs
%
%   GH                  3D coordinates of GH point
%
%  This function, like parts of the main code, uses elements of the GIBBON
%  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
%  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).
%
% References
%
% Meskers et al. (1998). In vivo estimation of the glenohumeral joint
% rotation center from scapular bony landmarks by linear regression. J
% Biomech, 31: 93-96.

    %% Check inputs
    if nargin < 2
        error('At least the glenoid surface points and shapes structure are required as inputs');
    end
    if nargin < 3
        generatePlots = false;
    end

    %% Fit a sphere to the glenoid surface

    %Create sphere fitting least squares function as per Meskers et al. (1998)
    f = @(M,x,y,z,r)sum((sqrt(((x - M(1)).^2) + ((y - M(2)).^2) + ((z - M(3)).^2)) - r).^2);

    %Set the known parameters
    %Points on glenoid surface
    x = glenoidSurfacePoints(:,1); y = glenoidSurfacePoints(:,2); z = glenoidSurfacePoints(:,3);
    %Radius of the humeral head sphere
    r = shapes.headSphere.radius;

    %Create standalone function to solve for sphere centre
    fun = @(M)f(M,x,y,z,r);

    %Minimise the function
    %Start at the humeral head centre
    x0 = shapes.headSphere.centre;
    GH = fminsearch(fun,x0);

    %Visualise
    if generatePlots
       cFigure; hold on;
       gpatch(surfaceF,surfaceV,'bw','none',0.7);
       axisGeom; camlight headlight;
       drawSphere(GH,r)
    end

%% %%%%% ----- End of createGH.m ----- %%%%%

end

%% Below code is initial attempt at Meskers et al. regression
%  It is complete except for the fact of converting the GH back to the
%  original coordinate system for export

% % % %% Create local scapula coordinate system at AC
% % % 
% % % % Note this scapula coordinate system is slightly different to the ISB
% % % % recommended, and follows that of Meskers et al. (1998)
% % % 
% % % %Create line that extends from AC, along the same direction as TS to AC
% % % %This represents the Z-axis of the scapula
% % % %First get direction of TS to AC
% % % TS_AC = createLine3d(landmarks.TS,landmarks.AC);
% % % %Keep the direction but start the line from the AA origin
% % % scapulaCS.Xc = createLine3d(landmarks.AC,TS_AC(4),TS_AC(5),TS_AC(6));
% % % 
% % % %Define the plane passing through AC, TS and AI
% % % scapBodyPlane = createPlane(landmarks.TS,landmarks.AC,landmarks.AI);
% % % %Get normal of the plane
% % % scapBodyNormal = planeNormal(scapBodyPlane);
% % % %Create line at AC origin with normal of this plane
% % % scapulaCS.Zc = createLine3d(landmarks.AC,scapBodyNormal(1),scapBodyNormal(2),scapBodyNormal(3));
% % % 
% % % %Take cross product of other two axes
% % % Ytemp = crossProduct3d(scapulaCS.Xc(4:end),scapulaCS.Zc(4:end))*-1;
% % % scapulaCS.Yc = [landmarks.AC,Ytemp];
% % % 
% % % %% Orient the data to the scapula coordinate system
% % % 
% % % %Create a plane based on the AC origin and use the forward facing X-axis as
% % % %the normal
% % % planeAC = createPlane(landmarks.AC,scapulaCS.Xc(4:end));
% % % 
% % % %Create the transformation matrix to shift the landmark data to the scaupla
% % % %coordinate system
% % % yzPlane = createPlane([0 0 0],[0 1 0],[0 0 1]);
% % % transformAC = createBasisTransform3d('global',planeAC);
% % % 
% % % %Transform the landmarks. In doing so create standalone variables to work with
% % % AC = transformPoint3d(landmarks.AC,transformAC);
% % % AA = transformPoint3d(landmarks.AA,transformAC);
% % % AI = transformPoint3d(landmarks.AI,transformAC);
% % % TS = transformPoint3d(landmarks.TS,transformAC);
% % % PC = transformPoint3d(landmarks.PC,transformAC);
% % % 
% % % %% Reset the coordinate system with the new landmarks
% % % 
% % % %Create line that extends from AC, along the same direction as TS to AC
% % % %This represents the Z-axis of the scapula
% % % %First get direction of TS to AC
% % % TS_AC = createLine3d(TS,AC);
% % % %Keep the direction but start the line from the AA origin
% % % Xc = createLine3d(AC,TS_AC(4),TS_AC(5),TS_AC(6));
% % % 
% % % %Define the plane passing through AC, TS and AI
% % % scapBodyPlane = createPlane(TS,AC,AI);
% % % %Get normal of the plane
% % % scapBodyNormal = planeNormal(scapBodyPlane);
% % % %Create line at AC origin with normal of this plane
% % % Zc = createLine3d(AC,scapBodyNormal(1),scapBodyNormal(2),scapBodyNormal(3));
% % % 
% % % %Take cross product of other two axes
% % % Ytemp = crossProduct3d(Xc(4:end),Zc(4:end))*-1;
% % % Yc = [AC,Ytemp];
% % % 
% % % %% Realign the coordinate system so that the X-axis is the global X-axis
% % % 
% % % %Create rotation vector
% % % rotVec = createRotationVector3d(Xc(4:end),[1 0 0]);
% % % 
% % % %Apply rotation vector to landmarks and coordinate systems
% % % AC = transformPoint3d(AC,rotVec);
% % % AA = transformPoint3d(AA,rotVec);
% % % AI = transformPoint3d(AI,rotVec);
% % % TS = transformPoint3d(TS,rotVec);
% % % PC = transformPoint3d(PC,rotVec);
% % % Xc(4:end) = transformPoint3d(Xc(4:end),rotVec);
% % % Yc(4:end) = transformPoint3d(Yc(4:end),rotVec);
% % % Zc(4:end) = transformPoint3d(Zc(4:end),rotVec);
% % % 
% % % %% Align the Z coordinate system to the global z-axis
% % % 
% % % %Create rotation vector
% % % rotVec2 = createRotationVector3d(Zc(4:end),[0 0 1]);
% % % 
% % % %Apply rotation vector to landmarks and coordinate systems
% % % AC = transformPoint3d(AC,rotVec2);
% % % AA = transformPoint3d(AA,rotVec2);
% % % AI = transformPoint3d(AI,rotVec2);
% % % TS = transformPoint3d(TS,rotVec2);
% % % PC = transformPoint3d(PC,rotVec2);
% % % Xc(4:end) = transformPoint3d(Xc(4:end),rotVec2);
% % % Yc(4:end) = transformPoint3d(Yc(4:end),rotVec2);
% % % Zc(4:end) = transformPoint3d(Zc(4:end),rotVec2);
% % % 
% % % %% Rotate 90 degrees about the X-axis to match Meskers et al. (1998) coordinate system
% % % 
% % % %Create rotation matrix
% % % xRot = createRotationOx(deg2rad(-90));
% % % 
% % % %Apply rotation vector to landmarks and coordinate systems
% % % AC = transformPoint3d(AC,xRot);
% % % AA = transformPoint3d(AA,xRot);
% % % AI = transformPoint3d(AI,xRot);
% % % TS = transformPoint3d(TS,xRot);
% % % PC = transformPoint3d(PC,xRot);
% % % 
% % % %% Calculate GH landmark according to Meskers et al. (1998)
% % % 
% % % % Note that the Y-axis is inverted in our example (i.e. +ve is down) vs.
% % % % that of Meskers et al. There is therefore a need to invert Y-axes values
% % % 
% % % %Determine length from AI to AA
% % % lenAI_AA = distancePoints3d(AI,AA);
% % % 
% % % %Determine length from AC to AA
% % % lenAC_AA = distancePoints3d(AC,AA);
% % % 
% % % %Determine length from TS to PC
% % % lenTS_PC = distancePoints3d(TS,PC);
% % % 
% % % %Determine length from AC to PC
% % % lenAC_PC = distancePoints3d(AC,PC);
% % % 
% % % %Calculate regression equations
% % % GH(1) = 18.9743 + (PC(1) * 0.2434) + (AI(1) * 0.2341) + ...
% % %     (lenAI_AA * 0.1590) + (PC(2) * 0.0558);
% % % GH(2) = -3.8791 + (lenAC_AA * -0.3940) + ...
% % %     (PC(2) * 0.1732) + (AI(1) * 0.1205) + (lenAC_PC * -0.1002);
% % % GH(3) = 9.2629 + (PC(3) * 1.0255) + (PC(2) * -0.2403) + (lenTS_PC * 0.1720);