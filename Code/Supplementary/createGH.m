function createGH()

%% This function ...


%%%%% Function looks like it works --- but still needs to be converted back
%%%%% to the original coordinate system (i.e. revert back, which could be
%%%%% problematic...)


% References
%
% Meskers et al. (1998). In vivo estimation of the glenohumeral joint
% rotation center from scapular bony landmarks by linear regression. J
% Biomech, 31: 93-96.


%inputs - landmarks


%% Create local scapula coordinate system at AC

% Note this scapula coordinate system is slightly different to the ISB
% recommended, and follows that of Meskers et al. (1998)

%Create line that extends from AC, along the same direction as TS to AC
%This represents the Z-axis of the scapula
%First get direction of TS to AC
TS_AC = createLine3d(landmarks.TS,landmarks.AC);
%Keep the direction but start the line from the AA origin
scapulaCS.Xc = createLine3d(landmarks.AC,TS_AC(4),TS_AC(5),TS_AC(6));

%Define the plane passing through AC, TS and AI
scapBodyPlane = createPlane(landmarks.TS,landmarks.AC,landmarks.AI);
%Get normal of the plane
scapBodyNormal = planeNormal(scapBodyPlane);
%Create line at AC origin with normal of this plane
scapulaCS.Zc = createLine3d(landmarks.AC,scapBodyNormal(1),scapBodyNormal(2),scapBodyNormal(3));

%Take cross product of other two axes
Ytemp = crossProduct3d(scapulaCS.Xc(4:end),scapulaCS.Zc(4:end))*-1;
scapulaCS.Yc = [landmarks.AC,Ytemp];

%% Orient the data to the scapula coordinate system

%Create a plane based on the AC origin and use the forward facing X-axis as
%the normal
planeAC = createPlane(landmarks.AC,scapulaCS.Xc(4:end));

%Create the transformation matrix to shift the landmark data to the scaupla
%coordinate system
yzPlane = createPlane([0 0 0],[0 1 0],[0 0 1]);
transformAC = createBasisTransform3d('global',planeAC);

%Transform the landmarks. In doing so create standalone variables to work with
AC = transformPoint3d(landmarks.AC,transformAC);
AA = transformPoint3d(landmarks.AA,transformAC);
AI = transformPoint3d(landmarks.AI,transformAC);
TS = transformPoint3d(landmarks.TS,transformAC);
PC = transformPoint3d(landmarks.PC,transformAC);


%% Reset the coordinate system with the new landmarks

%Create line that extends from AC, along the same direction as TS to AC
%This represents the Z-axis of the scapula
%First get direction of TS to AC
TS_AC = createLine3d(TS,AC);
%Keep the direction but start the line from the AA origin
Xc = createLine3d(AC,TS_AC(4),TS_AC(5),TS_AC(6));

%Define the plane passing through AC, TS and AI
scapBodyPlane = createPlane(TS,AC,AI);
%Get normal of the plane
scapBodyNormal = planeNormal(scapBodyPlane);
%Create line at AC origin with normal of this plane
Zc = createLine3d(AC,scapBodyNormal(1),scapBodyNormal(2),scapBodyNormal(3));

%Take cross product of other two axes
Ytemp = crossProduct3d(Xc(4:end),Zc(4:end))*-1;
Yc = [AC,Ytemp];

%% Realign the coordinate system so that the X-axis is the global X-axis

%Create rotation vector
rotVec = createRotationVector3d(Xc(4:end),[1 0 0]);

%Apply rotation vector to landmarks and coordinate systems
AC = transformPoint3d(AC,rotVec);
AA = transformPoint3d(AA,rotVec);
AI = transformPoint3d(AI,rotVec);
TS = transformPoint3d(TS,rotVec);
PC = transformPoint3d(PC,rotVec);
Xc(4:end) = transformPoint3d(Xc(4:end),rotVec);
Yc(4:end) = transformPoint3d(Yc(4:end),rotVec);
Zc(4:end) = transformPoint3d(Zc(4:end),rotVec);

%% Align the Z coordinate system to the global z-axis

%Create rotation vector
rotVec2 = createRotationVector3d(Zc(4:end),[0 0 1]);

%Apply rotation vector to landmarks and coordinate systems
AC = transformPoint3d(AC,rotVec2);
AA = transformPoint3d(AA,rotVec2);
AI = transformPoint3d(AI,rotVec2);
TS = transformPoint3d(TS,rotVec2);
PC = transformPoint3d(PC,rotVec2);
Xc(4:end) = transformPoint3d(Xc(4:end),rotVec2);
Yc(4:end) = transformPoint3d(Yc(4:end),rotVec2);
Zc(4:end) = transformPoint3d(Zc(4:end),rotVec2);

%% Calculate GH landmark according to Meskers et al. (1998)

%Determine length from AI to AA
lenAI_AA = distancePoints3d(AI,AA);

%Determine length from AC to AA
lenAC_AA = distancePoints3d(AC,AA);

%Determine length from TS to PC
lenTS_PC = distancePoints3d(TS,PC);

%Determine length from AC to PC
lenAC_PC = distancePoints3d(AC,PC);

%Calculate regression equations
GH(1) = 18.9743 + (PC(1) * 0.2434) + (AI(1) * 0.2341) + ...
    (lenAI_AA * 0.1590) + (PC(2) * 0.0558);
GH(2) = -3.8791 + (lenAC_AA * -0.3940) + ...
    (PC(2) * 0.1732) + (AI(1) * 0.1205) + (lenAC_PC * -0.1002);
GH(3) = 9.2629 + (PC(3) * 1.0255) + (landmarks.PC(2) * -0.2403) + (lenTS_PC * 0.1720);



%%

%Rotate head points
for pp = 1:length(scapulaV)
    V(pp,:) = transformPoint3d(V(pp,:),rotVec2);    
end
clear pp


%Visualise axes
csLabels = [{'Xc'}; {'Yc'}; {'Zc'}];
csColours = [{'r'}; {'g'}; {'b'}];
if generatePlots
    cFigure; hold on
    %Plot scapula
% % %     gpatch(scapulaF,V,'kw','none',0.5);
    axisGeom; camlight headlight;
    %Plot scapula origin point (use x-axes coordinates --- same across other axes)
    scatter3(Xc(1),Xc(2),Xc(3),1e3,'y.');
    %Plot axes
    for cc = 1:length(csLabels)
        %Get points and vectors
        if cc == 1
            x = Xc(1); y = Xc(2); z = Xc(3);
            u = Xc(4); v = Xc(5); w = Xc(6);
        elseif cc == 2
            x = Yc(1); y = Yc(2); z = Yc(3);
            u = Yc(4); v = Yc(5); w = Yc(6);            
        elseif cc == 3
            x = Zc(1); y = Zc(2); z = Zc(3);
            u = Zc(4); v = Zc(5); w = Zc(6);
        end
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


