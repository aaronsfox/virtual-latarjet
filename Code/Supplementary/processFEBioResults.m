function [processedOutputs] = processFEBioResults(febioFileNamePart,glenoidMeshOutput,headMeshOutput,...
    landmarks,exportAnimation,generatePlots)

    %% This function serves to process in the results of an FEBio run and
    %  collate the relevant output metrics from the simulation.
    %
    %  Inputs:
    %
    %   febioFileNamePart   name the FEBio file outputs were stored under.
    %                       This must be the prefix to the FEBio file (i.e.
    %                       'febioFileNamePart.feb')
    %   glenoidMeshOutput   structure containing glenoid mesh output for
    %                       the FEBio simulation (i.e. not the original)
    %   headMeshOutput      structure containing head mesh output for
    %                       the FEBio simulation (i.e. not the original)
    %   landmarks           structure containing imported landmark points
    %                       on scapula/humerus
    %   exportAnimation     flag whether to export the animated gif of the
    %                       FEA simulation via GIBBON   
    %   generatePlots       flag whether to generate figures from the
    %                       processing throughout the function (default = false)
    %   
    %  Outputs
    %
    %   processedOutputs    structure containing the distance to
    %                       dislocation and stability ratio calculations
    %
    %  This function, like parts of the main code, uses elements of the GIBBON
    %  toolbox (https://www.gibboncode.org/) and the geom3d MATLAB package
    %  (https://au.mathworks.com/matlabcentral/fileexchange/24484-geom3d).
    
    warning off
    
    %%
    
    %%%%%%%% TODO: INSTEAD OF A SIM TYPE IT'S PROBABLY BETTER JUST TO
    %%%%%%%% INCORPORATE BOTH FEBIO RESULT SETS IN THIS FUNCTION

    %% Check inputs

    %%%%% TODO: add more input checks if necessary
    
    %Check for filename
    if nargin < 1
        error('An FEBio file name part is required. This is the prefix name for the .feb file')
    end
    
    %% Set-up
    
    %Set filename variable
    fileNameTranslate = [febioFileNamePart,'_translation'];
    fileNameForce = [febioFileNamePart,'_force'];
    
    %Import the output log file for displacements from translation sim
    [timeMatTranslate,nodeDispMatTranslate] = importFEBio_logfile([fileNameTranslate,'_dispHead_out.txt']); %Nodal displacements
    
    %Import the output log file for displacements from force sim
    [timeMatForce,nodeDispMatForce] = importFEBio_logfile([fileNameForce,'_dispHead_out.txt']); %Nodal displacements
    
    %% Calculate dislocation point data from the translation simulation
    
    %%%%% Need to just consider if this is correct - the current
    %%%%% dislocation specification is identifying the point where it gets
    %%%%% to the lip and the starts dropping down, but is still in contact
    %%%%% with the glenoid -- is this correct? It probably is a point of no
    %%%%% return in reality...
    
    %Find where anterior translation starts
    indStartTrans = find(diff(nodeDispMatTranslate(:,2)) > 0,1);
    
    %Find where anterior translation equates to radius of glenoid length as
    %a section to fit the spline curve to
    
    %Create a plane to project the superior and inferior glenoid points on
    %to so that they can be measured along the same axis
    xyPlane = [0 0 0 1 0 0 0 1 0];
    
    %Project the superior and inferior glenoid points onto this plane
    projSupGlenoid = projPointOnPlane(landmarks.SupGlenoid,xyPlane);
    projInfGlenoid = projPointOnPlane(landmarks.InfGlenoid,xyPlane);
    
    %Calculate the distance between the points for glenoid length
    glenoidLength = distancePoints3d(projSupGlenoid,projInfGlenoid);
    
    %Find point where translation is 0.75 of glenoid length
    indEndTrans = find(nodeDispMatTranslate(:,2) > (glenoidLength*0.75),1);
    
    %Smooth values
    windowSize = 15; %%%%% seems to result in appropriately smoothed curve
    smoothVals = smooth(nodeDispMatTranslate(:,2),nodeDispMatTranslate(:,4),windowSize);
    
    %Visualise translation vs. compression curve
    if generatePlots
        figure; hold on
        plot(nodeDispMatTranslate(indStartTrans:indEndTrans,2),...
            smoothVals(indStartTrans:indEndTrans,1),...
            '-b','LineWidth',2)
        plot(nodeDispMatTranslate(indStartTrans:indEndTrans,2),...
            nodeDispMatTranslate(indStartTrans:indEndTrans,4),'r.')
        xlabel('Anterior Translation');
        ylabel('Compressive Translation');
    end
    
    %Find peaks of smoothed curve
    [~,dislocDist] = findpeaks(smoothVals(indStartTrans:indEndTrans,1),...
        nodeDispMatTranslate(indStartTrans:indEndTrans,2));    
    
% % %     %Find the peak of normal points
% % %     [~,dislocDist] = findpeaks(nodeDispMat(indStartTrans:indEndTrans,4),...
% % %         nodeDispMat(indStartTrans:indEndTrans,2));
    
    %Check for singular peak
    if size(dislocDist,1) ~= 1
        error('More or less than one peak found in anterior-compression translation curve.')
    end
    
    %Find the nearest index to this for animation purposes
    dislocInd = find(nodeDispMatTranslate(1:end,2) > dislocDist,1);
    
    %% Calculate the force required to get to the dislocation distance
    
    %Identify the index where the dislocation distance is exceeded
	dislocForceInd = find(nodeDispMatForce(:,2) > dislocDist,1) - 1;
    
    %Get the time where dislocation occurred
    dislocTime = timeMatForce(dislocForceInd);

    %Identify the translation force being applied at the point of dislocation
    %This considers that at 2.1 seconds 10N of force was being applied, and
    %force increases at a rate of 10N per second in the simulation set-up
    forceStartTime = 2.1;
    forceStartMagnitude = 10;
    forcePerSecond = 10;

    %Identify difference between force start and dislocation times
    dislocTimeDiff = dislocTime - forceStartTime;

    %Calculate dislocation force based on time difference and force rate
    dislocForce = forceStartMagnitude + (forcePerSecond * dislocTimeDiff);

    %Calculate stability ratio, considering compressive force of 100N
    compressForce = 100;
    stabilityRatio = dislocForce / compressForce * 100;
    
    %% Animate the simulation from start up to dislocation (+1 sec)
    
    if exportAnimation
    
        %Set scaled RGB bone colour
        boneRGB = [251/255 244/255 231/255];
        
        %Get starting point displacement based on when translation starts
        DNstart = nodeDispMatTranslate(indStartTrans,2:4);
        
        %Get the displaced nodal coordinates, factoring in starting position
        V_defStart = headMeshOutput.headVolV + DNstart;
        
        %Create basic view and store graphics handle to initiate animation    
        hf = cFigure;
        gtitle([febioFileNamePart,' Animation']);
        hp1 = gpatch(headMeshOutput.headVolFb,V_defStart,boneRGB,'none');
        hp2 = gpatch(glenoidMeshOutput.glenoidF,glenoidMeshOutput.glenoidVolV,boneRGB,'none'); %Add graphics object to animate
        axisGeom; camlight headlight

        %Set animation time
        %Go to the dislocation point plus a second
        animEndInd = find(timeMatTranslate > (timeMatTranslate(dislocInd) + 1.0),1);
        animStruct.Time = timeMatTranslate(indStartTrans:animEndInd);

        %Set the axes to stay in their current position
        ax = gca;
        axis([ax.XLim(1) ax.XLim(2) ax.YLim(1) ax.YLim(2) ax.ZLim(1) ax.ZLim(2)]);

        %Loop through time steps to set animation
        for tt = indStartTrans:1:animEndInd

            %Get the current displacement value
            DN = nodeDispMatTranslate(tt,2:4);

            %Get the displaced nodal coordinates, factoring in starting position
            V_def = V_defStart + DN;

            %Set entries in animation structure
            animStruct.Handles{tt-(indStartTrans-1)} = [hp1]; %[hp1 hp1 hp2]; %Handles of objects to animate
            animStruct.Props{tt-(indStartTrans-1)} = {'Vertices'}; %{'Vertices','CData','Vertices'}; %Properties of objects to animate
            animStruct.Set{tt-(indStartTrans-1)} = {V_def};%{V_def,CF,V_def}; %Property values for to set in order to animate

        end
        clear tt

        %Initiate the animation
        anim8(hf,animStruct);

        %Export animation as GIF
        exportGifAnim8()

        %Close figure
        close all

        %Move and rename GIF
        cd('efw');
        f = dir('*.gif');
        copyfile(f.name,[febioFileNamePart,'.gif']);
        delete(f.name);
        movefile([febioFileNamePart,'.gif'],['..\',febioFileNamePart,'.gif'])
        cd('..');
        rmdir('efw');
    
    end
        
    %% Outputs
    
    %Output the dislocation distance stability ratio and
    processedOutputs.stabilityRatio = stabilityRatio;
    processedOutputs.distanceToDislocation = dislocDist;
        
%% %%%%% ----- End of importFEBioResults.m ----- %%%%%
end