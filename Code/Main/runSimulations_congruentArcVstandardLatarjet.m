function runSimulations_congruentArcVstandardLatarjet()

%%%% TODO: function inputs like participant ID

%%%%% enter notes...references to relevant papers where necessary

%%%%% following a lot of the sphere sliding demo in GIBBON (still?)
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0007_sphere_sliding.html
%%%%% a fair bit for slicing and fixing up the hole also comes from the foot insole demo
    %%%%% https://www.gibboncode.org/html/DEMO_febio_0050_foot_insole_01.html

%%%%% TODO: add generatePlots logical input

%%%%% TODO: add notes on starting from own directory
    
    generatePlots = false;
    participantID = 'id12984';

    %%%%% TODO: mesh convergence...

    %%%%% TODO: add notes on GIBBON toolbox and geom3d being needed
    
    % References
    %
    % Itoi et al. (2000). Effect of a glenoid defect on anteroinferior
    % stability of the shoulder after Bankart repair. J Bone Joint Surg,
    % 82-A: 35-46.
    %
    % Klemt et al. (2019). The critical size of a defect in the glenoid
    % causing anterior instability of the shoulder after a Bankart repair,
    % under physiological joint loading. Bone Joint J, 101-B: 68-74.
    %
    % Walia et al. (2013). Theoretical model of the effect of combined
    % glenohumeral bone defects on anterior shoulder instability: A finite
    % element approach. J Orthop Res, 31: 601-607.
    %
    % Walia et al. (2015). Influence of combined Hill-Sachs and bony
    % Bankart defects on range of motion in anterior instability of the
    % shoulder in a finite element model. Arthroscopy, 31: 2119-2127.
    
    %% Current test code
    %  Building up piece by piece to get to final simulation

    %% Set-up

    %Turn warnings off as these can get annoying
    warning off

    %Add supplementary code folder to path
    addpath(genpath('..\Supplementary'));

    %% Load and create relevant surface and volumetric meshes

    % This step takes in relevant surface meshes exported from the 3matic
    % processing to create a base system to run FEA simulations and adapt with
    % bone defects to run subsequent simulations. As part of this the glenoid
    % section is extracted from the scapula, and the surfaces are rotated and
    % aligned so that the deep glenoid point is the origin and the Z-axis runs
    % away from the face (-ive direction).
    %
    % Relevant files from the participants segmentation directory are used in
    % the function that creates the base FEA system, so this directory is
    % passed as an unput to the function.

    %Navigate to the segmentation directory
    cd('..\..\Segmentation\');
    segmentationDir = [pwd,'\'];
    
    %Set the current participant segmentation directory
    currSegDir = [segmentationDir,participantID,'\'];
    
    %Run import and create surfaces function
    [scapulaMesh,humerusMesh,glenoidMesh,headMesh,...
        scapulaCS,humerusCS,landmarks,shapes] = ...
        importAndCreateSurfaces(currSegDir,generatePlots);
    
    %% Mesh refinement...
    
    %%%%% TODO: include an external function that does mesh refinement to
    %%%%% determine the mesh size needed for surfaces
    
    %% Create anterior glenoid bone defects
    
    %%%%% DEPRECATED: using only anteroinferior defects
    
    % This step takes in the extracted base glenoid surface and creates
    % anterior Bankart bone defects of varying sizes. The process included
    % in the function follows that of Klemt et al. (2019), whereby the
    % defects are created by performing a simulated osteotomy at various
    % distances relative to the lengh of the glenoid. The length of the
    % glenoid is determined based on landmarks added to the scapula surface
    % during segmentation.
    [anteriorDefectsMesh] = createAnteriorBankartDefects(glenoidMesh,landmarks,generatePlots);
    
    %%%%%% TODO: remeshing of glenoids could be necessary if we want to fix
    %%%%%% up the flat surfaces, but it might not matter?
    
    %Generate list of anterior defect meshes
    anteriorDefectsList = fieldnames(anteriorDefectsMesh);

    %Visualise anterior defects on a subplot
    if generatePlots
        cFigure;
        for dd = 1:length(anteriorDefectsList)
            %Plot current glenoid
            subplot(2,round(length(anteriorDefectsList)/2),dd)
            gpatch(anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidF,...
                anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidV,...
                'kw','none')
            axisGeom; camlight headlight
            %Get percentage number for title
            splitStr = strsplit((anteriorDefectsList{dd}),'_');
            per = sscanf(splitStr{2},'%f');
            title([num2str(per),'% Anterior Defect'])
            %Cleanup
            clear splitStr per
        end
        clear dd
    end
    
    %Create volumetric meshes of anterior defect glenoids
    
    %Create generic aspects of tetgen input structure
    tetGenInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
    tetGenInputStruct.holePoints = []; %Interior points for holes
    
    %Setup figure to visualise meshes
    if generatePlots
        hFig = cFigure;
    end
    
    %Start the waitbar
    wbar = waitbar(0/(length(anteriorDefectsList)),'Creating anterior defect volumetric meshes...');
    
    %Loop through anterior defects list
    for dd = 1:length(anteriorDefectsList)
        
        %Update waitbar
        splitStr = strsplit((anteriorDefectsList{dd}),'_');
        per = sscanf(splitStr{2},'%f');
        waitbar(dd/(length(anteriorDefectsList)),wbar,...
            ['Creating ',num2str(per),'% anterior defect volumetric mesh...']);
        
        %Mesh specific tetgen inputs
        tetGenInputStruct.Faces = anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidF; %Boundary faces
        tetGenInputStruct.Nodes = anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidV; %Nodes of boundary
        tetGenInputStruct.regionPoints= getInnerPoint(anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidF,...
            anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidV); %Interior points for regions
        tetGenInputStruct.regionA = tetVolMeanEst(anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidF,...
            anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidV); %Desired tetrahedral volume for each region
        
        %Mesh model
        [tetGenOutput] = runTetGen(tetGenInputStruct); %Run tetGen
        
        %Store tetgen outputs
        anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidVolE = tetGenOutput.elements; %The elements
        anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidVolV = tetGenOutput.nodes; %The vertices or nodes
        anteriorDefectsMesh.(anteriorDefectsList{dd}).glenoidVolFb = tetGenOutput.facesBoundary; %The boundary faces
        
        %Visualise current mesh on subplot figure
        if generatePlots
            %Set subplot
            hs = subplot(2,round(length(anteriorDefectsList)/2),dd);
            %Set title
            title(['Cut View of ',num2str(per),'% Anterior Defect Mesh'],...
                'FontSize',12);
            %Plot mesh
            optionStruct.hFig = [hFig hs];
            meshView(tetGenOutput,optionStruct);
        end
        
    end
    clear dd
    
    %Close wbar
    close(wbar);
    
    %%%%% TODO: internal anterior defect meshes look a little iffy...
    
    %% Create anteroinferior bone defects
    
% % %     %%%%% NOTE: code can create multiple anteroinferior defect sizes, but
% % %     %%%%% given simulation time with each of these a decision to focus on
% % %     %%%%% one, the 37.5% size was made    
% % %     
% % %     % This step takes in the extracted base glenoid surface and creates
% % %     % anterior Bankart bone defects of varying sizes. The process included
% % %     % in the function follows that of Walia et al. (2013, 2015) which is an
% % %     % adaptation of Itoi et al. (2000), whereby the defects are created by
% % %     % performing a simulated osteotomy at warious distances relative to the
% % %     % length of the glenoid on a 45 degree orientation. The length of the
% % %     % glenoid is determined based on landmarks added to the scapula surface
% % %     % during segmentation.
% % %     [anteroInferiorDefectsMesh] = createAnteroInferiorBankartDefects(glenoidMesh,landmarks,generatePlots);
% % %     
% % %     %%%%%% TODO: remeshing of glenoids could be necessary if we want to fix
% % %     %%%%%% up the flat surfaces, but it might not matter?
% % %     
% % %     %Generate list of anterior defect meshes
% % %     anteroInferiorDefectsList = fieldnames(anteroInferiorDefectsMesh);
% % % 
% % %     %Visualise anterior defects on a subplot
% % %     if generatePlots
% % %         if length(anteroInferiorDefectsList) > 1
% % %             cFigure;
% % %             for dd = 1:length(anteroInferiorDefectsList)
% % %                 %Plot current glenoid
% % %                 subplot(2,round(length(anteroInferiorDefectsList)/2),dd)
% % %                 gpatch(anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidF,...
% % %                     anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidV,...
% % %                     'kw','none')
% % %                 axisGeom; camlight headlight
% % %                 %Get percentage number for title
% % %                 splitStr = strsplit((anteroInferiorDefectsList{dd}),'_');
% % %                 per = sscanf(splitStr{2},'%f');
% % %                 title([num2str(per),'% Glenoid Radius Anteroinferior Defect'])
% % %                 %Cleanup
% % %                 clear splitStr per
% % %             end
% % %             clear dd
% % %         else
% % %             %Plot current glenoid
% % %             cFigure;
% % %             gpatch(anteroInferiorDefectsMesh.(anteroInferiorDefectsList{1}).glenoidF,...
% % %                 anteroInferiorDefectsMesh.(anteroInferiorDefectsList{1}).glenoidV,...
% % %                 'kw','none')
% % %             axisGeom; camlight headlight
% % %             %Get percentage number for title
% % %             splitStr = strsplit((anteroInferiorDefectsList{1}),'_');
% % %             per = sscanf(splitStr{2},'%f');
% % %             title([num2str(per),'% Glenoid Radius Anteroinferior Defect'])
% % %             %Cleanup
% % %             clear splitStr per
% % %         end
% % %     end
% % %     
% % %     %Create volumetric meshes of anterior defect glenoids
% % %     
% % %     %Create generic aspects of tetgen input structure
% % %     tetGenInputStruct.stringOpt = '-pq1.2AaY'; %Tetgen options
% % %     tetGenInputStruct.holePoints = []; %Interior points for holes
% % %     
% % %     %Setup figure to visualise meshes
% % %     if generatePlots
% % %         hFig = cFigure;
% % %     end
% % %     
% % %     %Start the waitbar
% % %     wbar = waitbar(0/(length(anteroInferiorDefectsList)),'Creating anteroinferior defect volumetric meshes...');    
% % %     
% % %     %Loop through anterior defects list
% % %     for dd = 1:size(anteroInferiorDefectsList,2)
% % %         
% % %         %Update waitbar
% % %         splitStr = strsplit((anteroInferiorDefectsList{dd}),'_');
% % %         per = sscanf(splitStr{2},'%f');
% % %         waitbar(dd/(length(anteroInferiorDefectsList)),wbar,...
% % %             ['Creating ',num2str(per),'% anteroinferior defect volumetric mesh...']);
% % %         
% % %         %Mesh specific tetgen inputs
% % %         tetGenInputStruct.Faces = anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidF; %Boundary faces
% % %         tetGenInputStruct.Nodes = anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidV; %Nodes of boundary
% % %         tetGenInputStruct.regionPoints= getInnerPoint(anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidF,...
% % %             anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidV); %Interior points for regions
% % %         tetGenInputStruct.regionA = tetVolMeanEst(anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidF,...
% % %             anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidV); %Desired tetrahedral volume for each region
% % %         
% % %         %Mesh model
% % %         [tetGenOutput] = runTetGen(tetGenInputStruct); %Run tetGen
% % %         
% % %         %Store tetgen outputs
% % %         anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidVolE = tetGenOutput.elements; %The elements
% % %         anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidVolV = tetGenOutput.nodes; %The vertices or nodes
% % %         anteroInferiorDefectsMesh.(anteroInferiorDefectsList{dd}).glenoidVolFb = tetGenOutput.facesBoundary; %The boundary faces
% % %         
% % %         %Visualise current mesh on subplot figure
% % %         if generatePlots
% % %             if length(anteroInferiorDefectsList) > 1
% % %                 %Set subplot
% % %                 hs = subplot(2,round(length(anteroInferiorDefectsList)/2),dd);
% % %                 %Set title
% % %                 title(['Cut View of ',num2str(per),'% Glenoid Length Anteroinferior Defect Mesh'],...
% % %                     'FontSize',11);
% % %                 %Plot mesh
% % %                 optionStruct.hFig = [hFig hs];
% % %                 meshView(tetGenOutput,optionStruct);
% % %             else
% % %                 %Set title
% % %                 title(['Cut View of ',num2str(per),'% Glenoid Length Anteroinferior Defect Mesh'],...
% % %                     'FontSize',11);
% % %                 %Plot mesh
% % %                 optionStruct.hFig = [hFig hs];
% % %                 meshView(tetGenOutput,optionStruct);
% % %             end
% % %         end
% % %         
% % %     end
% % %     clear dd
% % %     
% % %     %Close wbar
% % %     close(wbar);
    
    %% Run FEBio simulations
    
    %%%%% TODO: consider Hill-Sachs defect incorporation here too?
    
    %%%%% TODO: need to simplify mesh around edges, the computation time
    %%%%% for the smaller, gritty edges is quite large
    
    %Navigate to the FEA directory
    cd('..\FEA');
    
    %Create directory for current participant
    mkdir(participantID); cd(participantID);
    
    %Make study directory
    mkdir('congruentVstandard'); cd('congruentVstandard');
    
    %Create a list of the meshes to loop through. This includes the
    %baseline mesh and anteroinferior defect meshes
    glenoidTests = ['intact'; anteriorDefectsList];
    
    %Create list of glenoid rotations to test
    %(i.e. different translation directions)
    %Currently 3 o clock, 4 o clock, and 5 o clock
    glenoidRotations = [0; 30; 60];

    %Loop through the different meshes
    for gg = 1:size(glenoidTests,1)
        
        %Loop through glenoid rotation
        for tt = 1:size(glenoidRotations,1)
            
            %Determine clock face based on rotation
            rotationVal = glenoidRotations(tt);
            switch rotationVal
                case -30
                    clockPos = 2;
                case 0
                    clockPos = 3;
                case 30
                    clockPos = 4;
                case 60
                    clockPos = 5;
            end
            
            %Create a simulation name based on different components
            simName = [glenoidTests{gg},'_',num2str(clockPos),'oClockTranslate'];
                
            %Create a directory to store current simulation data
            mkdir(simName); cd(simName);

            %Create the FEBio run file for the distance to
            %dislocation testing
            simType = 'translation';
            if gg == 1
                %Use the baseline glenoid mesh structure
                [feaMeshOutputs.(char(simName)).glenoidMeshOutput,...
                    feaMeshOutputs.(char(simName)).headMeshOutput] = createFEBioRunFile_congruentArcVstandardLatarjet(glenoidMesh,...
                    simName,glenoidRotations(tt),humerusCS,shapes,simType,generatePlots);
            else
                %Use the defect glenoid mesh structure
                [feaMeshOutputs.(char(simName)).glenoidMeshOutput,...
                    feaMeshOutputs.(char(simName)).headMeshOutput] = createFEBioRunFile_congruentArcVstandardLatarjet(anteriorDefectsMesh.(glenoidTests{gg}),...
                    simName,glenoidRotations(tt),humerusCS,shapes,simType,generatePlots);
            end

            %Set FEBio analysis details
            febioAnalysis.run_filename = [simName,'_',simType,'.feb']; %The input file name
            febioAnalysis.run_logname = [simName,'_',simType,'.txt']; %The name for the log file
            febioAnalysis.disp_on = 1; %Display information on the command window
            febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
            febioAnalysis.runMode = 'internal';
            febioAnalysis.t_check = 0.25; %Time for checking log file (dont set too small)
            febioAnalysis.maxtpi = 1e99; %Max analysis time
            febioAnalysis.maxLogCheckTime = 60; %Max log file checking time
            %Specifying path required on lab PC for some reason???
            %Also works on laptop anyway...
            febioAnalysis.FEBioPath = 'C:\Program Files\FEBio2.8.3\bin\FEBio2.exe';

            %Run the simulation in FEBio
            clc
            %%%%% TODO: write up display outputs better to track progress
            [runFlag(gg)] = runMonitorFEBio(febioAnalysis);

            %Create the FEBio run file for the force simulation
            simType = 'force';
            if gg == 1
                %Use the baseline glenoid mesh structure
                [feaMeshOutputs.(char(simName)).glenoidMeshOutput,...
                    feaMeshOutputs.(char(simName)).headMeshOutput] = createFEBioRunFile_congruentArcVstandardLatarjet(glenoidMesh,...
                    simName,glenoidRotations(tt),humerusCS,shapes,simType,generatePlots);
            else
                %Use the defect glenoid mesh structure
                [feaMeshOutputs.(char(simName)).glenoidMeshOutput,...
                    feaMeshOutputs.(char(simName)).headMeshOutput] = createFEBioRunFile_congruentArcVstandardLatarjet(anteriorDefectsMesh.(glenoidTests{gg}),...
                    simName,glenoidRotations(tt),humerusCS,shapes,simType,generatePlots);
            end

            %Set FEBio analysis details
            febioAnalysis.run_filename = [simName,'_',simType,'.feb']; %The input file name
            febioAnalysis.run_logname = [simName,'_',simType,'.txt']; %The name for the log file
            febioAnalysis.disp_on = 1; %Display information on the command window
            febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
            febioAnalysis.runMode = 'internal';
            febioAnalysis.t_check = 0.25; %Time for checking log file (dont set too small)
            febioAnalysis.maxtpi = 1e99; %Max analysis time
            febioAnalysis.maxLogCheckTime = 60; %Max log file checking time
            %Specifying path required on lab PC for some reason???
            febioAnalysis.FEBioPath = 'C:\Program Files\FEBio2.8.3\bin\FEBio2.exe';

            %Run the simulation in FEBio
            clc
            %%%%% TODO: write up display outputs better to track progress
            [runFlag(gg)] = runMonitorFEBio(febioAnalysis);

            %Navigate back up a directory
            cd('..');
                    
        end
        clear tt
       
    end
    clear gg
    
    %% Import FEBio results
    
    %Set flag to export animated gifs
    %Note that this adds a bit of time to the process function. It also
    %generates a number of warnings despite wanting these to be turned
    %off.
    exportAnimation = true;
    
    %Loop through glenoid tests variable
    for gg = 1:size(glenoidTests,1)
        
        %Check for successful run
        
        %%%% TODO: need better characterisation of runFlag variable
        
        if runFlag(gg) == 1

            %Loop through glenoid rotation
            for tt = 1:size(glenoidRotations,1)
                
                %Determine clock face based on rotation
                rotationVal = glenoidRotations(tt);
                switch rotationVal
                    case -30
                        clockPos = 2;
                    case 0
                        clockPos = 3;
                    case 30
                        clockPos = 4;
                    case 60
                        clockPos = 5;
                end

                %Create corresponding simulation name for current iteration
                simName = [glenoidTests{gg},'_',num2str(clockPos),'oClockTranslate'];

                %Navigate to run directory
                cd(simName);

                %Run process function to collate results from both simulations 
                [processedResults.(char(simName))] = processFEBioResults(simName,...
                    feaMeshOutputs.(char(simName)).glenoidMeshOutput,...
                    feaMeshOutputs.(char(simName)).headMeshOutput,...
                    landmarks,exportAnimation,generatePlots);

                %Return up a directory
                cd('..');

            end
            clear tt
            
        end
    end
    clear gg
    
    %% Compile simulation results
    
    
    
    
    %%
    
    %%

end
