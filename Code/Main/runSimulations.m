function runSimulations()

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
    % Klemt et al. (2019). The critical size of a defect in the glenoid
    % causing anterior instability of the shoulder after a Bankart repair,
    % under physiological joint loading. Bone Joint J, 101-B: 68-74.

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
    
    %Visualise anterior defects on a subplot
    anteriorDefectsList = fieldnames(anteriorDefectsMesh);
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
    
    %%%%%% TODO: tetgen on anterior defect glenoids
    
    %%
    
    %% Run baseline FEBio simulation
    
    %%%%% TODO: appropriate spot to store FEBio run file
    
    createFEBioRunFile()
    
    %%%%%% TODO: set this up better...
    
    febioAnalysis.run_filename = 'baselineSim.feb'; %The input file name
    febioAnalysis.run_logname = 'baselineSim.txt'; %The name for the log file
    febioAnalysis.disp_on = 1; %Display information on the command window
    febioAnalysis.disp_log_on = 1; %Display convergence information in the command window
    febioAnalysis.runMode = 'external';%'internal';
    febioAnalysis.t_check = 0.25; %Time for checking log file (dont set too small)
    febioAnalysis.maxtpi = 1e99; %Max analysis time
    febioAnalysis.maxLogCheckTime = 60; %Max log file checking time
    
    %%%%% log checking seems to take a lomg time on laptop --- it's having
    %%%%% a slow day though...although it worked quicker another time...
    
    %%%%%% may need to setup FEBio better on the lab computer...
    
    clc
    [runFlag] = runMonitorFEBio(febioAnalysis);%START FEBio NOW!!!!!!!!
    
    
    %%%%% Surprisingly, this all works. Need to therefore work on the FEBio
    %%%%% creation function so that it is suitably adaptable to the
    %%%%% different scenarios.
    
    
    
    %%
    %%
    %%   
    %% Import FEBio results
    
    %%%%% TODO: edit appropriately, cleanup to display better...
    %%%%% Needs to factor in which displacements belong to which nodes...
    %%%%% Perhaps it would be best to separate the displacement collection
    %%%%% into two files?

    %Check for successful run
    if runFlag == 1

        %Importing nodal displacements from a log file
        %%%%% TODO: comment this better and adjust variable names
        %%%%% TODO: might be a better way to extract results too...

        %%%%% TODO: anim8 is not working right, there is some disconnect
        %%%%% between the elements it's animating and the data...

        %Import the output logfile for displacements
        [time_mat, N_disp_mat,~] = importFEBio_logfile('baselineSim_disp_out.txt'); %Nodal displacements

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

end
