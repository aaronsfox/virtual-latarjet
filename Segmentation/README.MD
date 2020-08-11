# Segmentation & Surface Meshing Processes

The following sections outline the details and processes used for segmenting and extracting surface meshes from the [NMDID Dataset](https://nmdid.unm.edu/). Each participant has their own folder which includes the Mimics segmentation file, the 3matic project file, and the extracted surface meshes of the humeru and scapula.

## Segmentation in Mimics

Segmentation of the CT data was performed in Materialise Mimics Research v21.0. The major tool used was the threshold feature to create masks of the bony elements of the scans. Mimics splitting tool was used to extract indiviual parts. The various filling tools were used to fill out the masks and create subsequent complete surface parts of the relevant bodies. Manual refinement of the masks was performed using the other editing tools (e.g. threshold painting) to finalise the parts.

The masks of the scapula and humerus were converted to parts and copied across to a 3matic project file, where some further processing (see notes below) were undertaken before exporting as an ASCII STL file.

## Surface Meshing in 3matic

Some processing of the meshes prior to exporting was performed in 3matic. This included some surface smoothing as well as marking of landmarks for processing later on. Materialise 3matic Research v13.0 was used to edit and export surface meshes. 

### Meshing Steps

The following section outlines the steps for creating the meshes used for analysis:

- The 3matic 'Smooth' feature was used with a smooth factor of 0.7 to smooth the scapula and humerus surfaces.
- The 3matic 'Adaptive Remesh' feature was used to remesh the triangles of the surfaces. A target edge length of 0.5mm and maximum edge length of 1.0mm were used while maintaining surface contours. This resulted in quite a fine mesh but this will be re-assessed later with mesh refinement processes.
- A sphere was fit to the humeral head ('HeadSphere') by marking selected triangles with the brush mark and fitting an analytical sphere to these. The triangles were copiedand retained in a new surface part to be exported later.
- A series of landmarks were marked according to the definitions provided by Wu et al. (2005).
- A 'GHJ' point was created at the centre of the aforementioned sphere.
- The epicondyle points of the humerus were marked and labelled as 'EL' and 'EM'.
- The angulus acromialis ('AA'), trigonum scapulae ('TS') and angulus inferior ('AI') were marked on the scapula.
- The supraglenoid ('SGT') and infraglenoid ('IGT') tubercles of the scapula were marked.
- A plane was fit to the glenoid surface ('GlenoidPlane') by marking triangles on the glenoid face. These triangles were copied and retained as a separate part.
- This plane was used to identify and create a point at the central deepest point of the glenoid surface ('DeepGlenoidPoint').
- The humerus and scapula were exported as ASCII .stl files.
- All landmarks and objects were exported as 'Curves and Primitives' to separate .xml files.
- These files were used as inputs in a subsequent Matlab code (see: ***TODO: add .m filename***).

## References

Wu et al. (2005). ISB recommendation on definitions of joint coordinate systems of various joints for the reporting of human joint motion - Part II: Shoulder, elbow, wrist and hand. *J Biomech*, 38: 981-992.