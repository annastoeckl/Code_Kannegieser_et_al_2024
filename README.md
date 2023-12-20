# Code_Kannegieser_et_al_2024
Analysis files for the publication "Visual guidance fine-tunes probing movements of an insect appendage" , DOI follows soon. Processes body and proboscis positions of hawkmoths probing artificial flower patterns.  
The raw data for this publication is deposited at https://doi.org/10.6084/m9.figshare.22639981.v1.

The scripts and functions in this repository build an analysis pipeline, to extract and process the data shown in the manuscript from csv files containing the curated tracking data (for further details, see Manuscript Methods). 

# Data extraction and pre-processing

0) readH5fromDLC reads the tracking points from DeepLabCut, which are packaged in an .h5 file, and extracts the x/y positions of the points, as well as their labels, to feed them into a .csv file which we can further process with the Hedrick's Lab's dltv software for data point curation, and for marking the flower and pattern outlines https://biomech.web.unc.edu/dltdv/

1) load_excel_file_proboscis_track.m extracts the proboscis, head and thorax coordinates, as well as outline of the flower in the videos, and that of potential flower patterns, and stores these as Matlab variables in an output .mat file
1a) load_excel_file_proboscis_track_head.m does the same as the above, but also includes three more head keypoints, the base of the proboscis, left and right antenna.

2) append_and_rotate_data loads the tracked data from the .mat files generated in step 1, adds together data from multiple files (if multiple trials were recorded for the same animal and condition), and rotates the data to align the pattern's main axis with the cardinal coordinates. In addition to identifying trials from multiple files, it also identifies individual approaches of the moths to the flower within one datafile (videofile), as flight bouts. It saves the output as allRotData_filename.mat, which provides the input for all subsequent analysis
2a) append_and_rotate_data_heads does the same as the above, but also includes three more head keypoints, the base of the proboscis, left and right antenna.

The data repository contains these rotated and annotated data files.

3) plot_flower_animal_track plots the head, thorax and proboscis positions in one data file (uses output from load_excel_file_proboscis_track or append_and_rotate_data) with a colour code to denote time.

# Temporal movement analysis of body and proboscis

1) analyseProboscisMotion analyse the movement of the proboscis and the body in relation to each other, and to the pattern. It generates Fast Fourier Transforms of body and proboscis movements, relative to the pattern axes, and relative to the animals' body axis (see Methods of Manuscript for more details on this metric).
1a) analyseProboscisMotion_heads does the same as the above, but also includes three more head keypoints, the base of the proboscis, left and right antenna. Uses these to calculate coherence between movement of different head keypoints, and head angle relative to the body, and proboscis movement.

# Spatial movement analysis of body and proboscis
Two functions with partly overlapping functionality were used in this analysis
1) plot_heat_map analyses the distribution of proboscis contacts on the flower, relative  to the pattern. It plots a heatmap of proboscis contacts and generates histograms of the contact distributions along the cardinal axes of the pattern, and in the outer and inner thirds of the pattern. It also calculates contact scores, which count and compare the number of proboscis contacts along either cardinal axis of the flower and pattern, as well as within the inner and outer thirds of the pattern. For line patterns, it compares these scores to a hypothetical pattern rotated 90 degrees from the originial, for cross patterns 45 degrees rotated.
2) plot_heat_map_proboscis_track analyses the distribution of proboscis contacts on the flower, relative  to the pattern. 
It plots a heatmap of proboscis contacts and the paths as line plots. 
It also calculates contact scores, which count and compare the number of proboscis contacts along either cardinal axis of the flower and pattern, and within the pattern and on the background, and includes scales for the area of pattern and background. 
From the proboscis tracks, the duration, speed and tortuosity of probing are extracted and plotted, as well as the mean direction and vector length of the probing tracks, relative to the orientation of the pattern.
Moreover, the body orientation (defined by the head-thorax axis) and position of the proboscis, head and thorax relative to the long axis of a line pattern are extracted and plotted.


# Analysis of flight relative to pattern

1) flight_orientation extracts key descriptors of the animals' flight during the approach, the probing phase and the departure. These are forward and rotational flight speed and predominant rotation direction, as well as body angle relative to the pattern axis. The output is saved as a flightORI_filename.mat file.
2) Uses the output of the above script, and plots the descriptors of flight across all selected data files (generally used to plot data from multiple animals within one condition).

# Summary comparisons and statistics

1) compAcrossConds plots and statistically compares individually saved proboscis data from separate conditions (data generated with plot_heat_map_proboscis_track and analyseProboscisMotion)
2) compAcrossConds_HW plots and statistically compares individually saved data of halfwidth of proboscis contact distributions from separate conditions (data generated with plot_heat_map_proboscis_track)
