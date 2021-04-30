# SURF SUrvival of Rosettes in Flow
Tracking Algorithm for cell aggregates


If you need to trace assembled objects along a channel and want to detect change in shape and size and loss of components, SURF may help you. Here we describe a MATLAB based tracking algoritm named **SURF** (**SU**rvival of **R**osettes in **F**low). The script is designed to analyze high speed micrographs of cell aggregates passing a microfluidic channel. A standard measurement generates a video of 74531 frames measuring 1280 x 120 Pixels at 2000 frames per second. You can download an [example video] to run the code.
 
 [![bsptracemitQR](https://user-images.githubusercontent.com/83273863/116454640-819fea80-a860-11eb-9f53-be79e97065e2.png)](https://youtu.be/rmezGUd0p08)
The image shows an overlay of a few examplary frames. SURF reads the videos, indentifies cell aggregates and follows them along the channel, to analyze their deformation and events as the loss of a cell.  All figures here are taken from an exemplary data set to visualize the function of the code. The exemplary measurement uses a cascadic channel with three stenoses of decreasing diameter of  d=11 µm, 7 µm und 5 µm and rosetting red blood cells from blood group A. %BSP Ergebnis Trajektorien im Kaskadenkanal BGA 25 mulhr 20180928 155155.avi
  
 
## Application: Rosetting in Malaria
In Malaria pathogenesis, the role of rosetting, the formation of red blood cell aggregates,  is unclear. We study rosette stability using a set of microfluidic stenotic channels, and analyze the influence of the shape of the stenosis. The results created using the algorithm SURF are published in [Biomicrofluidics 2020].  We find reduced ability of a rosette to pass a stenosis, the larger the rosette and the longer and the narrower the stenosis is. 


## Tracing
![verfolgung](https://user-images.githubusercontent.com/83273863/116463278-ce88be80-a86a-11eb-89aa-d0ff305263a4.png)
The main function SURF_main.m contains a list of the file paths of the videos of interest. For each video first the function SURF_video_reader.m is called, which returns a textfile, that contains a tabular list of every object and its size and position on every frame. This textfile is the input to the function SURF_textfile_reader.m, which contains the actual tracking algorithm and identifies and connects rosettes from frame to frame. 

## Events
The preprocessing of the trajectories includes the identification of "cell loss" or "cell gain" events. "cell loss" is further distinguished into "real rupter" and "rupture and reconnect", "cell gain" can be "connect only" or "pass by". The following scheme shows how the events are categorized.
![scheme](https://user-images.githubusercontent.com/83273863/116463192-b1ec8680-a86a-11eb-9730-c1050970d5d9.png)

# How to use the code
- Download the [example video]
- Start with SURF_main.m
- Adjust path and filename of the video
- Open SURF_textfile_reader.m and adjust timesteps (time delay between frames in the video) and measurements of the channel for correct calculation of the flow rate

Running SURF_main.m will call the functions SURF_video_reader.m and SURF_textfile_reader.m in that order.

The results will be exported to the folder, which contains the video. Each "trace" will be saved as a textfile and a plot. The events will be listed in event_counter_grid.txt.

## Further post-processing
The SURF_trace n.txt-files and event_counter_grid.txt are used for further post-processing by the remaining three scripts:
### SURF_sort_traces_by_rosette_size.m
- Adjust directory and folders
- will sort traces by rosette size, calculate flowrate in each size class and export statistics of events sorted by rosette size.
- rosettes_by_size.txt contains counts, relative frequencies, and a list of normalized sizes (size in pixel divided by single cell size in pixel) of all rosettes in each class
- event_grid_by_size.txt contains tables for each kind of event. the six lines of the tables correspond to the six classes of rosettes, the columns of the table correspond to the x-coordinates where the events were registered.

### SURF_sort_traces_by_rosette_fate.m
- Adjust directoy, folders and add x-location in pixels of the end of the stenosis.
- will sort by rosette size, as before, and calculate what events happenend before and after the "end of the stenosis" x-coordinate. This result can also be read manually from event_grid_by_size.txt, which is exported from SURF_sort_traces_by_rosette_size.m.
- what is extra is damaged_area.txt, the export of lost area or lost cells per rosette.

### SURF_elastic_modulus.m
- Adjust directoy, folders and add x-location in pixels of the end of the stenosis.
- Adjust length and width of the stenosis.
- calculates velocity and deformation in each class for the undamaged rosettes only.
- plotting stress (calculated from velocity) against deformation, gives the elastic modulus of the deformed object.






[Biomicrofluidics 2020]: https://doi.org/10.1063/1.5125038
[example video]: https://drive.google.com/file/d/1EDqv4EtH839AH-NfndiahIGv0mBsXkVL/view?usp=sharing
