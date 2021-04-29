# SURF SUrvival of Rosettes in Flow
Tracking Algorithm for cell aggregates


If you need to trace assembled objects along a channel and want to detect change in shape and size and loss of components, SURF may help you. Here we describe a MATLAB based tracking algoritm named **SURF** (**SU**rvival of **R**osettes in **F**low). The script is designed to analyze high speed micrographs of cell aggregates passing a microfluidic channel. A standard measurement generates a video of 74531 frames measuring 1280 x 120 Pixels at 2000 frames per second. You can download an example video to run the code here: https://drive.google.com/file/d/1EDqv4EtH839AH-NfndiahIGv0mBsXkVL/view?usp=sharing
 
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
*Start with SURF_main.m
*Adjust path and filename of the video
Running SURF_main.m will call the functions SURF_video_reader.m and SURF_textfile_reader.m in that order.






[Biomicrofluidics 2020]: https://doi.org/10.1063/1.5125038
