%adjust file path
directory = "C:\Users\joettean\Desktop\Paper AJ 2021\SURF Code\";
folders = [ "SURF_example_video\";
];
video_filenames = [ "SURF_example_video.avi";
];

for i = 1:size(folders,1)
folder =directory + folders(i);
folder = char(folder);
video_filename = directory + folders(i) + video_filenames (i);
video_filename = char(video_filename);
t = datetime('now');
disp(['## Start Videofile ' num2str(i) ': ' datestr(t)]);
disp ([ video_filename ]);

video_read_successful = SURF_video_reader(folder, video_filename);
if(video_read_successful)
    t = datetime('now');
    disp(['Videoread successful ' datestr(t)]);
end
text_read_successful = SURF_textfile_reader(folder, video_filename);
if(text_read_successful)
    t = datetime('now');
    disp(['Textfileread successful ' datestr(t)]);
end
close all
end