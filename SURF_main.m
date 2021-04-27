directory = "E:\Stockholm March 2019\BGA\";
folders = [ "BGA 7mum Standard 25mulh_20190318_152513\";
"BGA 5mum Standard 25mulh_20190319_133319\";
...
"BGA 5mum 3EL 25mulh_20190319_150141\"];
video_filenames = [ "BGA 7mum Standard 25mulh_20190318_152513.avi";
"BGA 5mum Standard 25mulh_20190319_133319.avi";
...
"BGA 5mum 3EL 25mulh_20190319_150141.avi"];

for i = 1:size(folders,1)
folder =directory + folders(i);
folder = char(folder);
video_filename = directory + folders(i) + video_filenames (i);
video_filename = char(video_filename);
t = datetime('now');
disp(['## Start Videofile ' num2str(i) ': ' datestr(t)]);
disp ([ video_filename ]);

video_read_successful = SURF_video_reader_nofig3(folder, video_filename);
if(video_read_successful)
    t = datetime('now');
    disp(['Videoread successful ' datestr(t)]);
end
text_read_successful = SURF_textfile_reader_(folder, video_filename);
if(text_read_successful)
    t = datetime('now');
    disp(['Textfileread successful ' datestr(t)]);
end
close all
end