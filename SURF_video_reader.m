function rosette_tracker_part_1 = SURF_video_reader(folder, video_filename)
rosette_tracker_part_1 = 0;
Videostored = VideoReader(video_filename );
video_filename = Videostored.Name;
image_height = Videostored.Height;
image_width = Videostored.Width;
image_count = Videostored.NumberOfFrames;
ImageArray = cell(image_count, 1);
for i=1:image_count
    img_i = read(Videostored,i);
    ImageArray{i} = img_i;
end
%%calculate background
filecount = 1000; %bei ca. 100.000 frames
bg = zeros(image_height, image_width);
vals = zeros(1,filecount);
for a=1:image_width
    for b=1:image_height
        for i=1:filecount
            vals(i) = ImageArray{i}(b,a);
        end
        bg(b,a) = int8(median(vals));
    end
end

figure(1)
imshow(uint8(bg));
outputname = [folder 'BG.png'];
imwrite(uint8(bg),  outputname);
black_bg = zeros(image_height, image_width);

%% find object locations
object_locations = cell(image_count, 1);
for n=1:image_count
    img = round(255*im2double(ImageArray {n}));
    img_subtr_med = (((img) - ( bg ))); %subtract background
    img_rescaled=uint8(abs(img_subtr_med ));
    img_smooth = imguidedfilter(img_rescaled);
    ImageArray {n} = img_smooth;
end

%% find tresh value for contour plot
tresh_value = 0;
for k = 1:10
    figure (2)
    [C, h] = imcontour(ImageArray {k}, 1); 
    [x,y,z] = C2xyz(C);
    tresh_value = tresh_value+z(1);
end
tresh_value =0.75* tresh_value/10;
%% use contour plot to draw binary images and find objects
for n = 1:image_count
    img_bw = black_bg;
    C = contourc(double(ImageArray {n}), [tresh_value tresh_value]); 
    if ~isempty(C)
        [x,y,~] = C2xyz(C);
        for a = 1:size(x,2) %number of contours on image
            for b = 1:size(x{a},2) %length of contour
                img_bw(round(y{a}(b)), round(x{a}(b))) = 255;
            end
        end
    end
    img_clear = imclearborder(img_bw);
    img_filled = imfill( img_clear, 'holes');
    img_removed = bwareaopen(img_filled, 20);%Remove objects containing fewer than 20 pixels using bwareaopen function.
    ImageArray {n} = img_removed;
    s = regionprops (img_removed, 'Area', 'Centroid', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation');
    anzahl_objekte_bildi = size(s);
    for j=1:anzahl_objekte_bildi(1)
        object_locations{n}(j,1)= n;%Frame
        object_locations{n}(j,2)= j;%Objektnummer
        object_locations{n}(j,3)=round(s(j).Centroid(1));%x coordinate
        object_locations{n}(j,4)=round(s(j).Centroid(2));%y coordinate
        object_locations{n}(j,5)=s(j).Area;%area in pixels
        object_locations{n}(j,6)=s(j).MajorAxisLength; %surrounding ellipse
        object_locations{n}(j,7)=s(j).MinorAxisLength;
        object_locations{n}(j,8)=s(j).Eccentricity; 
        object_locations{n}(j,9)=s(j).Orientation;
        
    end
end

%% text file object locations
ObjectLocationsFileName = [folder 'SURF_object_locations.txt'];
fid = fopen(ObjectLocationsFileName,'w');
if fid ~= -1
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\t%s%s\r\n','results rosette tracker', folder , video_filename);
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', '#frame','#object','x Pos (pixel)','y Pos (pixel)','Area', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity', 'Orientation');
    fprintf(fid,'%s\r\n', '===...===');
    for n = 1:image_count
        for k = 1:(size( object_locations{n}, 1))
            fprintf(fid,'%s\t%s\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\r\n',num2str(object_locations{n}(k,1)),num2str(object_locations{n}(k,2)),num2str(object_locations{n}(k,3)),num2str(object_locations{n}(k,4)),num2str(object_locations{n}(k,5)),num2str(object_locations{n}(k,6)),num2str(object_locations{n}(k,7)),num2str(object_locations{n}(k,8)),num2str(object_locations{n}(k,9)));
        end
    end
end
fclose(fid);
rosette_tracker_part_1 = 1;
