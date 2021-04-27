function rosette_tracker_part_2 = SURF_textfile_reader(folder, video_filename)
rosette_tracker_part_2 = 0;
%%read txtfile object locations
file_objectlocations = sprintf('%sobject_locations.txt', folder);
objectlocations_complete = dlmread(file_objectlocations,'', 5, 0);
nr_of_frames =  objectlocations_complete(end,1); %number of current trace
object_locations = cell(nr_of_frames, 1);
for k = 1:nr_of_frames
    indexlist_k = find(objectlocations_complete(:,1)==k);
    if ~isempty(indexlist_k)
        object_locations{k}=objectlocations_complete(indexlist_k(1):indexlist_k(end),1:end);
    end
    if isempty(indexlist_k)
        object_locations{k}=0*objectlocations_complete(1,1:end);
    end
end
%%read background image
background_filename = sprintf('%sBG.png', folder);
bg = imread(background_filename);
dt = 0.0005; %delay frames (s)
image_height = size(bg, 1);
image_width = size(bg, 2);
canal_depth = 8*10^-6;
meter_per_pixel = 0.55*10^-6; %neue Photron 20x: 90Px = 50mum
canal_crosssection = 90*meter_per_pixel*canal_depth;

%% reduce to large objects (rosettes = >1.5* median_object_area) and define input_delay_x
rosette_locations = cell(nr_of_frames, 1);
input_delay_x=50;
medians_object_area = zeros(nr_of_frames,1);
for n=1:nr_of_frames
    medians_object_area (n) = median(object_locations{n}(:,5));
end
median_object_area = median (medians_object_area);
for n = 1:nr_of_frames
    rosette_locations{n}(:,1)=object_locations{n}(:,3);
    rosette_locations{n}(:,2)=object_locations{n}(:,4);
    rosette_locations{n}(:,3)=object_locations{n}(:,5);
    rosette_locations{n}(:,4)=object_locations{n}(:,1);
    objects_before = size (object_locations{n},1);
    for j=objects_before:-1:1
        if object_locations{n}(j,5) <= 1.5* median_object_area || object_locations{n}(j,3) < input_delay_x
            rosette_locations{n}(j,:) = [];
            
        end
    end
end

%% TRACING
traces = cell(nr_of_frames, 1);
bridges = cell(nr_of_frames, 1);
%target area in flow direction
target_dy = 45;
target_delay_x = 0; %minimum displacement in x (good for large displacements per frame)
max_distance =150; %maximum displacement in x (good for small displacements per frame)

for rosette_n = 1:(nr_of_frames-10)
    frame_n = rosette_n;
    if ~isempty(rosette_locations{frame_n}) %find first object
        if  rosette_locations{frame_n}(1,1) < 1.5*input_delay_x %start of trace always between input_delay_x and 1.5*input_delay_x
            traces{rosette_n}(frame_n-rosette_n+1,:) = rosette_locations{frame_n}(1,:); %set first coodinates
        else
            continue
        end
    else
        continue
    end
    if rosette_n > 1
        %ignore trace if same as previous trace
        if ~isempty(traces{rosette_n-1})
            if size(traces{rosette_n-1},1) ==1
                continue
            end
            if traces{rosette_n}(1,1)==traces{rosette_n-1}(2,1) && traces{rosette_n}(1,2)==traces{rosette_n-1}(2,2)
                continue
            end
        end
    end
    while 1
        if isempty(rosette_locations{frame_n+1}) %break if no rosettes in following frame
            break
        end
        max_x_coord_next_frame_index = size (rosette_locations{frame_n+1},1);
        max_x_coord_next_frame = rosette_locations{frame_n+1}(max_x_coord_next_frame_index,1);
        if max_x_coord_next_frame >  traces{rosette_n}(frame_n-rosette_n+1,1)- 2*max_distance %find rosettes in target region
            rosettes_downstream = rosette_locations{frame_n+1};
            rosettes_left_downstream = rosette_locations{frame_n+1};
            %delete rosettes to the left
            for k = size(rosettes_downstream,1):-1:1
             %reduce to rosettes strict downstream
                if rosettes_downstream(k,1) < (traces{rosette_n}(frame_n-rosette_n+1,1)+target_delay_x) || (rosettes_downstream(k,1) > traces{rosette_n}(frame_n-rosette_n+1,1) + max_distance)
                    rosettes_downstream(k,:) = [];
                end
            end
            for k = size(rosettes_left_downstream,1):-1:1 
            %reduce to rosettes 2*sqrt(traces{rosette_n}(frame_n-rosette_n+1,3)/pi())) upstream
                if (rosettes_left_downstream(k,1) < (traces{rosette_n}(frame_n-rosette_n+1,1)+target_delay_x - 2*sqrt(traces{rosette_n}(frame_n-rosette_n+1,3)/pi()))) || (rosettes_left_downstream(k,1) > (traces{rosette_n}(frame_n-rosette_n+1,1)+target_delay_x ))
                    rosettes_left_downstream(k,:) = [];
                end
            end
            for k = size(rosettes_downstream,1):-1:1
                if  abs(rosettes_downstream(k,2)-traces{rosette_n}(frame_n-rosette_n+1,2)) > target_dy
                    rosettes_downstream(k,:) = [];
                end
            end
            for k = size(rosettes_left_downstream,1):-1:1
                if  abs(rosettes_left_downstream(k,2)-traces{rosette_n}(frame_n-rosette_n+1,2)) > target_dy
                    rosettes_left_downstream(k,:) = [];
                end
            end
            
            distances_downstream = zeros(size(rosettes_downstream,1),1);
            for i = 1:size(rosettes_downstream,1)
                distances_downstream(i) = (rosettes_downstream(i,1) - traces{rosette_n}(frame_n-rosette_n+1,1))^2+(rosettes_downstream(i,2) - traces{rosette_n}(frame_n-rosette_n+1,2))^2;
            end
            [Min_value_downstream, Min_index_downstream] = sort(distances_downstream(distances_downstream>0));
            
            distances_left_downstream = zeros(size(rosettes_left_downstream,1),1); 
            for i = 1:size(rosettes_left_downstream,1)
                distances_left_downstream(i) = (rosettes_left_downstream(i,1) - traces{rosette_n}(frame_n-rosette_n+1,1))^2+(rosettes_left_downstream(i,2) - traces{rosette_n}(frame_n-rosette_n+1,2))^2;
            end
            [Min_value_left_downstream, Min_index_left_downstream] = sort(distances_left_downstream);
            
            areas_left_downstream = zeros(size(rosettes_left_downstream,1),1); %absolute difference of area
            for i = 1:size(rosettes_left_downstream,1)
                areas_left_downstream(i) = abs(rosettes_left_downstream(i,3) - traces{rosette_n}(frame_n-rosette_n+1,3))/traces{rosette_n}(frame_n-rosette_n+1,3);
            end
            [Min_value_areas, Min_index_areas] = sort(areas_left_downstream);  
            %Min_value too far down stream ?
            if ~isempty(Min_value_downstream)
                if Min_value_downstream(1) > max_distance^2
                    %do not lose rosette to the left at end of trace
                    if ~isempty(Min_value_left_downstream)
                        if Min_value_left_downstream > max_distance^2
                            break
                        end
                    end
                    if ~isempty(Min_index_left_downstream)
                        traces{rosette_n}(frame_n+1-rosette_n+1,:) = rosettes_left_downstream((Min_index_left_downstream(1)),:);
                        bridge_n = size (bridges{rosette_n},1)+1;
                        bridges{rosette_n}(bridge_n,:) = rosettes_left_downstream((Min_index_left_downstream(1)),:);
                    end
                    if isempty(Min_index_left_downstream)
                        break
                    end
                end
            end
            %matching rosette found >> add to trace:
            if ~isempty(Min_index_downstream)
                traces{rosette_n}(frame_n+1-rosette_n+1,:) = rosettes_downstream((Min_index_downstream(1)),:);
            end
            %change in area too big? find bulk rosette instead of offspring
            %(smallest area difference)
            if  ~isempty(Min_value_downstream) && ~isempty(Min_value_left_downstream)
                if abs (traces{rosette_n}(frame_n+1-rosette_n,3) - rosettes_downstream((Min_index_downstream(1)),3)) / traces{rosette_n}(frame_n+1-rosette_n,3) > 0.5 && Min_value_areas(1) < 0.5
                    traces{rosette_n}(frame_n+1-rosette_n+1,:) = rosettes_left_downstream((Min_index_areas(1)),:);
                    bridge_n = size (bridges{rosette_n},1)+1;
                    bridges{rosette_n}(bridge_n,:) = rosettes_left_downstream((Min_index_left_downstream(1)),:);
                end
            end
            %no rosette found down stream ?
            if isempty(Min_index_downstream)
                %do not lose rosette to the left at end of trace
                if ~isempty(Min_value_left_downstream)
                    if Min_value_left_downstream(1) > max_distance^2
                        break
                    end
                end
                if ~isempty(Min_index_left_downstream)
                    traces{rosette_n}(frame_n+1-rosette_n+1,:) = rosettes_left_downstream((Min_index_left_downstream(1)),:);
                    bridge_n = size (bridges{rosette_n},1)+1;
                    bridges{rosette_n}(bridge_n,:) = rosettes_left_downstream((Min_index_left_downstream(1)),:);
                end
                if isempty(Min_index_left_downstream)
                    break
                end
            end
        else
            break
            
        end
        if frame_n == (nr_of_frames-1)
            break
        end
        frame_n = frame_n+1;
        
    end
end

%%filter traces
min_trace_length = 5; %trace_length > min_trace_length 
number_vel_points = 0;
number_filtered_traces = 0;
for n = 1:rosette_n
    if size(traces{n},1)>min_trace_length
        number_vel_points = number_vel_points + size(traces{n},1) - 1;
        number_filtered_traces = number_filtered_traces +1;
    end
end
%%velocity matrix
vel_points = zeros(number_vel_points, 10);

%%Output All traces on Background
outputname_traces = [folder 'traces.png'];
figure (4)
clf;
imshow(uint8(bg));
title(sprintf('%d traces (length > %d )',  number_filtered_traces, min_trace_length));
hold on
for n = 1:rosette_n
    if size(traces{n},1)>min_trace_length
        plot(traces{n}(:,1),traces{n}(:,2),'LineStyle', '-', 'Marker', '.');
    end
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_traces);

traces_filtered = cell(number_filtered_traces, 2);
bridges_filtered = cell(number_filtered_traces, 1);
k = 1;
for n = 1:rosette_n
    if size(traces{n},1)>min_trace_length
        traces_filtered{k,1}=traces{n};
        for j = 1:size(traces{n},1)
            l = find(object_locations{traces{n}(j,4)}(:,3)==traces{n}(j,1));
            if size(l,1)==1 %if there is only one object with matching x value
                traces_filtered{k,1}(j,5)=object_locations{traces{n}(j,4)}(l,6);
                traces_filtered{k,1}(j,6)=object_locations{traces{n}(j,4)}(l,7);
                traces_filtered{k,1}(j,7)=object_locations{traces{n}(j,4)}(l,8);
                traces_filtered{k,1}(j,8)=object_locations{traces{n}(j,4)}(l,9);
                traces_filtered{k,1}(j,9)=0;
            end
            if size(l,1)>1 %if there is more than one object with matching x value, check y
                for m = 1:size(l,1)
                    if object_locations{traces{n}(j,4)}(l(m),4)==traces{n}(j,2)
                        traces_filtered{k,1}(j,5)=object_locations{traces{n}(j,4)}(l(m),6);
                        traces_filtered{k,1}(j,6)=object_locations{traces{n}(j,4)}(l(m),7);
                        traces_filtered{k,1}(j,7)=object_locations{traces{n}(j,4)}(l(m),8);
                        traces_filtered{k,1}(j,8)=object_locations{traces{n}(j,4)}(l(m),9);
                        traces_filtered{k,1}(j,9)=0;
                    end
                end
            end
        end
        bridges_filtered{k} = bridges{n};
        k=k+1;
    end
end

%% neighbouring objects
surroundings = cell(number_filtered_traces, 1); %
for n = 1:number_filtered_traces
    end_x = size( traces_filtered{n}, 1);
    surroundings{n} = zeros(end_x, 5); %matrix to save number of objects in sourroundings of rosette now, previous picture and difference
    for k = 1:end_x-1
        %frame traces_filtered{n}(k,4)
        %x positions column 2 in object_locations{traces_filtered{n}(k,4)}(:,2)
        deltax = traces_filtered{n}(k+1,1) -  traces_filtered{n}(k,1);
        surroundings{n}(k, 1) = traces_filtered{n}(k,1);
        surroundings{n}(k, 2) = deltax;
        radius = 115; %10*(deltax); Breite "Einheitszelle Kanal" ueber eine Stenose, dickste Stelle bis dickste Stelle
        surroundings{n}(k, 3) = size (find ((object_locations{traces_filtered{n}(k,4)}(:,3) >  traces_filtered{n}(k,1) - radius) & (object_locations{traces_filtered{n}(k,4)}(:,3) <  traces_filtered{n}(k,1) + radius)),1);
        surroundings{n}(k, 4) = size (find ((object_locations{traces_filtered{n}(k+1,4)}(:,3) >  traces_filtered{n}(k+1,1) - radius) & (object_locations{traces_filtered{n}(k+1,4)}(:,3) <  traces_filtered{n}(k+1,1) + radius)),1);
        surroundings{n}(k, 5) = surroundings{n}(k, 4) - surroundings{n}(k, 3);
    end
    %what happens after trace ends
    if  traces_filtered{n}(k,4) == nr_of_frames
        continue
    else
        if  traces_filtered{n}(end_x,4) == nr_of_frames
            continue
        else
            surroundings{n}(end_x, 1) = traces_filtered{n}(end_x,1);
            surroundings{n}(end_x, 2) = deltax;
            surroundings{n}(end_x, 3) = size (find ((object_locations{traces_filtered{n}(end_x,4)}(:,3) >  traces_filtered{n}(end_x,1) - radius) & (object_locations{traces_filtered{n}(end_x,4)}(:,3) <  traces_filtered{n}(end_x,1) + radius)),1);
            surroundings{n}(end_x, 4) = size (find ((object_locations{traces_filtered{n}(end_x,4)+1}(:,3) >  traces_filtered{n}(end_x,1) + deltax - radius) & (object_locations{traces_filtered{n}(end_x,4)+1}(:,3) <  traces_filtered{n}(end_x,1) + deltax + radius)),1);
            surroundings{n}(end_x, 5) = surroundings{n}(end_x, 4) - surroundings{n}(end_x, 3);
        end
    end
end

outputname_traces_endings = [folder 'traces endings.png'];
figure (5)
clf;
imshow(uint8(bg));
title(sprintf('%d traces (length > %d )',  number_filtered_traces, min_trace_length));
hold on
for n = 1:number_filtered_traces
    plot(traces_filtered{n}(:,1),traces_filtered{n}(:,2),'LineStyle', '-');
    start_x = 1;
    last_x = size(traces_filtered{n},1);
    text(traces_filtered{n}(start_x,1), (traces_filtered{n}(start_x,2)), sprintf('o'), 'Color', 'w', 'FontSize', 10 , 'FontWeight', 'bold');
    text(traces_filtered{n}(last_x,1), (traces_filtered{n}(last_x,2)), sprintf('x'), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
    for k = 1:last_x-1
        if surroundings{n}(k, 5) > 0
            text(traces_filtered{n}(k,1), (traces_filtered{n}(k,2)), sprintf('%d', surroundings{n}(k, 5)), 'Color', 'g', 'FontSize', 10 , 'FontWeight', 'bold');
        end
        if surroundings{n}(k, 5) < 0
            text(traces_filtered{n}(k,1), (traces_filtered{n}(k,2)), sprintf('%d', surroundings{n}(k, 5)), 'Color', 'y', 'FontSize', 10 , 'FontWeight', 'bold');
        end
    end
    text(traces_filtered{n}(last_x,1), (traces_filtered{n}(last_x,2)), sprintf('%d', surroundings{n}(last_x, 5)), 'Color', 'm', 'FontSize', 10 , 'FontWeight', 'bold');  
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_traces_endings);

%%velocity matrix complete
j=0;
for n = 1:number_filtered_traces
    for k = 1:(size( traces_filtered{n}, 1)-1)
        vel_points(j+k, 1) = n; %number of trace
        vel_points(j+k, 2) = 0.5*(traces_filtered{n}(k,1)+traces_filtered{n}(k+1,1)); % x position between x1 and x2
        vel_points(j+k, 3) = 0.5*(traces_filtered{n}(k,2)+traces_filtered{n}(k+1,2)); % y position between y1 and y2
        vel_points(j+k, 4) = (traces_filtered{n}(k+1,1) - traces_filtered{n}(k,1))*meter_per_pixel/dt; % x velocity
        vel_points(j+k, 5) = (traces_filtered{n}(k+1,2) - traces_filtered{n}(k,2))*meter_per_pixel/dt; % y velocity
        if ~isempty(bridges_filtered{n})
            isbridge = find (bridges_filtered{n} (: ,1) == traces_filtered{n}(k+1,1) & bridges_filtered{n} (: ,2) == traces_filtered{n}(k+1,2) & bridges_filtered{n} (: ,3) == traces_filtered{n}(k+1,3), 1);
            if ~isempty(isbridge)
                vel_points(j+k, 4) = 0; % x velocity invalid because of bridging
                vel_points(j+k, 5) = 0; % y velocity invalid because of bridging
            end
        end
        vel_points(j+k, 6) = sqrt(vel_points(j+k, 4)^2+vel_points(j+k, 5)^2); % absolute velocity
        vel_points(j+k, 7) = (traces_filtered{n}(k+1,3) - traces_filtered{n}(k,3));
        vel_points(j+k, 8) = (traces_filtered{n}(k+1,3) - traces_filtered{n}(k,3))/traces_filtered{n}(k,3); % change in area from pos 1 to 2
        vel_points(j+k, 9 ) = surroundings{n}(k, 5);
        if  (vel_points(j+k, 7) < 0) && (vel_points(j+k, 7) < -0.75*median_object_area)
            vel_points(j+k, 10 ) = surroundings{n}(k, 5);
        else
            vel_points(j+k, 10 ) = 0;
        end
    end
    j=j+(size( traces_filtered{n}, 1)-1);
end

outputname_vel = [folder 'traces velocities.png'];
figure (6)
clf;
imshow(uint8(bg));
title(sprintf('%d traces, flowrate %.4f mul/h',  number_filtered_traces, flowrate));
hold on
quiverc(vel_points(:,2),vel_points(:,3),vel_points(:,4),vel_points(:,5));
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
bar = colorbar;
barpos = bar.Position;
export_fig(outputname_vel);

outputname_vel_bar = [folder 'traces velocities bar.png'];
figure (7) %vel plot
clf;
hold on
quiverc(vel_points(:,2),vel_points(:,3),vel_points(:,4),vel_points(:,5));
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
vmax = round(max(vel_points(:, 6))*1000)/10 ; %entspricht 1 im colorbar wegen normierung in quiverc
vmax_third = round(max(vel_points(:, 6))*1000/3)/10 ;
vmax_twothird = round(max(vel_points(:, 6))*2000/3)/10 ;
c = colorbar('Ticks',[0,1/3,2/3,1],'TickLabels',{'0',num2str(vmax_third),num2str(vmax_twothird),num2str(vmax)});
c.Label.String = 'velocity (cm/s)';
barpos(1)=1.03*barpos(1);
set(c,'position',barpos);
export_fig(outputname_vel_bar);
%mean_velocity = mean(nonzeros(first10vel_points(:,6))); unused, avoid stenosis:
%% calculate flowrate ONLY FROM FIRST 10 TRACING STEPS (avoid stenosis)
first10vel_points = zeros(10*number_filtered_traces, 6);
velocity_avgfirst10 = zeros(number_filtered_traces, 1);
j=0;
for n = 1:number_filtered_traces
    for k = 1:(min_trace_length-1)
        first10vel_points(j+k, 1) = n; %number of trace
        first10vel_points(j+k, 2) = 0.5*(traces_filtered{n}(k,1)+traces_filtered{n}(k+1,1)); % x position between x1 and x2
        first10vel_points(j+k, 3) = 0.5*(traces_filtered{n}(k,2)+traces_filtered{n}(k+1,2)); % y position between y1 and y2
        first10vel_points(j+k, 4) = (traces_filtered{n}(k+1,1) - traces_filtered{n}(k,1)); % x velocity
        first10vel_points(j+k, 5) = (traces_filtered{n}(k+1,2) - traces_filtered{n}(k,2)); % y velocity
        if ~isempty(bridges_filtered{n})
            isbridge = find (bridges_filtered{n} (: ,1) == traces_filtered{n}(k+1,1) & bridges_filtered{n} (: ,2) == traces_filtered{n}(k+1,2) & bridges_filtered{n} (: ,3) == traces_filtered{n}(k+1,3), 1);
            if ~isempty(isbridge)
                first10vel_points(j+k, 4) = 0; % x velocity invalid because of bridging
                first10vel_points(j+k, 5) = 0; % y velocity invalid because of bridging
            end
        end
        first10vel_points(j+k, 6) = (meter_per_pixel/dt)*sqrt(first10vel_points(j+k, 4)^2+first10vel_points(j+k, 5)^2); % absolute velocity   
    end
    velocity_avgfirst10(n,1)=(meter_per_pixel/(min_trace_length*dt))*sqrt(((traces_filtered{n}(min_trace_length,1) - traces_filtered{n}(1,1)))^2+((traces_filtered{n}(min_trace_length,2) - traces_filtered{n}(1,2)))^2);
    j=j+min_trace_length;
end

mean_velocity = mean( velocity_avgfirst10);
mean_xdisplacement = mean(nonzeros(first10vel_points(:,4)));
flowrate = canal_crosssection*mean_velocity*3600*10^9; %canal_crosssection(m^2)
disp([num2str(number_filtered_traces) ' filtered traces, flow rate ' num2str(flowrate) ' mul/h']);

%% basic analysis objects size distribution objects (rosettes = >1.5* median_object_area)
objects_inlet = objectlocations_complete(:,5);%column 5 = area in pixels
objects_inlet_foll = objectlocations_complete(:,5);
for  n = size(objectlocations_complete,1):-1:1
    if (objectlocations_complete(n,3) < input_delay_x) || (objectlocations_complete(n,3) > input_delay_x+mean_xdisplacement)
        objects_inlet(n,:) = [];
    end
    if (objectlocations_complete(n,3) < input_delay_x+mean_xdisplacement) || (objectlocations_complete(n,3) > input_delay_x+2*mean_xdisplacement)
        objects_inlet_foll(n,:) = [];
    end
    
end

Xinlet=objects_inlet/median_object_area;
Xinletfoll=objects_inlet_foll/median_object_area;
hist_edges = [0 1.5 2.22 3.33 4.44 5.56 7.78 30];
hist_mitte = hist_edges+0.555;
[Ninlet,edgesinlet,bininlet]  = histcounts(Xinletfoll,hist_edges);
[Ninletfoll,edgesinletfoll,binfoll]  = histcounts(Xinletfoll,hist_edges);

filename_size = sprintf('size_dist.txt');
SizeInletName  = [folder filename_size];
fid = fopen(SizeInletName,'w');
if fid ~= -1
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\t%s\r\n','results rosette tracker', folder );
    fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(number_filtered_traces) );
    fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
    fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(median_object_area) );
    fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
    fprintf(fid,'%s\t%s\r\n','counts: ', num2str(Ninlet) );
    fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(Ninlet/sum(Ninlet)) );
    fprintf(fid,'%s\t%s\r\n','counts: ', num2str(Ninletfoll) );
    fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(Ninletfoll/sum(Ninlet)) );
    fprintf(fid,'%s\r\n', '===...===');
end
fclose(fid);

%%plot area vs x position
outputname_area = [folder 'area vs x position.png'];
figure (8)
clf;
title(sprintf('area vs x position'));
hold on
for n = 1:number_filtered_traces
    plot(traces_filtered{n}(:,1),traces_filtered{n}(:,3),'LineStyle', '-', 'Marker', '.');
    last_x = size(traces_filtered{n},1);
    text(traces_filtered{n}(last_x,1), traces_filtered{n}(last_x,3), sprintf('trace %d', n), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
    bridges_n = size(bridges_filtered{n},1);
    for k = 1:bridges_n
        text(bridges_filtered{n}(k,1),bridges_filtered{n}(k,3), sprintf('%d', k), 'Color', 'g', 'FontSize', 10 , 'FontWeight', 'bold');
    end
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_area);

%%plot normalized area vs x position
outputname_normarea = [folder 'norm area vs x position.png'];
figure (9)
clf;
title(sprintf('normalized area vs x position'));
hold on
for n = 1:number_filtered_traces
    plot(traces_filtered{n}(:,1),(traces_filtered{n}(:,3)/traces_filtered{n}(1,3)),'LineStyle', '-', 'Marker', '.');
    last_x = size(traces_filtered{n},1);
    text(traces_filtered{n}(last_x,1), (traces_filtered{n}(last_x,3)/traces_filtered{n}(1,3)), sprintf('%d', n), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_normarea);

%%plot area change vs x position
changes_along_traces = cell(number_filtered_traces, 1); %identical values to vel_points, but one matrix for each trace
for n = 1:number_filtered_traces
    changes_along_traces{n} = zeros(size( traces_filtered{n}, 1)-1, 6);
    for k = 1:(size( traces_filtered{n}, 1)-1)
        changes_along_traces{n}(k, 1) = 0.5*(traces_filtered{n}(k,1)+traces_filtered{n}(k+1,1));
        changes_along_traces{n}(k, 2) = 0.5*(traces_filtered{n}(k,2)+traces_filtered{n}(k+1,2));
        changes_along_traces{n}(k, 3) = (traces_filtered{n}(k+1,1) - traces_filtered{n}(k,1))/dt;
        changes_along_traces{n}(k, 4) = (traces_filtered{n}(k+1,2) - traces_filtered{n}(k,2))/dt;
        changes_along_traces{n}(k, 5) = sqrt(changes_along_traces{n}(k, 3)^2+changes_along_traces{n}(k, 4)^2);
        changes_along_traces{n}(k, 6) = (traces_filtered{n}(k+1,3) - traces_filtered{n}(k,3)); %/traces_filtered{n}(k,3) %>0 = groesser, <0= kleiner
    end
end

%%identify events - rupture and connect
endings = zeros (number_filtered_traces, 5);
end_rupture_events = 0;
end_connect_events = 0;
event_counter = zeros (number_filtered_traces, 4); %#trace #ruptures #connects ending
connect_counter = 0; %#trace #frame x y #objdiff
rupture_counter = 0;

for n = 1:number_filtered_traces
    event_counter(n,1)=n;
    for k = 1: size( traces_filtered{n}, 1)-1 %trace n, step k in trace and in sourroundings
        current_frame = traces_filtered{n}(k,4); %current_frame in video
        %rupture event?
        if surroundings{n}(k, 5) > 0 && changes_along_traces{n}(k, 6) < 0 && abs(changes_along_traces{n}(k, 6)) > 0.65*median_object_area %%if #objects and area changes along traces
            event_counter(n,2)=event_counter(n,2)+1;
            rupture_counter = rupture_counter +1;
            traces_filtered{n}(k,9)=changes_along_traces{n}(k, 6);
            ruptures( rupture_counter, 1) = n;
            ruptures( rupture_counter, 2) = current_frame;
            ruptures( rupture_counter, 3) = traces_filtered{n}(k,1);
            ruptures( rupture_counter, 4) = traces_filtered{n}(k,2);
            ruptures( rupture_counter, 5) =  surroundings{n}(k, 5);
            ruptures( rupture_counter, 6) =   changes_along_traces{n}(k, 6) ;
        end
        
        %connect event?
        if surroundings{n}(k, 5) < 0 && changes_along_traces{n}(k, 6) > 0 && abs(changes_along_traces{n}(k, 6)) > 0.65*median_object_area %%if #objects and area changes along traces
            event_counter(n,3)=event_counter(n,3)+1;
            connect_counter = connect_counter+1;
            traces_filtered{n}(k,9)=changes_along_traces{n}(k, 6);
            connects( connect_counter, 1) = n;
            connects( connect_counter, 2) = current_frame;
            connects( connect_counter, 3) = traces_filtered{n}(k,1);
            connects( connect_counter, 4) = traces_filtered{n}(k,2);
            connects( connect_counter, 5) = surroundings{n}(k, 5);
            connects( connect_counter, 6) = changes_along_traces{n}(k, 6) ;
        end
    end
    
    %end of trace?
    k =   size( traces_filtered{n}, 1);
    current_frame = traces_filtered{n}(k,4);
    endings (n,1) = n;
    endings (n,2) = current_frame;
    endings (n,3) = traces_filtered{n}(k,1);
    endings (n,4) = traces_filtered{n}(k,2);
    endings (n,5) =  surroundings{n}(k, 5);
    %end of trace rupture event?
    if surroundings{n}(k, 5) > 0 %%if #objects and area changes along traces
        end_rupture_events = end_rupture_events+1;
    end
    
    %end of trace: connect event?
    if surroundings{n}(k, 5) < 0  %if #objects and area changes along traces
        end_connect_events = end_connect_events +1;
    end
end

end_no_events = number_filtered_traces-end_rupture_events- end_connect_events;

outputname_events = [folder 'events.png'];
figure (10)
clf;
imshow(uint8(bg));
title(sprintf('%d traces (length > %d )',  number_filtered_traces, min_trace_length));
hold on
for n = 1:number_filtered_traces
    plot(traces_filtered{n}(:,1),traces_filtered{n}(:,2),'LineStyle', '-', 'Color', 'w');
end
for n = 1:number_filtered_traces
    for k = 1:size(changes_along_traces{n},1)
        if changes_along_traces{n}(k,6)<-2*median_object_area
            text(traces_filtered{n}(k,1), traces_filtered{n}(k,2), sprintf('x'), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
        end
    end
end
for i = 1:size(ruptures, 1)
    text(ruptures(i,3), ruptures(i,4), sprintf('%d', ruptures(i,5)), 'Color', 'g', 'FontSize', 10 , 'FontWeight', 'bold');
end
for i = 1:size(connects, 1)
    text(connects(i,3), connects(i,4), sprintf('%d', connects(i,5)), 'Color', 'y', 'FontSize', 10 , 'FontWeight', 'bold');
end
for i = 1:size(endings , 1)
    text(endings (i,3), endings (i,4), sprintf('%d', endings (i,5)), 'Color', 'm', 'FontSize', 10 , 'FontWeight', 'bold');
end

hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_events);

outputname_darea = [folder 'change in area vs x position.png'];
figure (11)
clf;
title(sprintf('change in area vs x position'));
hold on
for n = 1:number_filtered_traces
    plot( changes_along_traces{n}(:,1), changes_along_traces{n}(:,6),'LineStyle', '-', 'Marker', '.');
    last_x = size(changes_along_traces{n},1);
    text(changes_along_traces{n}( last_x,1),changes_along_traces{n}( last_x,6), sprintf('%d', n), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_darea);

%% categorize events along traces - real rupture and pass by
for n = 1:number_filtered_traces
    for k = 1:(size( traces_filtered{n}, 1))
        traces_filtered{n}(k,10) = 0; %clear column
    end
end

for n = 1:number_filtered_traces
    event_column = traces_filtered{n}(:,9);
    l = find(traces_filtered{n}(:,9));
    if size(l,1)>1
        for m = 1:(size(l,1)-1)
            if l(m)+1 == l(m+1) && ( event_column(l(m+1)) * event_column(l(m)) > 0) %aufeinanderfolgende Indizes und selbes Vorzeichen -> wird zusammengefasst
                event_column(l(m+1))= event_column(l(m+1)) + event_column(l(m));
                event_column(l(m)) = 0;
            end
        end
    end
    ev = find(event_column);
    while ~isempty(ev)
        if traces_filtered{n}(ev(1),9)<0 %rupture
            if size(ev ,1) == 1
                traces_filtered{n}(ev(1),10) = 1; %real rupture
                ev(1)=[];
                continue
            end
            if size(ev ,1) > 1
                if (ev(2)-ev(1)) > 25 || traces_filtered{n}(ev(2),9) < 0
                    traces_filtered{n}(ev(1),10) = 1; %real rupture
                end
                if (ev(2)-ev(1)) <=25 && traces_filtered{n}(ev(2),9) > 0
                    for p = ev(1):ev(2)
                        traces_filtered{n}(p,10) = 2; %connect after rupture
                    end
                    ev(2)=[];
                end
                ev(1)=[];
                continue
            end
        end
        if traces_filtered{n}(ev(1),9)>0 %connect
            if size(ev ,1) == 1
                traces_filtered{n}(ev(1),10) = 3; %connect only
                ev(1)=[];
                continue
            end
            if size(ev ,1) > 1
                if (ev(2)-ev(1)) > 25 || traces_filtered{n}(ev(2),9) > 0
                    traces_filtered{n}(ev(1),10) = 3; %connect only
                end
                if (ev(2)-ev(1)) <=25 && traces_filtered{n}(ev(2),9) < 0
                    for p = ev(1):ev(2)
                        traces_filtered{n}(p,10) = 4; %pass by
                    end
                    ev(2)=[];
                end
                ev(1)=[];
                continue
            end
        end
    end
end

%%Farben Uni Augsburg
uni_gruen = [0 101/255 97/255];
uni_lila = [173/255 0 124/255]; %RGB: R:173 G:0 B:124
uni_gelb = [246/255 168/255 0]; %R:246 G:168 B:0
uni_orange = [ 235/255 105/255 11/255]; %R:235 G:105 B:11
uni_rot = [ 212/255 0 45/255];%R:212 G:0 B:45
uni_hellblau =[0 174/255 207/255];% R:0 G:174 B:207
uni_dunkelblau =[21/255 85/255 129/255];%R:21 G:85 B:129
uni_hellgruen = [72/255 147/255 36/255];%R:72 G:147 B:36
uni_grau = [179/255 179/255 179/255];%R:179 G:179 B:179

%%Plots
%one for each trace
for n = 1:number_filtered_traces
    filename = sprintf('trace %d.png', n);
    outputname_trace = [folder filename];
    figure (12)
    clf
    hold on
    title(sprintf('trace %d.png', n));
    plot(traces_filtered{n}(:,1),traces_filtered{n}(:,3),traces_filtered{n}(:,1),traces_filtered{n}(:,5),traces_filtered{n}(:,1),traces_filtered{n}(:,6), traces_filtered{n}(:,1),100*traces_filtered{n}(:,7), changes_along_traces{n}(:,1), changes_along_traces{n}(:,6));
    legend('Area', 'MajorAxis', 'MinorAxis', 'Eccentricity', 'Change in Area');
    for j = 1:size(changes_along_traces{n},1)
        if traces_filtered{n}(j,9) ~= 0
            text( traces_filtered{n}(j+1,1),traces_filtered{n}(j+1,3), sprintf('%d', traces_filtered{n}(j,9)), 'Color', 'r', 'FontSize', 10 , 'FontWeight', 'bold');
        end
    end
    if size(traces_filtered{n},2)==10
        for j =1:(size(traces_filtered{n},1)-1)
            if traces_filtered{n}(j,10) ~= 0
                text( traces_filtered{n}(j,1),traces_filtered{n}(j,3), sprintf('%d', traces_filtered{n}(j,10)), 'Color', 'g', 'FontSize', 10 , 'FontWeight', 'bold');
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(outputname_trace);
end

%one subplot for each event category
figure (13)
clf;
filename = sprintf('categorized events along traces.png');
outputname_cat_events = [folder filename];
hold on
subplot(5, 1, 1);
imshow(uint8(bg));
title(sprintf('real ruptures'));
for n = 1:30%number_filtered_traces
    realruptures = find(traces_filtered{n}(:,10)==1);
    if ~isempty(realruptures)
        for r = 1:size(realruptures)
            rel_fontsize = round(-0.1*traces_filtered{n}(realruptures(r),9));
            text(traces_filtered{n}(realruptures(r),1), (traces_filtered{n}(realruptures(r),2)), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
    end
end
subplot(5, 1, 2);
imshow(uint8(bg));
title(sprintf('rupture & reconnect'));
for n = 1:30%number_filtered_traces
    reconnects = find(traces_filtered{n}(:,10)==2);
    if ~isempty(reconnects)
        plotstartvalues = 0*reconnects;
        for x = 1:(size(reconnects,1)-1)
            if reconnects(x+1)==reconnects(x)+1
                plotstartvalues(x+1)=1;
            end
        end
        plotstartvalues = find(plotstartvalues==0);
        hold on
        for r = 1:(size(plotstartvalues,1)-1)
            plot(traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(plotstartvalues(r+1)-1),1), (traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(plotstartvalues(r+1)-1),2)),  'Color', 'w');
            text(traces_filtered{n}(reconnects(plotstartvalues(r)),1), (traces_filtered{n}(reconnects(plotstartvalues(r)),2)), sprintf('x'), 'Color', 'r', 'FontSize', 8 , 'FontWeight', 'bold');
            text(traces_filtered{n}(reconnects(plotstartvalues(r+1)-1),1), (traces_filtered{n}(reconnects(plotstartvalues(r+1)-1),2)), sprintf('o'), 'Color', 'y', 'FontSize', 8 , 'FontWeight', 'bold');
        end
        r = size(plotstartvalues,1);
        plot(traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(size(reconnects)),1), (traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(size(reconnects)),2)),  'Color', 'w');
        text(traces_filtered{n}(reconnects(plotstartvalues(r)),1), (traces_filtered{n}(reconnects(plotstartvalues(r)),2)), sprintf('x'), 'Color', 'r', 'FontSize', 8 , 'FontWeight', 'bold');
        text(traces_filtered{n}(reconnects(size(reconnects,1)),1), (traces_filtered{n}(reconnects(size(reconnects,1)),2)), sprintf('o'), 'Color', 'y', 'FontSize', 8 , 'FontWeight', 'bold');
        hold off
    end
end
subplot(5, 1, 3);
imshow(uint8(bg));
title(sprintf('connect only'));
for n = 1:30%number_filtered_traces
    connectsonly = find(traces_filtered{n}(:,10)==3);
    if ~isempty(connectsonly)
        for r = 1:size(connectsonly)
            text(traces_filtered{n}(connectsonly(r),1), (traces_filtered{n}(connectsonly(r),2)), sprintf('o'), 'Color', 'y', 'FontSize', 10 , 'FontWeight', 'bold');
        end
    end
end
subplot(5, 1, 4);
imshow(uint8(bg));
title(sprintf('pass by'));
for n = 1:30%number_filtered_traces
    
    passby = find(traces_filtered{n}(:,10)==4);
    if ~isempty(passby)
        plotstartvalues = 0*passby;
        for x = 1:(size(passby,1)-1)
            if passby(x+1)==passby(x)+1
                plotstartvalues(x+1)=1;
            end
        end
        plotstartvalues = find(plotstartvalues==0);
        hold on
        if size(plotstartvalues,1)>1
            for r = 1:(size(plotstartvalues,1)-1)
                plot(traces_filtered{n}(passby(plotstartvalues(r)):passby(plotstartvalues(r+1)-1),1), (traces_filtered{n}(passby(plotstartvalues(r)):passby(plotstartvalues(r+1)-1),2)),  'Color', 'g');
                text(traces_filtered{n}(passby(plotstartvalues(r)),1), (traces_filtered{n}(passby(plotstartvalues(r)),2)), sprintf('o'), 'Color', 'y', 'FontSize', 8 , 'FontWeight', 'bold');
                text(traces_filtered{n}(passby(plotstartvalues(r+1)-1),1), (traces_filtered{n}(passby(plotstartvalues(r+1)-1),2)), sprintf('x'), 'Color', 'w', 'FontSize', 8 , 'FontWeight', 'bold');
            end
        end
        r = size(plotstartvalues,1);
        plot(traces_filtered{n}(passby(plotstartvalues(r)):passby(size(passby,1)),1), (traces_filtered{n}(passby(plotstartvalues(r)):passby(size(passby,1)),2)),  'Color', 'g');
        text(traces_filtered{n}(passby(plotstartvalues(r)),1), (traces_filtered{n}(passby(plotstartvalues(r)),2)), sprintf('o'), 'Color', 'y', 'FontSize', 8 , 'FontWeight', 'bold');
        text(traces_filtered{n}(passby(size(passby,1)),1), (traces_filtered{n}(passby(size(passby,1)),2)), sprintf('x'), 'Color', 'w', 'FontSize', 8 , 'FontWeight', 'bold');
        hold off
    end
end
subplot(5, 1, 5);
imshow(uint8(bg));
title(sprintf('end of trace'));
for n = 1:30%number_filtered_traces
    endxv = size(traces_filtered{n},1);
    text(traces_filtered{n}(endxv,1), (traces_filtered{n}(endxv,2)), sprintf('x'), 'Color', 'm', 'FontSize', 10 , 'FontWeight', 'bold');
end

hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig(outputname_cat_events);

%grid to plot categorized events along traces per area
x_steps = input_delay_x:10:image_width;
x_steps_center = x_steps+5;
x_steps_center(:,size(x_steps,2))=[];
y_steps = -5:10:image_height;
y_steps_center = y_steps+5;
y_steps_center(:,size(y_steps,2))=[];
grid = zeros(size(y_steps_center,2),size(x_steps_center,2));
figure (14)
clf;
filename = sprintf('categorized events along traces per area.png');
outputname_cat_area = [folder filename];
hold on
counter_realruptures = grid;
sum_realruptures = grid;
means_realruptures = grid;
for n = 1:number_filtered_traces
    realruptures = find(traces_filtered{n}(:,10)==1);
    if ~isempty(realruptures)
        for r = 1:size(realruptures)
            i = find (abs(traces_filtered{n}(realruptures(r),1) - x_steps_center(1, :)) <=5);
            i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
            j = find (abs(traces_filtered{n}(realruptures(r),2) - y_steps_center(1, :)) <=5);
            j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
            counter_realruptures(j,i)=counter_realruptures(j,i)+1;
            sum_realruptures(j,i)=sum_realruptures(j,i)+traces_filtered{n}(realruptures(r),9);
        end
    end
end

subplot(6, 1, 1);
imshow(uint8(bg));
title(sprintf('real ruptures (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_realruptures(j,i)>0
            means_realruptures(j,i) = sum_realruptures(j,i)/counter_realruptures(j,i);
            rel_fontsize = round(0.5*counter_realruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
    end
end
hold off

subplot(6, 1, 2);
imshow(uint8(bg));
title(sprintf('real ruptures (mean values per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_realruptures(j,i)>0
            rel_fontsize = round(-0.1*means_realruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
        end
    end
end
hold off

counter_ruptures = grid;
counter_reconnects = grid;
subplot(6, 1, 3);
imshow(uint8(bg));
title(sprintf('rupture & reconnect (counts per area)'));
hold on
for n = 1:number_filtered_traces
    reconnects = find(traces_filtered{n}(:,10)==2);
    if ~isempty(reconnects)
        plotstartvalues = 0*reconnects;
        for x = 1:(size(reconnects,1)-1)
            if reconnects(x+1)==reconnects(x)+1
                plotstartvalues(x+1)=1;
            end
        end
        plotstartvalues = find(plotstartvalues==0);
        
        for r = 1:(size(plotstartvalues,1)-1)
            plot(traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(plotstartvalues(r+1)-1),1), (traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(plotstartvalues(r+1)-1),2)),  'Color', 'w');
            i = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r)),1) - x_steps_center(1, :)) <=5);
            i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
            j = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r)),2) - y_steps_center(1, :)) <=5);
            j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
            counter_ruptures(j,i)=counter_ruptures(j,i)+1;
            ii = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r+1)-1),1) - x_steps_center(1, :)) <=5);
            ii = ii(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
            jj = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r+1)-1),2) - y_steps_center(1, :)) <=5);
            jj = jj(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
            counter_reconnects(jj,ii)=counter_reconnects(jj,ii)+1;
        end
        r = size(plotstartvalues,1);
        plot(traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(size(reconnects)),1), (traces_filtered{n}(reconnects(plotstartvalues(r)):reconnects(size(reconnects)),2)),  'Color', 'w');
        i = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r)),1) - x_steps_center(1, :)) <=5);
        i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
        j = find (abs(traces_filtered{n}(reconnects(plotstartvalues(r)),2) - y_steps_center(1, :)) <=5);
        j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
        counter_ruptures(j,i)=counter_ruptures(j,i)+1;
        ii = find (abs(traces_filtered{n}(reconnects(size(reconnects,1)),1) - x_steps_center(1, :)) <=5);
        ii = ii(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
        jj = find (abs(traces_filtered{n}(reconnects(size(reconnects,1)),2) - y_steps_center(1, :)) <=5);
        jj = jj(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
        counter_reconnects(jj,ii)=counter_reconnects(jj,ii)+1;
        
    end
end
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_ruptures(j,i)>0
            rel_fontsize = round(0.5*counter_ruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
        if counter_reconnects(j,i)>0
            rel_fontsize = round(0.5*counter_reconnects(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', 'y');
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
    end
end
hold off
counter_connects = grid;
for n = 1:number_filtered_traces
    connectsonly = find(traces_filtered{n}(:,10)==3);
    if ~isempty(connectsonly)
        for r = 1:size(connectsonly)
            i = find (abs(traces_filtered{n}( connectsonly(r),1) - x_steps_center(1, :)) <=5);
            i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
            j = find (abs(traces_filtered{n}( connectsonly(r),2) - y_steps_center(1, :)) <=5);
            j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
            counter_connects(j,i)=counter_connects(j,i)+1;
        end
    end
end
subplot(6, 1, 4);
imshow(uint8(bg));
title(sprintf('connect only (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_connects(j,i)>0
            rel_fontsize = round(0.5*counter_connects(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', 'y');
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
    end
end
hold off

counter_passbystart = grid;
counter_passbyend = grid;
subplot(6, 1, 5);
imshow(uint8(bg));
title(sprintf('pass by (counts per area)'));
hold on
for n = 1:number_filtered_traces
    passby = find(traces_filtered{n}(:,10)==4);
    if ~isempty(passby)
        plotstartvalues = 0*passby;
        for x = 1:(size(passby,1)-1)
            if passby(x+1)==passby(x)+1
                plotstartvalues(x+1)=1;
            end
        end
        plotstartvalues = find(plotstartvalues==0);
        if size(plotstartvalues,1)>1
            for r = 1:(size(plotstartvalues,1)-1)
                plot(traces_filtered{n}(passby(plotstartvalues(r)):passby(plotstartvalues(r+1)-1),1), (traces_filtered{n}(passby(plotstartvalues(r)):passby(plotstartvalues(r+1)-1),2)),  'Color', [0.8 0.8 0.8]);
                i = find (abs(traces_filtered{n}(passby(plotstartvalues(r)),1) - x_steps_center(1, :)) <=5);
                i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
                j = find (abs(traces_filtered{n}(passby(plotstartvalues(r)),2) - y_steps_center(1, :)) <=5);
                j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
                counter_passbystart(j,i)=counter_passbystart(j,i)+1;
                ii = find (abs(traces_filtered{n}(passby(plotstartvalues(r+1)-1),1) - x_steps_center(1, :)) <=5);
                ii = ii(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
                jj = find (abs(traces_filtered{n}(passby(plotstartvalues(r+1)-1),2) - y_steps_center(1, :)) <=5);
                jj = jj(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
                counter_passbyend(jj,ii)=counter_passbyend(jj,ii)+1;
            end
        end
        r = size(plotstartvalues,1);
        plot(traces_filtered{n}(passby(plotstartvalues(r)):passby(size(passby,1)),1), (traces_filtered{n}(passby(plotstartvalues(r)):passby(size(passby,1)),2)),  'Color',  uni_grau);
        i = find (abs(traces_filtered{n}(passby(plotstartvalues(r)),1) - x_steps_center(1, :)) <=5);
        i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
        j = find (abs(traces_filtered{n}(passby(plotstartvalues(r)),2) - y_steps_center(1, :)) <=5);
        j = j(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
        counter_passbystart(j,i)=counter_passbystart(j,i)+1;
        ii = find (abs(traces_filtered{n}(passby(size(passby,1)),1) - x_steps_center(1, :)) <=5);
        ii = ii(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
        jj = find (abs(traces_filtered{n}(passby(size(passby,1)),2) - y_steps_center(1, :)) <=5);
        jj = jj(1,1); %falls y wert in zwei Felder faellt (auf der Kante liegt) wird er oben eingeordnet - untere Kante gehoert zum Feld
        counter_passbyend(jj,ii)= counter_passbyend(jj,ii)+1;
    end
end
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_passbystart(j,i)>0
            rel_fontsize = round(0.5*counter_passbystart(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', 'y');
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
        if counter_passbyend(j,i)>0
            rel_fontsize = round(0.5*counter_passbyend(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_orange);
            %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
        end
    end
end
hold off

counter_endings = grid(1,:);
j=1;
for n = 1:number_filtered_traces
    endxv = size(traces_filtered{n},1);
    i = find (abs(traces_filtered{n}(endxv,1) - x_steps_center(1, :)) <=5);
    i = i(1,1); %falls x wert in zwei Felder faellt (auf der Kante liegt) wird er links eingeordnet - rechte Kante gehoert zum Feld
    counter_endings(j,i)=counter_endings(j,i)+1;
end
subplot(6, 1, 6);
imshow(uint8(bg));
title(sprintf('end of trace'));
hold on
for i = 1:size(grid,2)
    if counter_endings(j,i)>0
        rel_fontsize = round(0.5*counter_endings(j,i));
        plot(x_steps_center(1, i), y_steps_center(1, 5),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_lila);
    end
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig (outputname_cat_area);

figure (15)
clf;
filename = sprintf('categorized events over area.png');
outputname_cat_x = [folder filename];
hold on
subplot(6, 1, 1);
title(sprintf('real ruptures (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_realruptures(j,i)>0
            rel_fontsize = round(0.5*counter_realruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
    end
end
hold off

subplot(6, 1, 2);
title(sprintf('real ruptures (mean values per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_realruptures(j,i)>0
            rel_fontsize = round(-0.1*means_realruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
    end
end
hold off
subplot(6, 1, 3);
title(sprintf('rupture & reconnect (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_ruptures(j,i)>0
            rel_fontsize = round(0.5*counter_ruptures(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
        if counter_reconnects(j,i)>0
            rel_fontsize = round(0.5*counter_reconnects(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
    end
end
hold off
subplot(6, 1, 4);
title(sprintf('connect only (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_connects(j,i)>0
            rel_fontsize = round(0.5*counter_connects(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
    end
end
hold off
subplot(6, 1, 5);
title(sprintf('pass by (counts per area)'));
hold on
for i = 1:size(grid,2)
    for j = 1:size(grid,1)
        if counter_passbystart(j,i)>0
            rel_fontsize = round(0.5*counter_passbystart(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
        if counter_passbyend(j,i)>0
            rel_fontsize = round(0.5*counter_passbyend(j,i));
            plot(x_steps_center(1, i), y_steps_center(1, j),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_orange);
            xlim([0 image_width]);
            ylim([0 image_height]);
        end
    end
end
hold off
subplot(6, 1, 6);
title(sprintf('end of trace'));
hold on
for i = 1:size(grid,2)
    if counter_endings(1,i)>0
        rel_fontsize = round(0.5*counter_endings(1,i));
        plot(x_steps_center(1, i), y_steps_center(1, 5),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_lila);
        xlim([0 image_width]);
        ylim([0 image_height]);
    end
end
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig (outputname_cat_x);
%% legende
figure (16)
clf;
filename = sprintf('legend norm 500.png');
outputname_legend = [folder filename];
subplot(6, 1, 1);
title(sprintf('bullet size corresponding to counts per 500 traces'));
bullets=[5 25 50 100];
y_bullets =[0.75  0.25 -0.25 -0.85];
hold on
for i = 1:size(bullets,2)
    rel_fontsize = round(bullets(i));
    plot((image_width/2)-10,y_bullets(i),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', 'k');
    text((image_width/2)+10,y_bullets(i), num2str(bullets(i)), 'FontSize', 10 , 'FontWeight', 'bold');
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
set(gca,'xtick',[]);
set(gca,'xcolor',[1 1 1]);
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig (outputname_legend);
%grid reduced to one dimension: categorized events along flow direction
figure (17)
clf;
filename = sprintf('categorized events along flow direction.png');
outputname_cat_vs_x = [folder filename];
hold on
subplot(6, 1, 1);
title(sprintf('real ruptures (counts)'));
xcounter_realruptures=sum(counter_realruptures);
hold on
for i = 1:size(grid,2)
    if xcounter_realruptures(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_realruptures(1,i));%statt 0.5 mult mit 500/number_fil rel_fontsize = round(0.5*xcounter_realruptures(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), y_steps_center(1, 1),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
subplot(6, 1, 2);
title(sprintf('real ruptures (mean values area change)'));
xmean_realruptures = sum(means_realruptures);
hold on
for i = 1:size(grid,2)
    if xcounter_realruptures(1,i)>0
        rel_fontsize = round(-0.1*xmean_realruptures(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), y_steps_center(1, 1),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
subplot(6, 1, 3);
title(sprintf('catch up (counts rupture (red) and reconnect (green))'));
xcounter_ruptures = sum(counter_ruptures);
xcounter_reconnects = sum(counter_reconnects);
hold on
for i = 1:size(grid,2)
    if xcounter_ruptures(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_ruptures(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), 0.25,  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_rot);
        end
    end
    if xcounter_reconnects(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_reconnects(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), -0.25,  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
subplot(6, 1, 4);
title(sprintf('connect only (counts)'));
xcounter_connects = sum(counter_connects);
hold on
for i = 1:size(grid,2)
    if xcounter_connects(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_connects(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), y_steps_center(1, 1),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
subplot(6, 1, 5);
title(sprintf('pass by (counts start (green) and end (orange))'));
xcounter_passbystart = sum(counter_passbystart);
xcounter_passbyend = sum(counter_passbyend);
hold on
for i = 1:size(grid,2)
    if xcounter_passbystart(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_passbystart(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), 0.25,  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_gruen);
        end
    end
    if xcounter_passbyend(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*xcounter_passbyend(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), -0.25,  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_orange);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
subplot(6, 1, 6);
title(sprintf('end of trace (counts)'));
hold on
for i = 1:size(grid,2)
    if counter_endings(1,i)>0
        rel_fontsize = round(500/number_filtered_traces*counter_endings(1,i));
        if rel_fontsize>0
            plot(x_steps_center(1, i), y_steps_center(1, 1),  'Marker', '.', 'MarkerSize', rel_fontsize, 'Color', uni_lila);
        end
    end
end
xlim([0 image_width]);
ylim([-1 1]);
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1]);
hold off
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'color','white')
export_fig (outputname_cat_vs_x);

%%txt export
%one for each trace
for n = 1:number_filtered_traces
    filename_trace = sprintf('rosette_tracker_trace %d.txt', n);
    TracesDataFileName = [folder filename_trace];
    fid = fopen(TracesDataFileName,'w');
    if fid ~= -1
        
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s%s\r\n','results rosette tracker', folder , video_filename);
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(number_filtered_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n','#trace', '#frame','x Pos (pixel)','y Pos (pixel)','Area','MajorAxis','MinorAxis','Eccentricity', 'dArea', 'dx', '#sur objects', '#sur objects fr+1', 'diff #objects', 'event', 'cat event');
        fprintf(fid,'%s\r\n', '===...===');
        
        for k = 1:(size( traces_filtered{n}, 1)-1)
            fprintf(fid,'%s\t%s\t\t%s\t\t%s\t\t%s\t%s\t\t%s\t\t%s\t\t%s\t%s\t%s\t\t%s\t\t\t%s\t\t%s\t%s\r\n',num2str(n), num2str(traces_filtered{n}(k,4)), num2str(traces_filtered{n}(k,1)),num2str(traces_filtered{n}(k,2)),num2str(traces_filtered{n}(k,3)),num2str(traces_filtered{n}(k,5)),num2str(traces_filtered{n}(k,6)),num2str(traces_filtered{n}(k,7)),num2str(changes_along_traces{n}(k,6)),num2str(surroundings{n}(k, 2)),num2str(surroundings{n}(k, 3)),num2str(surroundings{n}(k, 4)),num2str(surroundings{n}(k, 5)),num2str(traces_filtered{n}(k,9)),num2str(traces_filtered{n}(k,10)));
        end
        k = k+1;
        fprintf(fid,'%s\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\r\n',num2str(n), num2str(traces_filtered{n}(k,4)), num2str(traces_filtered{n}(k,1)),num2str(traces_filtered{n}(k,2)),num2str(traces_filtered{n}(k,3)),[],num2str(surroundings{n}(k, 2)),num2str(surroundings{n}(k, 3)),num2str(surroundings{n}(k, 4)),num2str(surroundings{n}(k, 5)));
        
    end
    
    fclose(fid);
end

%event_counter_grid.txt
EventCounterFileName = [folder 'event_counter_grid.txt'];
fid = fopen(EventCounterFileName,'w');
if fid ~= -1
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\t%s%s\r\n','results rosette tracker', folder , video_filename);
    fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(number_filtered_traces) );
    fprintf(fid,'%s\t%s\r\n','mean x displacement (px): ', num2str(mean_xdisplacement) );
    fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
    fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(median_object_area) );
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','grid x coordinates');
    fprintf(fid,'%s\r\n','grid y coordinates');
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','realruptures (counts per area)');
    fprintf(fid,'%s\r\n','realruptures (counts vs x)');
    fprintf(fid,'%s\r\n','realruptures (mean values area change per area)');
    fprintf(fid,'%s\r\n','realruptures (mean values area change vs x)');
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','catch up ruptures (counts per area)');
    fprintf(fid,'%s\r\n','catch up ruptures (counts vs x)');
    fprintf(fid,'%s\r\n','catch up reconnects (counts per area)');
    fprintf(fid,'%s\r\n','catch up reconnects (counts vs x)');
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','connect only (counts per area)');
    fprintf(fid,'%s\r\n','connect only (counts vs x)');
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','pass by start (counts per area)');
    fprintf(fid,'%s\r\n','pass by start (counts vs x)');
    fprintf(fid,'%s\r\n','pass by end (counts per area)');
    fprintf(fid,'%s\r\n','pass by end (counts vs x)');
    fprintf(fid,'%s\r\n', '===...===');
    fprintf(fid,'%s\r\n','end of trace (counts vs x)');
    fprintf(fid,'%s\r\n', '===...===');
end
fclose(fid);

dlmwrite(EventCounterFileName,x_steps_center,'-append','delimiter','\t','newline','pc');
dlmwrite(EventCounterFileName,y_steps_center,'-append','delimiter','\t','newline','pc','roffset',1);
dlmwrite(EventCounterFileName,counter_realruptures,'-append','delimiter','\t','newline','pc', 'roffset',2);
dlmwrite(EventCounterFileName,xcounter_realruptures,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,means_realruptures,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,xmean_realruptures,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_ruptures,'-append','delimiter','\t','newline','pc', 'roffset',2);
dlmwrite(EventCounterFileName,xcounter_ruptures,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_reconnects,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,xcounter_reconnects,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_connects,'-append','delimiter','\t','newline','pc', 'roffset',2);
dlmwrite(EventCounterFileName,xcounter_connects,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_passbystart,'-append','delimiter','\t','newline','pc', 'roffset',2);
dlmwrite(EventCounterFileName,xcounter_passbystart ,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_passbyend,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,xcounter_passbyend,'-append','delimiter','\t','newline','pc', 'roffset',1);
dlmwrite(EventCounterFileName,counter_endings ,'-append','delimiter','\t','newline','pc', 'roffset',2);

rosette_tracker_part_2 = 1;
