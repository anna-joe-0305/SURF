directory = "J:\Stockholm 2019\";
folders = [ "BGA\BGA 5mum Standard 50mulh_20190327_135904\"; ...];

for z=1:size(folders,1)
    folder =directory + folders(z);
    t = datetime('now');
    disp(['## Start Folder ' num2str(z) ': ' datestr(t)]);
    disp ([ folder ]);
    filename_eventgrid = sprintf('%sevent_counter_grid.txt',folder);
    fid = fopen(filename_eventgrid);
    firstline = fgetl(fid);
    secline = fgetl(fid);
    thirdline = fgetl(fid);
    meanx = fgetl(fid);
    flowrate = fgetl(fid);%prelimininary flowrate from textreader
    flowrate = erase(flowrate,"flow rate (mul/h): 	");
    flowrate = str2num(flowrate);
    single_cell_size  = fgetl(fid);%, [5 1 5 2]
    single_cell_size = erase(single_cell_size,"single cell size (px): 	");
    single_cell_size = str2num(single_cell_size);
    fclose(fid);
    
    %%read tracing files
    file_dir = char(folder + 'rosette_tracker_trace*.txt');
    txt_trace_list = dir (file_dir );
    nr_of_traces = numel(txt_trace_list); % Count, Anzahl Bilder im Ordner
    rosettes = zeros(nr_of_traces,3);
    rosettes_sorted = cell(6, 1);
    trace_imported = cell(nr_of_traces,1);
    for t = 1:nr_of_traces
        filename_trace_t = sprintf('%srosette_tracker_trace %d.txt',folder, t);
        trace_imported{t}  = dlmread(filename_trace_t,'', 7, 0);
        rosettes(t,1) =  trace_imported{t}(1,5)/single_cell_size; %area
        rosettes(t,2) = trace_imported{t}(end,3);%end x
        rosettes(t,3) = trace_imported{t}(end,4);%end y
    end
    
    %%categorize rosettes by size in classes 1-6
    X=rosettes(:,1);
    hist_edges = [1.11 2.22 3.33 4.44 5.56 7.78 30];
    hist_mitte = hist_edges+0.555;
    [N,edges,bin]  = histcounts(X,hist_edges);
    
    background_filename = sprintf('%sBG.png', folder);
    bg = imread(background_filename);
    
    figure(1)
    clf;
    filename = char(folder +'ends cat rosette size.png');
    hold on
    for bin_nr = 1:6
        subplot(6, 1,bin_nr);
        imshow(uint8(bg));
        bin_indices = find(bin ==bin_nr);% bin index = row in rosettes = nr of trace
        if ~isempty (bin_indices)
            for r = 1:size(bin_indices)
                text(rosettes(bin_indices(r),2), rosettes(bin_indices(r),3), sprintf('x'), 'Color', 'm', 'FontSize', 10 , 'FontWeight', 'normal');
                rosettes_sorted{bin_nr}(r,1)=bin_indices(r);
                rosettes_sorted{bin_nr}(r,2:4)=rosettes(bin_indices(r),:);
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white');
    export_fig (filename);
    
    %%calculate flowrate in each class
    dt = 0.0005; %delay frames (s)
    canal_depth = 8*10^-6;
    meter_per_pixel = 0.55*10^-6;
    canal_crosssection = 90*meter_per_pixel*canal_depth;
    flowrates = zeros (nr_of_traces,6);
    for t = 1:nr_of_traces
        flowrates(t,1)=t;
        flowrates(t,2)=bin(t);
        if size(trace_imported{t,1},1) >9 %at least 10 tracing steps
            step = 2;
            dx = 0;
            while dx < 200 && size(trace_imported{t,1},1) > step %find first 200 pixels of trace
                dx = trace_imported{t,1}(step,3)-trace_imported{t,1}(1,3); %letzter x-wert - erster x wert
                step = step+1;
            end
            flowrates(t,3)=dx;
            dy =trace_imported{t,1}(step,4)-trace_imported{t,1}(1,4); %letzter y wert - erster y wert
            flowrates(t,4)=dy;
            velocity = (meter_per_pixel/((step-1)*dt))*sqrt(dx^2+dy^2);
            flowrates(t,5)=velocity;
            flowrate = canal_crosssection*velocity*3600*10^9;
            flowrates(t,6)=flowrate;
        end
    end
    counter_flow = zeros(1,6);
    medians_flow = zeros(1,6);
    for bin_nr = 2:6
        bin_indices = find(bin ==bin_nr);% bin index = row in rosettes = nr of trace
        if ~isempty (bin_indices)
            for r = 1:size(bin_indices)
                if flowrates(bin_indices(r),6) ~= 0
                    counter_flow (1,bin_nr)=   counter_flow (1,bin_nr)+1;
                    medians_flow (1,bin_nr) =  medians_flow (1,bin_nr)+flowrates(bin_indices(r),6);
                end
            end
        end
    end
    total_median_flow = sum(medians_flow)/sum(counter_flow);
    median_2and3= (medians_flow(2)+medians_flow(3) )/(counter_flow(2)+counter_flow(3));
    medians_flow = medians_flow./counter_flow;
    
    
    SizeFileName = folder + 'flowrates_by_size.txt';
    fid = fopen(SizeFileName,'w');
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results rosette tracker', folder );
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(nr_of_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate_old (mul/h): ', num2str(flowrate_old) );
        fprintf(fid,'%s\t%s\r\n','flow rate_new (mul/h): ', num2str(total_median_flow) );
        fprintf(fid,'%s\t%s\r\n','flow rate_class 2 and 3 (mul/h): ', num2str(median_2and3) );
        fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(single_cell_size) );
        fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
        fprintf(fid,'%s\t%s\r\n','counts: ', num2str(N) );
        fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(N/sum(N)) );
        fprintf(fid,'%s\t%s\r\n','counts (length>9): ', num2str(counter_flow) );
        fprintf(fid,'%s\t%s\r\n','medians flow: ', num2str(medians_flow) );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', 'trace   class   dx   dy      velo        flowrate');
    end
    fclose(fid);
    dlmwrite(SizeFileName,flowrates,'-append','delimiter','\t','newline','pc','roffset',1);
    
    %%export rosettes sorted by size
    
    rosettes_sorted_columns =NaN (max(N), 6);
    for bin_nr = 1:6
        if N(bin_nr) > 0
            rosettes_sorted_columns(1:N(bin_nr),bin_nr) =      rosettes_sorted{bin_nr}(:,2);
        end
    end
    SizeFileName = folder + 'rosettes_by_size.txt';
    fid = fopen(SizeFileName,'w');
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results rosette tracker', folder );
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(nr_of_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
        fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(single_cell_size) );
        fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
        fprintf(fid,'%s\t%s\r\n','counts: ', num2str(N) );
        fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(N/sum(N)) );
        fprintf(fid,'%s\r\n', '===...===');
    end
    fclose(fid);
    dlmwrite(SizeFileName,rosettes_sorted_columns,'-append','delimiter','\t','newline','pc','roffset',1);
    
    
    
    %%count events for each rosette size
    x_steps = 0:10:size(bg,2);
    x_steps_center = x_steps+5;
    x_steps_center(:,size(x_steps,2))=[];
    grid = zeros(6,size(x_steps_center,2));%6 lines for 6 bins
    counter_endings = grid;
    for bin_nr = 1:6
        if ~isempty (rosettes_sorted{bin_nr})
            X=rosettes_sorted{bin_nr}(:,3);
            [counter_endings_binnr,edgesends,binends]  = histcounts(X,x_steps);
            counter_endings(bin_nr,:)=counter_endings_binnr;
        end
    end
    
    %%real ruptures
    counter_realruptures = grid;
    sum_area_realruptures = grid;
    mean_area_realruptures = grid;
    
    
    figure(2)
    clf;
    filenamerupt = char(folder +'real rupt cat rosette size.png');
    hold on
    for bin_nr = 1:6
        subplot(6, 1,bin_nr);
        imshow(uint8(bg));
        for r = 1:size(rosettes_sorted{bin_nr}, 1)
            rr_ind = find (trace_imported{rosettes_sorted{bin_nr}(r,1)}(:,15)==1);%column 15 = cat event, "1" = real rupture
            if ~isempty(rr_ind)
                for s=1:size(rr_ind)
                    x_koor = trace_imported{rosettes_sorted{bin_nr}(r,1)}(rr_ind(s),3);
                    y_koor = trace_imported{rosettes_sorted{bin_nr}(r,1)}(rr_ind(s),4);
                    lost = - trace_imported{rosettes_sorted{bin_nr}(r,1)}(rr_ind(s),14);%column 14 = area lost in pixels
                    text(x_koor, y_koor, sprintf('x'), 'Color', 'r', 'FontSize', lost/10 , 'FontWeight', 'normal');
                    counter_realruptures (bin_nr, floor(x_koor/10))=counter_realruptures(bin_nr, floor(x_koor/10))+1;
                    sum_area_realruptures(bin_nr, floor(x_koor/10))=sum_area_realruptures(bin_nr, floor(x_koor/10))+lost;
                end
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenamerupt);
    
    counter_undamaged = grid;
    relfreq_undamaged = grid;
    figure(3)
    clf;
    filenameund = char(folder +'undamaged cat rosette size.png');
    for bin_nr = 1:6
        
        for r = 1:size(rosettes_sorted{bin_nr}, 1)
            rr_ind = find (trace_imported{rosettes_sorted{bin_nr}(r,1)}(:,15)==1,1);%first rupture, column 15 = cat event, "1" = real rupture
            if isempty(rr_ind)
                rr_ind = size(trace_imported{rosettes_sorted{bin_nr}(r,1)});
            end
            x_koor = trace_imported{rosettes_sorted{bin_nr}(r,1)}(rr_ind(1),3);
            und_vector = grid(1,:);
            und_vector(:,1:floor(x_koor/10))=1;%value 1 as long as undamaged, then 0
            counter_undamaged(bin_nr,:)=counter_undamaged(bin_nr,:)+und_vector;
        end
    end
    for bin_nr = 1:6
        relfreq_undamaged(bin_nr,:)=counter_undamaged(bin_nr,:)/N(bin_nr);
    end
    plot(x_steps_center, relfreq_undamaged);
    legend('bin 1','bin 2','bin 3','bin 4','bin 5','bin 6')
    %plot matrix zeilenweise undamaged
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenameund);
    
    %%plots event counts
    figure(4)
    clf;
    filenameevents = char(folder +'event counts cat rosette size.png');
    subplot(3, 1, 1);
    title(sprintf('ends of traces (counts)'));
    hold on
    for i = 1:size(grid,2)
        for j = 1:size(grid,1)
            if counter_endings(j,i)>0
                rel_fontsize = round(counter_endings(j,i));
                plot(x_steps_center(1, i), - j,  'Marker', 'o', 'MarkerSize', rel_fontsize, 'Color', 'm');
                %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
            end
        end
    end
    xlim([0 size(bg,2)]);
    ylim([-7 0]);
    set(gca,'ytick',[-6 -5 -4 -3 -2 -1]);
    set(gca, 'yticklabels',{'bin 6', 'bin 5', 'bin 4', 'bin 3', 'bin 2', 'bin 1'});
    hold off
    
    subplot(3, 1, 2);
    title(sprintf('real ruptures (counts)'));
    hold on
    for i = 1:size(grid,2)
        for j = 1:size(grid,1)
            if counter_realruptures(j,i)>0
                mean_area_realruptures(j,i) = sum_area_realruptures(j,i)/counter_realruptures(j,i);
                rel_fontsize = round(counter_realruptures(j,i));
                plot(x_steps_center(1, i), - j,  'Marker', 'o', 'MarkerSize', rel_fontsize, 'Color', 'r');
                %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
            end
        end
    end
    xlim([0 size(bg,2)]);
    ylim([-7 0]);
    set(gca,'ytick',[-6 -5 -4 -3 -2 -1]);
    set(gca, 'yticklabels',{'bin 6', 'bin 5', 'bin 4', 'bin 3', 'bin 2', 'bin 1'});
    hold off
    
    subplot(3, 1, 3);
    title(sprintf('real ruptures (mean lost area)'));
    hold on
    for i = 1:size(grid,2)
        for j = 1:size(grid,1)
            if mean_area_realruptures(j,i) >0
                rel_fontsize = round(0.01*mean_area_realruptures(j,i) );
                plot(x_steps_center(1, i), - j,  'Marker', 'o', 'MarkerSize', rel_fontsize, 'Color', 'r');
                %text(x_steps_center(1, i), y_steps_center(1, j), sprintf('x'), 'Color', 'r', 'FontSize', rel_fontsize , 'FontWeight', 'bold');
            end
        end
    end
    xlim([0 size(bg,2)]);
    ylim([-7 0]);
    set(gca,'ytick',[-6 -5 -4 -3 -2 -1]);
    set(gca, 'yticklabels',{'bin 6', 'bin 5', 'bin 4', 'bin 3', 'bin 2', 'bin 1'});
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenameevents);
    
    EventFileName = folder + 'event_grid_by_size.txt';
    fid = fopen(EventFileName,'w');
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results rosette tracker', folder );
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(nr_of_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
        fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(single_cell_size) );
        fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
        fprintf(fid,'%s\t%s\r\n','counts: ', num2str(N) );
        fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(N/sum(N)) );
        fprintf(fid,'%s\r\n', '===...===');
    end
    fclose(fid);
    dlmwrite(EventFileName,counter_endings,'-append','delimiter','\t','newline','pc','roffset',1);
    dlmwrite(EventFileName,counter_realruptures,'-append','delimiter','\t','newline','pc','roffset',1);
    dlmwrite(EventFileName,mean_area_realruptures,'-append','delimiter','\t','newline','pc','roffset',1);
    dlmwrite(EventFileName,relfreq_undamaged,'-append','delimiter','\t','newline','pc','roffset',1);
    close all
    
end