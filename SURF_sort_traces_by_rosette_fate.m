directory = "C:\Users\joettean\Desktop\Paper AJ 2021\SURF Code\";
folders = [ "SURF_example_video\";
    ];
stenosis_end_x=[425;];

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
    flowrate = fgetl(fid);
    flowrate = erase(flowrate,"flow rate (mul/h): 	");
    flowrate = str2num(flowrate);
    single_cell_size  = fgetl(fid);%, [5 1 5 2]
    single_cell_size = erase(single_cell_size,"single cell size (px): 	");
    single_cell_size = str2num(single_cell_size);
    fclose(fid);
    
    %%read tracing files
    file_dir = char(folder + 'SURF_trace*.txt');
    txt_trace_list = dir (file_dir );
    nr_of_traces = numel(txt_trace_list); % Count, Anzahl trace*.txt im Ordner
    rosettes = zeros(nr_of_traces,3); %matrix mit allen getracten rosetten, spalten: rosettenflaeche, end x, end y
    rosettes_sorted_bin = cell(6, 1);%cell mit matrix pro bin, spalten: #trace, ros_area, end x, end y, 1 or 0 (1=survived, 0=trace too short)
    rosettes_sorted_undamaged =cell(6,1);
    trace_imported = cell(nr_of_traces,1);%cell mit kompletten traces aus den txt files
    counter_fate_traces = zeros(6,4);%Zeilen: Bins. Spalten: #all traces #survival #damaged #undamaged
    
    for t = 1:nr_of_traces
        filename_trace_t = sprintf('%sSURF_trace %d.txt',folder, t);
        trace_imported{t}  = dlmread(filename_trace_t,'', 7, 0);
        rosettes(t,1) =  trace_imported{t}(1,5)/single_cell_size; %area
        rosettes(t,2) = trace_imported{t}(end,3);%end x
        rosettes(t,3) = trace_imported{t}(end,4);%end y
    end
    
    %%categorize rosettes by size
    X=rosettes(:,1);
    %hist_edges=[200, 400, 600, 800, 1000, 1400, 5000];
    hist_edges = [1.11 2.22 3.33 4.44 5.56 7.78 30];
    hist_mitte = hist_edges+0.555;
    [N,edges,bin]  = histcounts(X,hist_edges);
    
    background_filename = sprintf('%sBG.png', folder);
    bg = imread(background_filename);
    
    figure(1)
    clf;
    filename = char(folder +'ends rosette size.png');
    hold on
    for bin_nr = 1:6
        subplot(6, 1,bin_nr);
        imshow(uint8(bg));
        bin_indices = find(bin ==bin_nr);% bin index = row in rosettes = nr of trace
        if ~isempty (bin_indices)
            for r = 1:size(bin_indices)
                text(rosettes(bin_indices(r),2), rosettes(bin_indices(r),3), sprintf('x'), 'Color', 'm', 'FontSize', 10 , 'FontWeight', 'normal');
                rosettes_sorted_bin{bin_nr}(r,1)=bin_indices(r);
                rosettes_sorted_bin{bin_nr}(r,2:4)=rosettes(bin_indices(r),:);
                counter_fate_traces (bin_nr,1)= counter_fate_traces (bin_nr,1)+1; %count all
                rosettes_sorted_bin{bin_nr}(r,5)=0;
                if rosettes(bin_indices(r),2) > stenosis_end_x(z) %trace longer than stenosis?
                    rosettes_sorted_bin{bin_nr}(r,5)=1;
                    counter_fate_traces (bin_nr,2)= counter_fate_traces (bin_nr,2)+1; %count survival
                end
                %rosettes_sorted (#trace area/single_cell_size end_x end_y)
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white');
    export_fig (filename);
    
    %%count events for each rosette size
    x_steps = 0:10:size(bg,2);
    x_steps_center = x_steps+5;
    x_steps_center(:,size(x_steps,2))=[];
    grid = zeros(6,size(x_steps_center,2));%6 Zeilen fuer 6 bins
    counter_endings = grid;
    for bin_nr = 1:6
        if ~isempty (rosettes_sorted_bin{bin_nr})
            X=rosettes_sorted_bin{bin_nr}(:,3);
            [counter_endings_binnr,edgesends,binends]  = histcounts(X,x_steps);
            counter_endings(bin_nr,:)=counter_endings_binnr;
        end
    end
    
    
    %%real ruptures
    counter_realruptures_dam = grid;
    sum_area_realruptures_dam = grid;
    counter_realruptures= grid;
    sum_area_realruptures= grid;
    mean_area_realruptures = grid;
    counter_undamaged = grid;
    relfreq_undamaged = grid;
    damaged_area = zeros(6,1);
    
    figure(2)
    clf;
    filenamerupt = char(folder +'real rupt cat rosette size.png');
    hold on
    for bin_nr = 1:6
        subplot(6, 1,bin_nr);
        imshow(uint8(bg));
        for r = 1:size(rosettes_sorted_bin{bin_nr}, 1)
            endstenosisline = find(  trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(:,3)> stenosis_end_x(z), 1 );
            rr_ind = find (trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(1:endstenosisline,15)==1);%column 15 = cat event, "1" = real rupture
            if ~isempty(rr_ind)
                for s=1:size(rr_ind)
                    x_koor = trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),3);
                    y_koor = trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),4);
                    lost = - trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),14);%column 14 = area lost in pixels
                    text(x_koor, y_koor, sprintf('x'), 'Color', 'r', 'FontSize', lost/10 , 'FontWeight', 'normal');
                    counter_realruptures (bin_nr, floor(x_koor/10))=counter_realruptures_dam(bin_nr, floor(x_koor/10))+1;
                    sum_area_realruptures(bin_nr, floor(x_koor/10))=sum_area_realruptures_dam(bin_nr, floor(x_koor/10))+lost;
                end
            end
            if rosettes_sorted_bin{bin_nr}(r,5)==1
                rr_ind = find (trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(1:endstenosisline,15)==1);%column 15 = cat event, "1" = real rupture
                if ~isempty(rr_ind)
                    for s=1:size(rr_ind)
                        x_koor = trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),3);
                        y_koor = trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),4);
                        lost = - trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(s),14);%column 14 = area lost in pixels
                        text(x_koor, y_koor, sprintf('x'), 'Color', 'r', 'FontSize', lost/10 , 'FontWeight', 'normal');
                        counter_realruptures_dam (bin_nr, floor(x_koor/10))=counter_realruptures_dam(bin_nr, floor(x_koor/10))+1;
                        sum_area_realruptures_dam(bin_nr, floor(x_koor/10))=sum_area_realruptures_dam(bin_nr, floor(x_koor/10))+lost;
                        %
                        %lost_area =  trace_imported{rosettes_sorted{bin_nr}(r,1)}(1,5)-trace_imported{rosettes_sorted{bin_nr}(r,1)}(endstenosisline,5);
                        if x_koor < stenosis_end_x(z)
                            damaged_area(bin_nr) =  damaged_area(bin_nr) + lost;
                        end
                    end
                    counter_fate_traces (bin_nr,3)=counter_fate_traces (bin_nr,3)+1;%counter damaged
                    
                end
                if isempty(rr_ind)
                    rr_ind = size(trace_imported{rosettes_sorted_bin{bin_nr}(r,1)});
                    counter_fate_traces (bin_nr,4)=counter_fate_traces (bin_nr,4)+1;%counter undamaged
                end
                x_koor = trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(rr_ind(1),3);
                und_vector = grid(1,:);
                und_vector(:,1:floor(x_koor/10))=1;%value 1 as long as undamaged, then 0
                counter_undamaged(bin_nr,:)=counter_undamaged(bin_nr,:)+und_vector;
                %Liste der Nummern von undamaged traces je bin:
                rosettes_sorted_undamaged{bin_nr}=[rosettes_sorted_undamaged{bin_nr}; rosettes_sorted_bin{bin_nr}(r,1)];
                %max(counter_fate_traces(:,4)){bin_nr}(:,size(rosettes_undamaged{bin_nr},2)+1)=trace_imported{rosettes_sorted_bin{bin_nr}(r,1)}(:,10);
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenamerupt);
    %cells for undamaged: v and epsilon
    rosettes_undamaged_epsilon = cell(6, max(counter_fate_traces(:,4)));%pro bin eine cell fuer undamaged, v und epsilon
    longest_trace = 10;
    
    for bin_nr = 1:6
        for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)%wie viele rosetten in diesem bin undamaged
            rosettes_undamaged_epsilon{bin_nr, rose_nr}=trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,3);%3: x koordinate
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,2)=trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,10);%10: dx, 6:major, 7:minor
            
            major = trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,6);
            minor = trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,7);
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,3)=major;
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,4)=minor;
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,5)=(major-minor)./(major+minor);
            %rupture and reconnect aussortieren
            for zeile = 2:size(rosettes_undamaged_epsilon{bin_nr, rose_nr},1)
                if  trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile,15) ==2 && trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile-1,15) ==2
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,2)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,3)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,4)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,5)=0;
                end
            end
            if  size(rosettes_undamaged_epsilon{bin_nr, rose_nr},1)>longest_trace
                longest_trace = size(rosettes_undamaged_epsilon{bin_nr, rose_nr},1);
            end
        end
    end
    velocities = cell(6,1);
    epsilons = cell(6,1);
    velocities_list = cell(6,1);
    velocities_sorted = cell(6,2);%(cell(:,1)=velocities_list_sorted, cell(:,2)=without duplicates
    epsilons_list = cell(6,1);
    epsilons_sorted =cell(6,2);
    
    for bin_nr = 1:6
        relfreq_undamaged(bin_nr,:)=counter_undamaged(bin_nr,:)/N(bin_nr);
    end
    
    lost_cells =sum_area_realruptures./counter_realruptures/single_cell_size ;
    lost_cells_dam =sum_area_realruptures_dam./counter_realruptures_dam/single_cell_size ;
    
    figure(3), clf;
    for bin_nr = 1:6
        subplot(6, 1,bin_nr);
        bar(x_steps_center, lost_cells_dam(bin_nr,:));
        text(stenosis_end_x(z), 0, sprintf('I'), 'Color', 'r', 'FontSize',20 , 'FontWeight', 'normal');
        legend_entry = sprintf('bin  %d',bin_nr);
        legend(legend_entry)
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    
    cells_lost_per_damaged_trace =damaged_area./counter_fate_traces(:,3)/single_cell_size;
    FateFileName = folder + 'damaged_area.txt';
    fid = fopen(FateFileName,'w');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(nr_of_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
        fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(single_cell_size) );
        fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
        fprintf(fid,'%s\t%s\r\n','counts: ', num2str(N) );
        fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(N/sum(N)) );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', 'Zeilen: Bins. Spalten: #all traces #survival #damaged #undamaged');
    end
    fclose(fid);
    dlmwrite(FateFileName,counter_fate_traces,'-append','delimiter','\t','newline','pc','roffset',1);
    fid = fopen(FateFileName,'a');%a=append
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','cells lost per damaged trace' );
    end
    fclose(fid);
    dlmwrite(FateFileName,cells_lost_per_damaged_trace,'-append','delimiter','\t','newline','pc','roffset',1);
    
    figure(4), clf;
    filenameund = char(folder +'undamaged cat rosette size.png');
    plot(x_steps_center, relfreq_undamaged);
    legend('bin 1','bin 2','bin 3','bin 4','bin 5','bin 6')
    %plot matrix zeilenweise undamaged
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenameund);
    
    
    
    
    
    cells_lost_per_damaged_trace =damaged_area./counter_fate_traces(:,3)/single_cell_size;
    FateFileName = folder + 'epsilon_fit_results.txt';
    fid = fopen(FateFileName,'w');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\t%s\r\n','total traces: ', num2str(nr_of_traces) );
        fprintf(fid,'%s\t%s\r\n','flow rate (mul/h): ', num2str(flowrate) );
        fprintf(fid,'%s\t%s\r\n','single cell size (px): ', num2str(single_cell_size) );
        fprintf(fid,'%s\t%s\r\n','edges of bins for rosette size (norm): ', num2str(hist_edges) );
        fprintf(fid,'%s\t%s\r\n','counts: ', num2str(N) );
        fprintf(fid,'%s\t%s\r\n','relative frequencies (counts/counts total): ', num2str(N/sum(N)) );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', 'Zeilen: Bins. Spalten: #all traces #survival #damaged #undamaged');
    end
    fclose(fid);
    dlmwrite(FateFileName,counter_fate_traces,'-append','delimiter','\t','newline','pc','roffset',1);
    fid = fopen(FateFileName,'a');%a=append
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','cells lost per damaged trace' );
    end
    fclose(fid);
    dlmwrite(FateFileName,cells_lost_per_damaged_trace,'-append','delimiter','\t','newline','pc','roffset',1);
  
    
end