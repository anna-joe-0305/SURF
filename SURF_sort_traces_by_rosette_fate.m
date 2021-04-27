directory = "J:\Stockholm 2019\";
folders = [ "BGA\BGA 5mum Standard 50mulh_20190327_135904\"; ...];
stenosis_end_x=[765; ...];

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
    file_dir = char(folder + 'rosette_tracker_trace*.txt');
    txt_trace_list = dir (file_dir );
    nr_of_traces = numel(txt_trace_list); % Count, Anzahl trace*.txt im Ordner
    rosettes = zeros(nr_of_traces,3); %matrix mit allen getracten rosetten, spalten: rosettenflaeche, end x, end y
    rosettes_sorted_bin = cell(6, 1);%cell mit matrix pro bin, spalten: #trace, ros_area, end x, end y, 1 or 0 (1=survived, 0=trace too short)
    rosettes_sorted_undamaged =cell(6,1);
    trace_imported = cell(nr_of_traces,1);%cell mit kompletten traces aus den txt files
    counter_fate_traces = zeros(6,4);%Zeilen: Bins. Spalten: #all traces #survival #damaged #undamaged
    
    for t = 1:nr_of_traces
        filename_trace_t = sprintf('%srosette_tracker_trace %d.txt',folder, t);
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
    
    %     %all_lost_cells_dam =sum_area_realruptures_dam/single_cell_size ;
     end_stenosis = floor(stenosis_end_x(z)/10);
     cells_lost_per_damaged_trace =damaged_area./counter_fate_traces(:,3)/single_cell_size;
     FateFileName = folder + 'damaged_area.txt';
    fid = fopen(FateFileName,'w');%w=write=overwrite existing file
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
    
    gauss = @(x, xdata)x(1)+x(4)/(x(3)*sqrt(pi/2))*exp( -2 * ((xdata - x(2))/x(3)) .^ 2) ;
    laplace =@(x, xdata)exp(-abs(xdata-x(1))/x(2))/(2*x(2))+x(3);
    x_initvel = [10, 500, 100 , 8000];
    %mu = x(1), s = x(2), ~ FWHM/2,35
    figure(5), clf;
    filenamevelo = char(folder +'velocities.png');
    for bin_nr = 2:6
        subplot(2, 5,bin_nr-1);
        for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
            plot(rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,1),rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,2), '.');
            hold on;
        end
        %text(stenosis_end_x(z), 0, sprintf('I'), 'Color', 'r', 'FontSize',20 , 'FontWeight', 'normal');
        legend_entry = sprintf('bin  %d',bin_nr);
        legend(legend_entry)
        axis([0 1500 0 120])
        title('single rosettes');
        xlabel('x (pixel)');
        ylabel('velocity(pixel/0.5ms)');
    end
    for bin_nr = 2:6
        subplot(2, 5,bin_nr+4);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                %letzte Zeile weglassen, v nicht moeglich
                %x berechnen =zwischenwert zwischen Bildern
                x_vel = rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,1)+0.5*rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,2);
                vel_list_here = [x_vel rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,2)];
                for zeile = size(vel_list_here,1):-1:1
                    if vel_list_here(zeile, 2)<=0 || vel_list_here(zeile, 1)>1100
                        vel_list_here(zeile, :)=[];
                    end
                end
                velocities_list{bin_nr}= [ velocities_list{bin_nr}; vel_list_here];
            end
            %plot(velocities_list{bin_nr,1}(:,1),velocities_list{bin_nr,1}(:,2),'.'); hold on;
            velocities_sorted{bin_nr,1} = sortrows(velocities_list{bin_nr,1});
            velocities_sorted{bin_nr,2}(:,1)=unique(velocities_sorted{bin_nr,1}(:,1), 'rows');
            for k=1:size(velocities_sorted{bin_nr,2},1)
                m = find (velocities_sorted{bin_nr,1}(:,1)==velocities_sorted{bin_nr,2}(k,1));
                velocities_sorted{bin_nr,2}(k,2)=mean(velocities_sorted{bin_nr,1}(m(1):m(end),2));
            end
            plot(velocities_sorted{bin_nr,2}(:,1),velocities_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
            hold on;
            vel_smooth=smoothdata(velocities_sorted{bin_nr,2}(:,2));
            velocities_sorted{bin_nr,2}(:,3)=vel_smooth;
            plot(velocities_sorted{bin_nr,2}(:,1), vel_smooth, '.-','Color', uni_lila);
            hold on;
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            axis([0 1500 0 120])
            title('smoothed curve');
            xlabel('x (pixel)');
            ylabel('velocity(pixel/0.5ms)');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenamevelo);
    
    counter_eps = zeros(6,1);
    figure(6), clf;
    filenameepsi = char(folder +'epsilons.png');
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                plot(rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,1),rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,5), '.');
                hold on;
                counter_eps(bin_nr,1)= counter_eps(bin_nr,1)+1;
            end
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            axis([0 1500 0 1])
            title('single rosettes');
            xlabel('x (pixel)');
            ylabel('strain');
        end
    end
    for bin_nr = 2:6
        subplot(2,5,bin_nr+4);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                eps_list_here = [rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,1) rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,5)];
                for zeile = size(eps_list_here,1):-1:1
                    if eps_list_here(zeile, 2)<=0 || eps_list_here(zeile, 1)>1100
                        eps_list_here(zeile, :)=[];
                    end
                end
                epsilons_list{bin_nr}= [ epsilons_list{bin_nr}; eps_list_here];
            end
            %plot(epsilons_list{bin_nr,1}(:,1),epsilons_list{bin_nr,1}(:,2),'.'); hold on;
            epsilons_sorted{bin_nr,1} = sortrows(epsilons_list{bin_nr,1});
            epsilons_sorted{bin_nr,2}(:,1)=unique(epsilons_sorted{bin_nr,1}(:,1), 'rows');
            for k=1:size(epsilons_sorted{bin_nr,2},1)
                m = find (epsilons_sorted{bin_nr,1}(:,1)==epsilons_sorted{bin_nr,2}(k,1));
                epsilons_sorted{bin_nr,2}(k,2)=mean(epsilons_sorted{bin_nr,1}(m(1):m(end),2));
            end
            plot(epsilons_sorted{bin_nr,2}(:,1),epsilons_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
            hold on;
            eps_smooth=smoothdata(epsilons_sorted{bin_nr,2}(:,2));
            epsilons_sorted{bin_nr,2}(:,3)=eps_smooth;
            plot(epsilons_sorted{bin_nr,2}(:,1), eps_smooth, '.-','Color', uni_lila);
            hold on;
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            axis([0 1500 0 1])
            title('smoothed curve');
            xlabel('x (pixel)');
            ylabel('strain');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenameepsi);
    %%stenosenlaenge 80 Pixel
    funcexpA = @(x,xdata)x(1)*exp(x(2)*xdata)+x(3);
    funcexpxc =@(x,xdata)exp(x(2)*(xdata-x(1)))+x(3);
    x_initdefvsx=[4*10^-6, 0.02, 0.2];
    x_initrelaxvsx=[3*10^8, -0.02, 0.2];
   
    %fit relax von ende luecke bis stenosis_end+50 und spiegel fuer deform
    xevenly=0:5:1500;
    xvel = zeros(6,4);
    xdefvst = zeros(6,3);
    xrelaxvst = zeros(6,3);
    velok=zeros(6,1);
    figure(7), clf;
    filenameepsfits = char(folder +'deformrelaxfits.png');
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            plot(velocities_sorted{bin_nr,2}(:,1),velocities_sorted{bin_nr,2}(:,3) , 'x-','Color', uni_grau);
            hold on;
            xvel(bin_nr, :) = lsqcurvefit(gauss,x_initvel,velocities_sorted{bin_nr,2}(:,1),velocities_sorted{bin_nr,2}(:,3));
            plot(xevenly,gauss(xvel(bin_nr,:),xevenly),'.', 'Color', uni_lila);
             if max(gauss(xvel(bin_nr,:),xevenly))>4* median(gauss(xvel(bin_nr,:),xevenly))
                 velok(bin_nr, 1) = 1;
                 legend_entry = 'velocity ok';
             else
                 legend_entry ='velocity fail';
             end
            axis([0 1500 0 100])
            legend(legend_entry);
            title( num2str(bin_nr));
            xlabel('x (pixel)');
            ylabel('velocity(pixel/0.5ms)');
        end
    end
    for bin_nr = 2:6
        subplot(2,5,bin_nr+4);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            plot(epsilons_sorted{bin_nr,2}(:,1),epsilons_sorted{bin_nr,2}(:,3) , 'x-','Color', uni_grau);
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            axis([0 1500 0 1])
            title( num2str(bin_nr));
            xlabel('x (pixel)');
            ylabel('strain');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenameepsfits);
    
    for bin_nr = 2:6
        if velok(bin_nr) && ~isempty(rosettes_sorted_undamaged{bin_nr})
            xsteps = epsilons_sorted{bin_nr,2}(:,1);
            tsteps = 0*xsteps;
            for i = 2:size(xsteps)
                deltat =( xsteps(i)-xsteps(i-1))/gauss(xvel(bin_nr,:),(xsteps(i-1)+0.5*( xsteps(i)-xsteps(i-1))));
                tsteps(i)=tsteps(i-1)+deltat;
            end
            epsilons_sorted{bin_nr,2}(:,4)=0.5*tsteps;
        end
    end

    figure(8), clf;
    filenameepsfits = char(folder +'deformrelaxfitstime.png');
    
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if velok(bin_nr) &&~isempty(rosettes_sorted_undamaged{bin_nr})
            plot(epsilons_sorted{bin_nr,2}(:,4),epsilons_sorted{bin_nr,2}(:,3) , '.','Color', uni_grau);
            hold on;
            [max_eps, max_idx] = max((epsilons_sorted{bin_nr,2}(:,3)));
            gapstart_idx = find(diff(epsilons_sorted{bin_nr,2}(:,1))>50);%index anfang luecke
            if size(gapstart_idx,1)>1
                [closest, closesttomax] = min(abs(gapstart_idx(:)-max_idx));
                gapstart_idx = gapstart_idx(closesttomax);
            end
            %DEFORMATION FIT
            med_def= median(epsilons_sorted{bin_nr,2}(1:max_idx, 3));
            if ~isempty(gapstart_idx)
                [max_eps, max_idx] = max((epsilons_sorted{bin_nr,2}(gapstart_idx-10:gapstart_idx,3)));
                max_idx=max_idx+gapstart_idx-10-1;
                Mdef_idx = max_idx;
            else
            [Mdef, Mdef_idx] = min(abs(0.9*max_eps-epsilons_sorted{bin_nr,2}(max_idx-10:max_idx,3)));
            Mdef_idx=Mdef_idx+max_idx-10-1;
            end
            Mdef =epsilons_sorted{bin_nr,2}(Mdef_idx,3); 
            xlinks = epsilons_sorted{bin_nr,2}(1:Mdef_idx, 4);
            x_initdefvst=[epsilons_sorted{bin_nr,2}(Mdef_idx, 4),2, med_def];
            k = find(abs(med_def-funcexpxc(x_initdefvst,xlinks))<0.05*med_def);
            Bdef_idx=k(end);
            Bdef =epsilons_sorted{bin_nr,2}(Bdef_idx,3);
            xdatadef = epsilons_sorted{bin_nr,2}(Bdef_idx:Mdef_idx, 4);%vs time
            ydatadef = epsilons_sorted{bin_nr,2}(Bdef_idx:Mdef_idx, 3);
            lb = [0.5*x_initdefvst(1),0.1,med_def];
            ub = [1.5*x_initdefvst(1),20,med_def];
            xdefvst(bin_nr, :) = lsqcurvefit(funcexpxc,x_initdefvst, xdatadef,ydatadef, lb, ub);
            plot( epsilons_sorted{bin_nr,2}(1:Mdef_idx, 4),funcexpxc(xdefvst(bin_nr, :),epsilons_sorted{bin_nr,2}(1:Mdef_idx, 4)),'-', 'Color', uni_gruen);
          
            %RELAXATION FIT
            med_relax= median(epsilons_sorted{bin_nr,2}(max_idx:end, 3));
            if ~isempty(gapstart_idx)
                gapend_idx = gapstart_idx+1;
                [max_eps, max_idx] = max((epsilons_sorted{bin_nr,2}(gapend_idx:gapend_idx+10,3)));
                max_idx = max_idx+gapend_idx-1;
                Mrelax_idx = max_idx;
            else
            [Mrelax, Mrelax_idx] = min(abs(0.9*max_eps-epsilons_sorted{bin_nr,2}(max_idx:max_idx+20,3)));
            Mrelax_idx = Mrelax_idx+max_idx-1;
            end
            Mrelax =epsilons_sorted{bin_nr,2}(Mrelax_idx,3);
            xrechts = epsilons_sorted{bin_nr,2}(Mrelax_idx:end, 4);
            x_initrelaxvst=[epsilons_sorted{bin_nr,2}(Mrelax_idx, 4), -2, med_relax];%[3*10^16, -1, 0.15];
            k = find(abs(med_relax-funcexpxc(x_initrelaxvst,xrechts))<0.05*med_relax);
            Brelax_idx=k(1)+Mrelax_idx-1;
            Brelax =epsilons_sorted{bin_nr,2}(Brelax_idx,3);
            xdatarelax = epsilons_sorted{bin_nr,2}(Mrelax_idx:Brelax_idx, 4);
            ydatarelax = epsilons_sorted{bin_nr,2}(Mrelax_idx:Brelax_idx, 3);  
            lb = [0.5*x_initdefvst(1),-20,med_relax];
            ub = [1.5*x_initdefvst(1),-0.1,med_relax];
            xrelaxvst(bin_nr, :) = lsqcurvefit(funcexpxc,x_initrelaxvst,xdatarelax,ydatarelax, lb, ub);
            plot(epsilons_sorted{bin_nr,2}(Mrelax_idx:end, 4),funcexpxc(xrelaxvst(bin_nr, :),epsilons_sorted{bin_nr,2}(Mrelax_idx:end, 4)),'-', 'Color', uni_lila);
            

            plot(xlinks,funcexpxc(x_initdefvst,xlinks),'--k');
            plot(xrechts,funcexpxc(x_initrelaxvst,xrechts),'--k');
            plot(epsilons_sorted{bin_nr,2}(Mdef_idx,4),epsilons_sorted{bin_nr,2}(Mdef_idx,3) , 'xk');
            plot(epsilons_sorted{bin_nr,2}(Bdef_idx,4),epsilons_sorted{bin_nr,2}(Bdef_idx,3) ,  'xk');
            plot(epsilons_sorted{bin_nr,2}(Mrelax_idx,4),epsilons_sorted{bin_nr,2}(Mrelax_idx,3) ,  'xk');
            plot(epsilons_sorted{bin_nr,2}(Brelax_idx,4),epsilons_sorted{bin_nr,2}(Brelax_idx,3) , 'xk');
            
            legend_entry = sprintf('smoothed curve');
            legend_def =  ['$\epsilon_0=' , num2str(xdefvst(bin_nr, 3)) ,'; \tau_{def}=',num2str(1/xdefvst(bin_nr, 2)),'$'];
            legend_relax =  ['$\epsilon_\infty=' , num2str(xrelaxvst(bin_nr, 3)) ,'; \tau_{relax}=',num2str(1/xrelaxvst(bin_nr, 2)),'$'];
            legend({legend_entry, legend_def , legend_relax}, 'Interpreter', 'Latex','Location', 'northoutside');
            axis([round(0.6*epsilons_sorted{bin_nr,2}(max_idx, 4)) round(1.4*epsilons_sorted{bin_nr,2}(max_idx, 4)) 0 1])
            %axis([0 20 0 1]);
            titlebin =sprintf('class  %d',bin_nr);
            title(titlebin);
            xlabel('t (ms)');
            ylabel('strain');
        end
    end
    subplot(2,5,7);
        xdefvst(xdefvst==0) = NaN;
    xrelaxvst(xrelaxvst==0) = NaN;
    bins = 2:6;
    bar(bins,counter_eps(2:6,1),'FaceColor', uni_gruen);
    titlebin =sprintf('Messung  %d',z);
    title(titlebin);
    xlabel('rosette class');
    ylabel('counts');
    subplot(2,5,8); 
    plot(bins, xdefvst(2:6,3), '--x', 'Color', uni_gruen), hold on;
    plot(bins, xrelaxvst(2:6,3), '--x', 'Color', uni_lila)
    legend_def =  ['$\epsilon_0$'];
    legend_relax =  ['$\epsilon_\infty$'];
    legend({legend_def , legend_relax}, 'Interpreter', 'Latex','Location', 'northoutside');
    axis([1.5 6.5 0 0.5])
    xlabel('rosette class');
    ylabel('strain');
    subplot(2,5,9);
    plot(bins,1./xdefvst(2:6, 2), '--x', 'Color', uni_gruen), hold on;
    plot(bins, -1./xrelaxvst(2:6,2), '--x', 'Color', uni_lila);
    legend_def =  ['$\tau_{def}$'];
    legend_relax =  ['$\tau_{relax}$'];
    legend({legend_def , legend_relax}, 'Interpreter', 'Latex','Location', 'northoutside' );
    axis([1.5 6.5 0 2])
    xlabel('rosette class');
    ylabel('time constant');
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white');
 annotation('textbox', [0, 0.5, 0, 0], 'string', folder);
    export_fig (filenameepsfits);
    
        fitresults = zeros(6, 5);
        fitresults(:,1)=counter_eps(:,1);
        fitresults(:,2)=xdefvst(:,3);%eps 0
        fitresults(:,3)=xrelaxvst(:,3);%eps inf
        fitresults(:,4)=1./xdefvst(:, 2);%taudef
        fitresults(:,5)=-1./xrelaxvst(:,2);%taurelax
    end_stenosis = floor(stenosis_end_x(z)/10);
    cells_lost_per_damaged_trace =damaged_area./counter_fate_traces(:,3)/single_cell_size;
    FateFileName = folder + 'epsilon_fit_results.txt';
    fid = fopen(FateFileName,'w');%w=write=overwrite existing file
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
    fid = fopen(FateFileName,'a');%a=append
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','epsilon fit results' );
        fprintf(fid,'%s\r\n', 'Zeilen: Bins. Spalten: #counts #eps0 #epsinf #taudef #taurelax');
    end
    fclose(fid);
    dlmwrite(FateFileName,fitresults,'-append','delimiter','\t','newline','pc','roffset',1);
   fid = fopen(FateFileName,'a');%a=append
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','Geschwindigkeit' );
        fprintf(fid,'%s\r\n', 'Zeilen: Bins. Spalte: Basisgeschwindigkeit aus GaussFit');
    end
    fclose(fid);
    dlmwrite(FateFileName,   xvel(:,1),'-append','delimiter','\t','newline','pc','roffset',1);
   
end