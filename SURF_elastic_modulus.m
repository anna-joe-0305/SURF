directory = "C:\Users\joettean\Desktop\Paper AJ 2021\SURF Code\";
folders = [ "SURF_example_video\";
    ];
stenosis_end_x = [425;
    ];
%%adjust length and width of stenosis
EL =108;%Standard Elongationslaenge 108 Pixel
d1 =88;%Breite Kanal in Pixel
d2 = 14;%Breite Stenose in Pixel

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

for z=1:size(folders,1)
    close all
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
    single_cell_size  = fgetl(fid);
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
    hist_edges = [1.11 2.22 3.33 4.44 5.56 7.78 30];
    hist_mitte =[1.86, 2.775, 3.885, 5, 6.67, 9];
    Apixel=hist_mitte*single_cell_size;
    rpixel=round(sqrt(Apixel/pi));
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
    % export_fig (filename);
    
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
            end
        end
    end
    hold off
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig (filenamerupt);
    
    %%cell rosettes undamaged: v and epsilon
    rosettes_undamaged_epsilon = cell(6, max(counter_fate_traces(:,4)));%pro bin eine cell fuer undamaged, v und epsilon
    longest_trace = 10;
    for bin_nr = 1:6
        for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)%wie viele rosetten in diesem bin undamaged
            rosettes_undamaged_epsilon{bin_nr, rose_nr}=trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,3);%3: x koordinate
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,2)=trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,10);%10: dx, 6:major, 7:minor
            %area? fuer r0 zeile auswaehlen die am naechsten am stenosenstart ist.
            startstenosisline = find(trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,3)> stenosis_end_x(z)-2*EL-75-100, 1 );
            area = trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(startstenosisline,5);
            r0=sqrt(area/pi);
            major = trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,6);
            minor = trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(:,7);
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,3)=major;
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,4)=minor;
            rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,5)=(major-2*r0)./(2*r0);
            %rupture and reconnect aussortieren
            for zeile = 2:size(rosettes_undamaged_epsilon{bin_nr, rose_nr},1)
                if  trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile,15) ==2 && trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile-1,15) ==2
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,2)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,3)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,4)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,5)=0;
                elseif trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile,15) ==4 && trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile-1,15) ==4
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,2)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,3)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,4)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile,5)=0;
                elseif trace_imported{rosettes_sorted_undamaged{bin_nr}(rose_nr,1)}(zeile,15) ==3
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile:end,2)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile:end,3)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile:end,4)=0;
                    rosettes_undamaged_epsilon{bin_nr, rose_nr}(zeile:end,5)=0;
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
    minor_list = cell(6,1);
    minor_sorted =cell(6,2);
    major_list = cell(6,1);
    major_sorted =cell(6,2);
    
    for bin_nr = 1:6
        relfreq_undamaged(bin_nr,:)=counter_undamaged(bin_nr,:)/N(bin_nr);
    end
    
    lost_cells =sum_area_realruptures./counter_realruptures/single_cell_size ;
    lost_cells_dam =sum_area_realruptures_dam./counter_realruptures_dam/single_cell_size ;
    
    gauss = @(x, xdata)x(1)+x(4)/(x(3)*sqrt(pi/2))*exp( -2 * ((xdata - x(2))/x(3)) .^ 2) ;
    x_initvel = [5, 300, 100 , 8000];
    figure(4)
    clf;
    filenamevelo = char(folder +'velocitiesneu.png');
    for bin_nr = 2:6
        subplot(2, 5,bin_nr-1);
        for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
            plot(rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,1),rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,2), '.');
            hold on;
        end
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
            velocities_sorted{bin_nr,1} = sortrows(velocities_list{bin_nr,1});
            velocities_sorted{bin_nr,2}(:,1) = [50:10:1090].'+5;
            for k=1:size(velocities_sorted{bin_nr,2},1)
                m = find ((velocities_sorted{bin_nr,1}(:,1)>=velocities_sorted{bin_nr,2}(k,1)-5));
                l=find(velocities_sorted{bin_nr,1}(:,1)<=velocities_sorted{bin_nr,2}(k,1)+4);
                g=intersect(m,l);
                if ~isempty(g)
                    velocities_sorted{bin_nr,2}(k,2)=mean(velocities_sorted{bin_nr,1}(g(1):g(end),2));
                else
                    velocities_sorted{bin_nr,2}(k,2)=NaN;
                end
            end
            plot(velocities_sorted{bin_nr,2}(:,1),velocities_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
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
    figure(5)
    clf;
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
            axis([0 1500 0 4])
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
            
            epsilons_sorted{bin_nr,1} = sortrows(epsilons_list{bin_nr,1});
            epsilons_sorted{bin_nr,2}(:,1) = [50:10:1090].'+5;
            for k=1:size(epsilons_sorted{bin_nr,2},1)
                m = find ((epsilons_sorted{bin_nr,1}(:,1)>=epsilons_sorted{bin_nr,2}(k,1)-5));
                l=find(epsilons_sorted{bin_nr,1}(:,1)<=epsilons_sorted{bin_nr,2}(k,1)+4);
                g=intersect(m,l);
                if ~isempty(g)
                    epsilons_sorted{bin_nr,2}(k,2)=mean(epsilons_sorted{bin_nr,1}(g(1):g(end),2));
                else
                    epsilons_sorted{bin_nr,2}(k,2)=NaN;
                end
            end
            plot(epsilons_sorted{bin_nr,2}(:,1),epsilons_sorted{bin_nr,2}(:,2),'x', 'Color', uni_grau), hold on;
            for zeile = size(epsilons_sorted{bin_nr,2},1):-1:1
                if epsilons_sorted{bin_nr,2}(zeile, 2)<0.1
                    epsilons_sorted{bin_nr,2}(zeile, 2)=NaN;
                end
            end
            plot(epsilons_sorted{bin_nr,2}(:,1),epsilons_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
            hold on;
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            axis([0 1500 0 4])
            title('smoothed curve');
            xlabel('x (pixel)');
            ylabel('strain');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenameepsi);
    
    funcexpminor =@(x,xdata)-exp(x(2)*(xdata-x(1)))+x(3);
    funcexpmajor =@(x,xdata)exp(x(2)*(xdata-x(1)))+x(3);
    x_initdefvsx=[4*10^-6, 0.02, 0.2];
    x_initrelaxvsx=[3*10^8, -0.02, 0.2];
    xevenly=0:5:1000;
    fitpara_vel = zeros(6,4);
    xvelexp = zeros(6,3);
    xvaluesvel = cell(6,1);
    fitpara_eps = zeros(6,3);
    fitpara_minor = zeros(6,3);
    fitpara_major = zeros(6,3);
    xvaluesdef=cell(6,1);
    xrelaxvsx = zeros(6,3);
    xvaluesrelax=cell(6,1);
    velok=zeros(6,1);
    velmax_x=zeros(6,1);
    %r=[9, 11.7, 13.8, 15.6,18.1, 27]*0.5*10^-6;%m
    v_results=zeros(105,30);
    eps_results=zeros(105, 40);
    figure(6)
    clf;
    filenamevfits = char(folder +'velocityvsx.png');
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            velocities_forfit =velocities_sorted{bin_nr,2};
            for zeile = size(velocities_forfit,1):-1:1
                if isnan(velocities_forfit(zeile, 2))
                    velocities_forfit(zeile, :)=[];
                end
            end
            plot(velocities_forfit(:,1),velocities_forfit(:,2) , 'x-','Color', uni_grau);
            hold on;
            [m,mind]=min(abs(stenosis_end_x(z)-(velocities_forfit(:,1))));
            [maxv,maxind]=max(velocities_forfit(1:mind,2));
            x_initvel(2)=velocities_forfit(maxind,1);
            x_initvel(1)=mean(velocities_forfit(1:5,2));
            lb = [x_initvel(1),0.9*x_initvel(2),0,10^3];
            ub = [x_initvel(1),1.1*x_initvel(2),300,10^5];
            fitpara_vel(bin_nr, :) = lsqcurvefit(gauss,x_initvel,velocities_forfit(1:mind,1),velocities_forfit(1:mind,2), lb, ub);
            plot(xevenly,gauss(fitpara_vel(bin_nr,:),xevenly),'.', 'Color', uni_lila);
            ydatavel = smoothdata(velocities_forfit(:, 2),'movmean',5);
            plot(velocities_forfit(:, 1),ydatavel ,'.', 'Color', uni_dunkelblau);
            
            [maxv,maxind]=max(gauss(fitpara_vel(bin_nr,:),xevenly));
            if maxv>(0.5*d1/d2)* median(gauss(fitpara_vel(bin_nr,:),xevenly))
                velok(bin_nr, 1) = 1;
                legend_entry = 'velocity ok';
                velmax_x(bin_nr, 1) = xevenly(maxind);
            else
                legend_entry ='velocity fail';
            end
            axis([0 1500 0 150])
            legend(legend_entry);
            title( num2str(bin_nr));
            xlabel('x (pixel)');
            ylabel('velocity(pixel/0.5ms)');
            v_results(:, (bin_nr-2)*6+1)=velocities_sorted{bin_nr,2}(:,1);
            v_results(:, (bin_nr-2)*6+2)=velocities_sorted{bin_nr,2}(:,2);
            v_results(:, (bin_nr-2)*6+3)=gauss(fitpara_vel(bin_nr,:),velocities_sorted{bin_nr,2}(:,1));
            [maxv, maxvind]= max(v_results(:, (bin_nr-2)*3+3));
            v_results(:, (bin_nr-2)*6+4)=v_results(:, (bin_nr-2)*3+1)-v_results(maxvind, (bin_nr-2)*3+1);
            v_baseline=mode(v_results(1:size(velocities_sorted{bin_nr,2}(:,1),1), (bin_nr-2)*3+3));
            v_results(:, (bin_nr-2)*6+5)=v_results(:, (bin_nr-2)*3+2)/v_baseline;
            v_results(:, (bin_nr-2)*6+6)=v_results(:, (bin_nr-2)*3+3)/v_baseline;
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenamevfits);
    %%major and minor fits
    counter_eps = zeros(6,1);
    figure(7)
    clf;
    filenameminor = char(folder +'minor.png');
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                plot(rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,1),rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,4), '.');
                hold on;
                counter_eps(bin_nr,1)= counter_eps(bin_nr,1)+1;
            end
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            title('single rosettes');
            xlabel('x (pixel)');
            ylabel('minor');
        end
    end
    for bin_nr = 2:6
        subplot(2,5,bin_nr+4);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                minor_list_here = [rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,1) rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,4)];
                for zeile = size(minor_list_here,1):-1:1
                    if minor_list_here(zeile, 2)<=0 || minor_list_here(zeile, 1)>1100
                        minor_list_here(zeile, :)=[];
                    end
                end
                minor_list{bin_nr}= [ minor_list{bin_nr}; minor_list_here];
            end
            
            minor_sorted{bin_nr,1} = sortrows(minor_list{bin_nr,1});
            minor_sorted{bin_nr,2}(:,1) = [50:10:1090].'+5;
            for k=1:size(minor_sorted{bin_nr,2},1)
                m = find ((minor_sorted{bin_nr,1}(:,1)>=minor_sorted{bin_nr,2}(k,1)-5));
                l=find(minor_sorted{bin_nr,1}(:,1)<=minor_sorted{bin_nr,2}(k,1)+4);
                g=intersect(m,l);
                if ~isempty(g)
                    minor_sorted{bin_nr,2}(k,2)=mean(minor_sorted{bin_nr,1}(g(1):g(end),2));
                else
                    minor_sorted{bin_nr,2}(k,2)=NaN;
                end
            end
            plot(minor_sorted{bin_nr,2}(:,1),minor_sorted{bin_nr,2}(:,2),'x', 'Color', uni_grau), hold on;
            plot(minor_sorted{bin_nr,2}(:,1),minor_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
            hold on;
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            title('smoothed curve');
            xlabel('x (pixel)');
            ylabel('minor');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenameminor);
    
    counter_eps = zeros(6,1);
    figure(8)
    clf;
    filenamemajor = char(folder +'major.png');
    for bin_nr = 2:6
        subplot(2,5,bin_nr-1);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                plot(rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,1),rosettes_undamaged_epsilon{bin_nr, rose_nr}(:,3), '.');%3=major, 4=minor, 5=epsilon
                hold on;
                counter_eps(bin_nr,1)= counter_eps(bin_nr,1)+1;
            end
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            title('single rosettes');
            xlabel('x (pixel)');
            ylabel('major');
        end
    end
    for bin_nr = 2:6
        subplot(2,5,bin_nr+4);
        if ~isempty(rosettes_sorted_undamaged{bin_nr})
            for rose_nr = 1:size(rosettes_sorted_undamaged{bin_nr},1)
                major_list_here = [rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,1) rosettes_undamaged_epsilon{bin_nr, rose_nr}(1:end-1,3)];
                for zeile = size(major_list_here,1):-1:1
                    if major_list_here(zeile, 2)<=0 || major_list_here(zeile, 1)>1100
                        major_list_here(zeile, :)=[];
                    end
                end
                major_list{bin_nr}= [ major_list{bin_nr}; major_list_here];
            end
            
            major_sorted{bin_nr,1} = sortrows(major_list{bin_nr,1});
            major_sorted{bin_nr,2}(:,1) = [50:10:1090].'+5;
            for k=1:size(major_sorted{bin_nr,2},1)
                m = find ((major_sorted{bin_nr,1}(:,1)>=major_sorted{bin_nr,2}(k,1)-5));
                l=find(major_sorted{bin_nr,1}(:,1)<=major_sorted{bin_nr,2}(k,1)+4);
                g=intersect(m,l);
                if ~isempty(g)
                    major_sorted{bin_nr,2}(k,2)=mean(major_sorted{bin_nr,1}(g(1):g(end),2));
                else
                    major_sorted{bin_nr,2}(k,2)=NaN;
                end
            end
            plot(major_sorted{bin_nr,2}(:,1),major_sorted{bin_nr,2}(:,2),'x', 'Color', uni_grau), hold on;
            plot(major_sorted{bin_nr,2}(:,1),major_sorted{bin_nr,2}(:,2),'x', 'Color', uni_gruen)
            hold on;
            legend_entry = sprintf('bin  %d',bin_nr);
            legend(legend_entry)
            title('smoothed curve');
            xlabel('x (pixel)');
            ylabel('major');
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(filenamemajor);
    
    %%
    figure(9)
    clf;
    filenameellipsefits = char(folder +'ellipsevsx.png');
    for bin_nr = 3:6
        subplot(2,5,bin_nr-1);
        if velok(bin_nr) &&~isempty(rosettes_sorted_undamaged{bin_nr})
            plot(minor_sorted{bin_nr,2}(:,1),minor_sorted{bin_nr,2}(:,2) , '.','Color', uni_grau);
            hold on;
            %DEFORMATION FIT MINOR AXIS
            minoraxisvsx = minor_sorted{bin_nr,2};
            for zeile = size(minoraxisvsx,1):-1:1
                if isnan(minoraxisvsx(zeile, 2))
                    minoraxisvsx(zeile, :)=[];
                end
            end
            [m,mind]=min(abs(stenosis_end_x(z)-(minoraxisvsx(:,1))));%values until stenosis end
            [mindef, Mdef_idx]=min(minoraxisvsx(1:mind,2));
            Mdef_idx=Mdef_idx-1;
            if minoraxisvsx(Mdef_idx,1)-minoraxisvsx(Mdef_idx-1,1)>20
                Mdef_idx=Mdef_idx-1;
            end
            Mdef =minoraxisvsx(Mdef_idx,2);
            med_def= mean(minoraxisvsx(1:5,2));
            xlinks = minoraxisvsx(1:Mdef_idx, 1);
            x_initdefvst=[minoraxisvsx(Mdef_idx, 1),0.01, med_def];
            [m,mind]=min(abs(stenosis_end_x(z)-75-50-2*EL-(minoraxisvsx(:,1))));
            Bdef_idx=mind;
            Bdef =minoraxisvsx(Bdef_idx,2);
            xdatadef = minoraxisvsx(Bdef_idx:Mdef_idx, 1);%vs x
            xvaluesdef{bin_nr}=xdatadef;
            ydatadef = minoraxisvsx(Bdef_idx:Mdef_idx, 2);
            
            lb = [0.5*x_initdefvst(1),0,med_def];
            ub = [1.5*x_initdefvst(1),0.2,med_def];
            fitpara_minor(bin_nr, :) = lsqcurvefit(funcexpminor,x_initdefvst, xdatadef,ydatadef, lb, ub);
            plot(minoraxisvsx(1:Mdef_idx, 1),funcexpminor(fitpara_minor(bin_nr, :),minoraxisvsx(1:Mdef_idx, 1)),'-', 'Color', uni_gruen);
            plot(xlinks,funcexpminor(x_initdefvst,xlinks),'--k');
            plot(minoraxisvsx(Mdef_idx,1),minoraxisvsx(Mdef_idx,2) , 'xk');
            plot(minoraxisvsx(Bdef_idx,1),minoraxisvsx(Bdef_idx,2) ,  'xk');
            
            legend_entry = sprintf('smoothed curve');
            legend_def =  ['$m_0=' , num2str(fitpara_minor(bin_nr, 3)) ,'; \tau_{def}=',num2str(1/fitpara_minor(bin_nr, 2)),'$'];
            legend({legend_entry, legend_def }, 'Interpreter', 'Latex','Location', 'northoutside');
            titlebin =sprintf('class  %d',bin_nr);
            title(titlebin);
            axis([0 800 0 50])
            xlabel('x (pixel)');
            ylabel('minor axis');
            eps_results(:, (bin_nr-2)*6+1)=minor_sorted{bin_nr,2}(:,1);%%ersatz fuer oben
            eps_results(:, (bin_nr-2)*6+2)=minor_sorted{bin_nr,2}(:,2);
            xvaluesminorfit=[55:10:minoraxisvsx(Mdef_idx, 1)].';
            eps_results(1:size(xvaluesminorfit,1), (bin_nr-2)*6+3)=funcexpminor(fitpara_minor(bin_nr, :),xvaluesminorfit);
            subplot(2,5,bin_nr+4);
            plot(major_sorted{bin_nr,2}(:,1),major_sorted{bin_nr,2}(:,2) , '.','Color', uni_grau);
            hold on;
            %DEFORMATION FIT MAJOR AXIS
            majoraxisvsx = major_sorted{bin_nr,2};
            for zeile = size(majoraxisvsx,1):-1:1
                if isnan(majoraxisvsx(zeile, 2))
                    majoraxisvsx(zeile, :)=[];
                end
            end
            [maxdef, Mdef_idx]=max(majoraxisvsx(1:Mdef_idx,2));
            Mdefmajor =majoraxisvsx(Mdef_idx,2);
            med_defmajor= mean(majoraxisvsx(1:5,2));
            xlinks = majoraxisvsx(1:Mdef_idx, 1);
            x_initdefvst=[majoraxisvsx(Mdef_idx, 1),0.01, med_defmajor];
            Bdef =majoraxisvsx(Bdef_idx,2);
            xdatadef = majoraxisvsx(Bdef_idx:Mdef_idx, 1);%vs x
            ydatadef = majoraxisvsx(Bdef_idx:Mdef_idx, 2);
            
            lb = [0.5*x_initdefvst(1),0,med_defmajor];
            ub = [1.5*x_initdefvst(1),0.2,med_defmajor];
            fitpara_major(bin_nr, :) = lsqcurvefit(funcexpmajor,x_initdefvst, xdatadef,ydatadef, lb, ub);
            plot(majoraxisvsx(1:Mdef_idx, 1),funcexpmajor(fitpara_major(bin_nr, :),majoraxisvsx(1:Mdef_idx, 1)),'-', 'Color', uni_gruen);
            plot(xlinks,funcexpmajor(x_initdefvst,xlinks),'--k');
            plot(majoraxisvsx(Mdef_idx,1),majoraxisvsx(Mdef_idx,2) , 'xk');
            plot(majoraxisvsx(Bdef_idx,1),majoraxisvsx(Bdef_idx,2) ,  'xk');
            
            legend_entry = sprintf('smoothed curve');
            legend_def =  ['$m_0=' , num2str(fitpara_major(bin_nr, 3)) ,'; \tau_{def}=',num2str(1/fitpara_major(bin_nr, 2)),'$'];
            legend({legend_entry, legend_def }, 'Interpreter', 'Latex','Location', 'northoutside');
            titlebin =sprintf('class  %d',bin_nr);
            title(titlebin);
            axis([0 800 0 100])
            xlabel('x (pixel)');
            ylabel('major axis');
            eps_results(:, (bin_nr-2)*6+4)=major_sorted{bin_nr,2}(:,2);
            xvaluesmajorfit=[55:10:majoraxisvsx(Mdef_idx, 1)].';
            eps_results(1:size(xvaluesmajorfit,1), (bin_nr-2)*6+5)=funcexpmajor(fitpara_major(bin_nr, :),xvaluesmajorfit);
            eps_results(1:size(xvaluesmajorfit,1), (bin_nr-2)*6+6)= (eps_results(1:size(xvaluesmajorfit,1), (bin_nr-2)*6+5)-2*rpixel(bin_nr))./(2*rpixel(bin_nr));
        end
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf,'color','white')
    export_fig(  filenameellipsefits);
    
    %% SPANNUNGS_DEHNUNGS_DIAGRAMM
    Lkanal = 1100;
    xraster=1:1:Lkanal ;
    L=80; %Stenosenlaenge
    Lgesamt =2*EL+L;
    FS=0*xraster; %Scherkraft
    FC=0*xraster; %Stroemungswiderstandskraft
    mu=0.001; %pas
    r=rpixel*0.55E-6;
    Ap=pi()*r.^2;
    Re=0*xraster; %Reynoldszahl
    h=8*10^-6;
    rho=1;
    FTepsresults=zeros(500,30);
    RateResults=zeros(6,6);
    tkomp=zeros(1,6);
    dsigmadt05=zeros(1,6);
    dsigmadt075=zeros(1,6);
    dsigmadtmax=zeros(1,6);
    %Geschwindigkeit aus bin 3 bzw. kleinstmoeglicher bin, schoenste Kurve, kleinster Fehler, fuer Bestimmung der Kraft F Stenose
    bin_nr=3;
    while ~velok(bin_nr) || isempty(rosettes_sorted_undamaged{bin_nr})
        bin_nr=bin_nr+1;
    end
    velocities_smooth =velocities_sorted{bin_nr,2};
    for zeile = size(velocities_smooth,1):-1:1
        if isnan(velocities_smooth(zeile, 2))
            velocities_smooth(zeile, :)=[];
        end
    end
    velocities_smooth(:, 2)= smoothdata(velocities_smooth(:, 2),'movmedian',5);
    velocities_smooth(:, 2)= smoothdata(velocities_smooth(:, 2),'movmean',5);
    vpixel=interp1(velocities_smooth(:,1), velocities_smooth(:, 2),xraster, 'spline');
    vpixel= smoothdata(vpixel,'movmean',25);
    
    v=vpixel*0.55*10^-6*2000;%pixel/0.5ms in m pro s
    vpixelgauss=gauss(fitpara_vel(bin_nr,:),xraster);
    [max_v, max_vidx] = max(vpixelgauss);
    RateResults(bin_nr,1)=fitpara_vel(bin_nr,1)*0.55*10^-6*2000;
    RateResults(bin_nr,2)=max_v;
    mitte=max_vidx;
    
    %% Kraft vs x und Spannungs-Dehnungs-Diagramm
    for bin_nr=3:5
        if velok(bin_nr) && ~isempty(rosettes_sorted_undamaged{bin_nr})
            startEL=mitte-round(0.5*Lgesamt);
            endEL = startEL+EL;
            deltat=1./v(startEL:endEL);%pixel/m/s
            tkomp(bin_nr)=sum(deltat)*0.55*10^-6;
            d=0*xraster+d1;
            for idx = startEL:endEL
                d(idx)=d1-(xraster(idx)-startEL)*(d1-d2)/EL;
                d(startEL+Lgesamt-(idx-startEL))=d(idx);
            end
            for idx = (endEL+1):(endEL+L-1)
                d(idx)=d2;
            end
            dmum=d*0.55*10^-6;
            
            for idx=1:size(xraster,2)
                FS(idx)=2*pi()*v(idx)*mu*r(bin_nr);
                Re(idx)=(4/3)*v(idx)*(rho/mu)*dmum(idx)*h/(dmum(idx)+h);
                CD=(24/Re(idx))*(1+0.15*Re(idx)^0.681)+0.407/(1+(8710/Re(idx)));
                FC(idx)=0.5*rho*(v(idx)^2)*CD*Ap(bin_nr);
            end
            FT=FC+FS;
            
            filename=sprintf('HysteresisClass%d.png',bin_nr);
            filename=char(folder +filename);
            figure,
            subplot(3, 4,1);
            plot(xraster,dmum, '.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('d(x) in m')
            subplot(3, 4,5);
            plot(xraster,v,'.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('v(x) in m/s')
            subplot(3, 4,2);
            plot(xraster,FS,'.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('FShear(x)')
            subplot(3, 4,9);
            plot(xraster,Re,'.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('Re(x)')
            subplot(3, 4,6);
            plot(xraster,FC,'.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('FCompressive(x)')
            subplot(3,4,10);
            plot(xraster,FT,'.');
            xlim([startEL-50 startEL+Lgesamt+50]);
            legend('FTotal(x)')
            epsfromfit=  nonzeros(eps_results(:, (bin_nr-2)*6+6));
            xdefevenly=eps_results(1, (bin_nr-2)*6+1):1:eps_results(size(epsfromfit), (bin_nr-2)*6+1);
            subplot(3,4,3);
            epsilondef=(funcexpmajor(fitpara_major(bin_nr, :),xdefevenly)-2*rpixel(bin_nr))./(2*rpixel(bin_nr));
            plot(xdefevenly, epsilondef,'.', 'Color', uni_gruen); hold on;
            legend('epsilon(x)')
            xlabel('x(pixel)');
            subplot(3,4,7);
            sigmadef=0*epsilondef;
            for i=1:size(epsilondef,2)
                ahalbeAggr = 0.5*funcexpmajor(fitpara_major(bin_nr, :),xdefevenly(i));
                bhalbeAggr = 0.5*funcexpminor(fitpara_minor(bin_nr, :),xdefevenly(i));
                sigmasum=0;
                ysum=0;
                for xstep=- round(ahalbeAggr)+1:1: round(ahalbeAggr)-1
                    y=  bhalbeAggr*sqrt(1-((xstep/ahalbeAggr)^2));
                    sigmasum=sigmasum+FT((xdefevenly(i)+xstep))/(4*Ap(bin_nr))*y;%faktor 4 weil Flaeche Kugel
                    ysum=ysum+y;
                end
                sigmadef(i)=sigmasum/ysum;
            end
            plot(xdefevenly,sigmadef , '.','Color', uni_gruen);hold on;
            legend('sigma(x)')
            xlabel('x(pixel)');
            subplot(3,4,4);
            plot(epsilondef,sigmadef, '.','Color', uni_gruen);hold on;
            ylabel('sigma');
            xlabel('epsilon');
            legend_entry = sprintf('Class  %d',bin_nr);
            legend(legend_entry)
            subplot(3,4,8);
            deltat=(0.55E-6)./v(xdefevenly(1:end-1));
            dsigmadt=diff(sigmadef)./deltat;
            taxis=cumsum(deltat);
            plot(taxis, dsigmadt); hold on;
            [eps05,eps05ind]=min(abs(0.5-epsilondef));%
            dsigmadt05(bin_nr)=dsigmadt(eps05ind);
            [eps075,eps075ind]=min(abs(0.75-epsilondef));%
            [dsigmadtmax(bin_nr),maxind]=max(dsigmadt);
            if eps075ind <= size(dsigmadt,2)
                dsigmadt075(bin_nr)=dsigmadt(eps075ind);
            elseif eps075ind > size(dsigmadt,2)
                dsigmadt075(bin_nr)=max(dsigmadt);
                eps075ind=maxind;
            end
            
            plot(taxis(eps05ind),dsigmadt(eps05ind),'x', 'Color', uni_gruen);
            plot(taxis(eps075ind), dsigmadt075(bin_nr),'x', 'Color', uni_hellblau);
            plot(taxis(maxind),dsigmadtmax(bin_nr),'x', 'Color', uni_rot);
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+1)=xdefevenly.';
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+2)=FC(xdefevenly);
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+3)=FS(xdefevenly);
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+4)=FT(xdefevenly);
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+5)=epsilondef.';
            FTepsresults(1:size(xdefevenly,2), (bin_nr-2)*6+6)=sigmadef;
            set(gcf, 'Position', get(0, 'Screensize'));
            set(gcf,'color','white')
            export_fig(filename);
        end
    end
    %% Export Ergebnisse
    FTFileName = folder + 'FT_vs_epsilon.txt';
    fid = fopen(FTFileName,'w');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\r\n','Spalte 1-6: Class 2, Spalte 7-12: Class 3,...' );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', ' Spalten: #x  #FC #FS #FT #epsilon #sigma');
    end
    fclose(fid);
    dlmwrite(FTFileName,FTepsresults,'-append','delimiter','\t','newline','pc','roffset',1);
    
    ExpFileName = folder + 'velocity_vs_x.txt';
    fid = fopen(ExpFileName,'w');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\r\n','FITS VELOCITY');
        fprintf(fid,'%s\r\n','   gauss = @(x, xdata)x(1)+x(4)/(x(3)*sqrt(pi/2))*exp( -2 * ((xdata - x(2))/x(3)) .^ 2) ;');
        fprintf(fid,'%s\r\n','Zeile = Class, Spalten: x(1) x(2) x(3) x(4)' );
    end
    fclose(fid);
    dlmwrite(ExpFileName,fitpara_vel,'-append','delimiter','\t','newline','pc','roffset',1);
    fid = fopen(ExpFileName,'a');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','KOMPRIMIERUNGSZEIT');
        fprintf(fid,'%s\r\n','Spalte = Class' );
    end
    fclose(fid);
    dlmwrite(ExpFileName,tkomp,'-append','delimiter','\t','newline','pc','roffset',1);
    fid = fopen(ExpFileName,'a');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', 'Spalten: #x  #v_exp #v_fit #xnorm #v_expnorm #v_fitnorm');
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n','Spalte 1-6: Class 2, Spalte 7-12: Class 3, Spalte 13-18: Class 4, ...' );
        
    end
    fclose(fid);
    dlmwrite(ExpFileName,v_results,'-append','delimiter','\t','newline','pc','roffset',1);
    
    ExpFileName = folder + 'epsilon_vs_x.txt';
    fid = fopen(ExpFileName,'w');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\r\n','FITS DEFORMATION MINOR');
        fprintf(fid,'%s\r\n','funcexpminor =@(x,xdata)-exp(x(2)*(xdata-x(1)))+x(3);');
        fprintf(fid,'%s\r\n','Zeile = Class, Spalten: x(1) x(2) x(3)' );
    end
    fclose(fid);
    dlmwrite(ExpFileName,fitpara_minor,'-append','delimiter','\t','newline','pc','roffset',1);
    fid = fopen(ExpFileName,'a');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\t%s\r\n','results SURF', folder );
        fprintf(fid,'%s\r\n','FITS DEFORMATION MAJOR');
        fprintf(fid,'%s\r\n','funcexpmajor =@(x,xdata)exp(x(2)*(xdata-x(1)))+x(3);');
        fprintf(fid,'%s\r\n','Zeile = Class, Spalten: x(1) x(2) x(3)' );
    end
    fclose(fid);
    dlmwrite(ExpFileName,fitpara_major,'-append','delimiter','\t','newline','pc','roffset',1);
    
    fid = fopen(ExpFileName,'a');%w=write=overwrite existing file
    if fid ~= -1
        fprintf(fid,'%s\r\n','Spalte 1-6: Class 2, Spalte 7-12: Class 3, Spalte 13-18: Class 4, ...' );
        fprintf(fid,'%s\r\n', '===...===');
        fprintf(fid,'%s\r\n', ' Spalten: #x  #minor_exp #minor_fit  #major_exp #major_fit #eps');
    end
    fclose(fid);
    dlmwrite(ExpFileName,eps_results,'-append','delimiter','\t','newline','pc','roffset',1);
end