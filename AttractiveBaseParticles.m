classdef AttractiveBaseParticles <handle
    properties
        c_b;c_p;c_t=0.75;
        r_d = 0.5;
        L = 50;
        p=0.75;
        dt = 0.01;steps;
        sigma = 1; % thermal coeffcient huge
        sigma_a = 0.02; k =5;
        pts_b;pts_p;
        save_pos;
        idx; F_max;
        tau = 1;
    end
    properties(Dependent)
        R_eff_vol,R_eff_attract,r_e;
    end
    
    methods
        function obj = AttractiveBaseParticles(c_b,c_p,steps)
            obj.c_b = c_b; obj.c_p = c_p; obj.steps = steps;
            N_t = floor(obj.L^2 / (pi * obj.r_d^2)*obj.c_t);
            N_b = floor(N_t*obj.c_t*c_b);
            N_p =  floor(N_t*obj.c_t*c_p);
            obj.pts_b = poissonDisc([obj.L,obj.L],2*obj.R_eff_vol,N_b);
            obj.pts_p = poissonDisc([obj.L,obj.L],2*obj.R_eff_vol,N_p);        
            save_pos = cell(steps,1);
            obj.idx = [1:N_p]';
            obj.F_max = obj.r_d /5/ obj.dt;%obj.r_d/5 / obj.dt;
        end 
        
        function traj_save = step(obj)        
            num = size(obj.pts_p,1);   
            mask = zeros(num,1);
           % traj_save = zeros(num * obj.steps,4);
            for m = 1:obj.steps      
                % random force             
                r_random = 2*rand(num,2)-1;
                F_random = obj.sigma * r_random./sqrt(sum(r_random.^2,2));
                F_attract = arrayfun(@(x,y) obj.get_attract_force(x,y),...
        obj.pts_p(:,1),obj.pts_p(:,2),'UniformOutput',false);
                F_volume = arrayfun(@(x,y) obj.get_volume_force(x,y),...
        obj.pts_p(:,1),obj.pts_p(:,2),'UniformOutput',false);
                
                delta_r = (F_random + cell2mat(F_attract) + cell2mat(F_volume)) * obj.dt;
                %delta_r = F_random  * obj.dt;
                %update pos
                mask = arrayfun(@(x,y,cnt) obj.judge_absorbed(x,y,cnt),...
                   obj.pts_p(:,1),obj.pts_p(:,2),mask); %absorbed
                %judge whether absorbed             
                %new_pos = obj.pts_p +(~mask).* delta_r;
                %mask2 = arrayfun(@(x,y) obj.judge_absorbed(x,y),...
                   %new_pos(:,1),new_pos(:,2));% new_absorbed
                   
                   
                obj.pts_p = obj.pts_p +(mask >= 0).*(~mask).* delta_r + -(mask < 0).* delta_r;
              %  traj_save((m-1)*num + 1 : m*num,:) = [obj.idx,ones(num,1)*m,...
               %     obj.pts_p];
                %
                if mod(m,100) == 0              
                    figure(m)
                    scatter(obj.pts_b(:,1),obj.pts_b(:,2),10,'filled')           
                    hold on
                    scatter(obj.pts_p(:,1),obj.pts_p(:,2),10)     
                    xlim([0,obj.L]);ylim([0,obj.L])
                    set(gcf,'Visible','off')               
                    saveas(gcf,['./results/',num2str(m),'.tif']);
                    disp('figure plot finished')
                    close(gcf)
                    disp(m)                   
                end
                
            end
        end
        
        
        function F = get_attract_force(obj,x,y)
            %1/r^2
         
            neighbor = obj.find_neighbor(x,y,obj.R_eff_attract,obj.pts_b);
           % r = bsxfun(@minus, [x,y],neighbor);
            r = [x,y] - neighbor;
            r_abs = sqrt(sum(r.^2,2));
            %r_unit = r ./r_abs;
            F = -obj.sigma_a * sum(heaviside(obj.R_eff_attract - r_abs)...
            * 1./max(r_abs.^2,obj.r_d) .* r./r_abs,1);         
            if sqrt(sum(F.^2,2)) > obj.F_max
            F = obj.F_max * F /sqrt(sum(F.^2,2)); 
            end
        end
        
        function F = get_volume_force(obj,x,y)
            neighbor = obj.find_neighbor(x,y,obj.R_eff_vol,obj.pts_p);
           % r = bsxfun(@minus, [x,y],neighbor);
            r = [x,y] - neighbor;
            r_abs = sqrt(sum(r.^2,2));
            %r_unit = r ./r_abs;
            F = obj.k/obj.r_d * sum(heaviside(obj.R_eff_vol - r_abs) ...
            .* (obj.R_eff_vol - r_abs) .* r./r_abs,1);
        end
        
        function cnt =judge_absorbed(obj,x,y,cnt)  
            % if absorbed 
            if cnt >0
                if mod(cnt,floor(obj.tau/obj.dt))==0 && rand(1) > obj.p %desorption
                    cnt = -1;
                else
                    cnt = cnt + 1;
                end
            elseif cnt ==-1 %desorption to free
                cnt = -double(~isempty(obj.find_neighbor(x,y,1.5*obj.r_d,obj.pts_b))); %in the circle
            else               
            % create  
                  cnt = double(~isempty(obj.find_neighbor(x,y,obj.r_e,obj.pts_b))....
                      &    rand(1) <= obj.p );       
                %create the absorption state            
            end
            
            %judge = (~isempty(obj.find_neighbor(x,y,obj.r_e,obj.pts_b)) & rand(1) <= obj.p);
            
%             if cnt 
%                  disp('absorbed')
%              end
        end
        
    
            
         function R_eff_vol = get.R_eff_vol(obj)
            R_eff_vol = 2 * obj.r_d;
         end
        function R_eff_attract = get.R_eff_attract(obj)
            R_eff_attract = 20 * obj.r_d;
        end
        function r_e = get.r_e(obj)
            r_e = 0.1 * obj.r_d;
        end
        
        function plot_every_step(obj)
        end
        function plot_save_trajectory(obj)
        end
         
    end  
    methods(Access = protected)
        function neighbor = find_neighbor(obj,x,y,R_eff,data_set)
%x array sorted
                idx = find(data_set(:,1)< x + R_eff & ...
                data_set(:,1) >= x - R_eff );
                temp = data_set(idx,:);
                idy = find(temp(:,2) < y + R_eff  & temp(:,2) >= y - R_eff);
                neighbor = temp(idy,:);
                neighbor((neighbor(:,1) == x & neighbor(:,2) ==y),:) = [];
        end
    end
end

