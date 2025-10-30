clear;      
clc;

%% Input Solution Field 
[data,~,~,~,~,~,~,~,~,~,~] = readnek('WING0.f00001');

[nele,node,~] = size(data);
count = 0;
arry_yz = zeros(2,nele*node);
arry_vw = zeros(2,nele*node);
arry_pt = zeros(2,nele*node);


x_loc = 21.0;       % Location of (Y,Z) plane
tol = 1.0;          % Tolerance should be such that at least one complete (Y, Z) plane is captured
for i = 1:node
    for j = 1:nele
        if (abs(data(j,i,1) - x_loc)<=tol)   
            count = count + 1;
            arry_yz(1,count) = data(j,i,2);
            arry_yz(2,count) = data(j,i,3);
            arry_vw(1,count) = data(j,i,5);
            arry_vw(2,count) = data(j,i,6);
            arry_pt(1,count) = data(j,i,7);
            arry_pt(2,count) = data(j,i,8);
        end
    end
end

arry_yz = arry_yz(:,1:count);
arry_vw = arry_vw(:,1:count);
arry_pt = arry_pt(:,1:count);

pts = arry_yz';

[~, idx] = sortrows(pts,[2 1]);

arry_yz = arry_yz(:,idx);
arry_vw = arry_vw(:,idx);
arry_pt = arry_pt(:,idx);

%% Output Solution Field 
[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('BASvortex_dir0.f00002');

[nele,node,~] = size(data);
count = 0;

for i = 1:nele
    for j = 1:node

        id = binarySearch(arry_yz(2,:),data(i,j,2));

        % Check for boundary points
        flag = 0;
        flag_bc = 0;
        if (abs(data(i,j,2) - arry_yz(2,id)) < 10.0e-8)
            flag_bc = 1;
          if (arry_yz(2,id) == arry_yz(2,length(arry_yz))) 
              id = length(arry_yz);
              if (abs(data(i,j,1) - arry_yz(1,id)) < 10.0e-8)
                      flag = 1;
              else
                  while(data(i,j,1) - arry_yz(1,id) < 0) 
                    id = id -1;
                    if (abs(data(i,j,1) - arry_yz(1,id)) < 10.0e-8)
                        flag = 1;
                        break
                    end
                  end
              end
              
              if (flag == 1)
                  data(i,j,3) = arry_vw(1,id);
                  data(i,j,4) = arry_vw(2,id);
                  data(i,j,5) = arry_pt(1,id);
                  data(i,j,6) = arry_pt(2,id);
              else
              % Interpolation
                  delta_x =  arry_yz(1,id+1) - arry_yz(1,id);
                  delX = data(i,j,1) - arry_yz(1,id);
                  data(i,j,3) = delX*(arry_vw(1,id+1) - arry_vw(1,id))/delta_x +  arry_vw(1,id);
                  data(i,j,4) = delX*(arry_vw(2,id+1) - arry_vw(2,id))/delta_x +  arry_vw(2,id);
                  data(i,j,5) = delX*(arry_pt(1,id+1) - arry_pt(1,id))/delta_x +  arry_pt(1,id);
                  data(i,j,6) = delX*(arry_pt(2,id+1) - arry_pt(2,id))/delta_x +  arry_pt(2,id);
              end
           elseif (arry_yz(2,id) == arry_yz(2,1)) 
              id = 1;
              if (abs(data(i,j,1) - arry_yz(1,id)) < 10.0e-8)
                      flag = 1;           
              else
                  while(data(i,j,1) - arry_yz(1,id) > 0) 
                      id = id +1;
                      if (abs(data(i,j,1) - arry_yz(1,id)) < 10.0e-8)
                          flag = 1;   
                          break
                      end
                  end 
              end
              if (flag == 1)
                  data(i,j,3) = arry_vw(1,id);
                  data(i,j,4) = arry_vw(2,id);
                  data(i,j,5) = arry_pt(1,id);
                  data(i,j,6) = arry_pt(2,id);
              else
              %Interpolation
                  delta_x =  arry_yz(1,id) - arry_yz(1,id-1);
                  delX = data(i,j,1) - arry_yz(1,id-1);
                  data(i,j,3) = delX*(arry_vw(1,id) - arry_vw(1,id-1))/delta_x +  arry_vw(1,id-1);
                  data(i,j,4) = delX*(arry_vw(2,id) - arry_vw(2,id-1))/delta_x +  arry_vw(2,id-1);
                  data(i,j,5) = delX*(arry_pt(1,id) - arry_pt(1,id-1))/delta_x +  arry_pt(1,id-1);
                  data(i,j,6) = delX*(arry_pt(2,id) - arry_pt(2,id-1))/delta_x +  arry_pt(2,id-1);
              end
          end
        end


        if (data(i,j,2) - arry_yz(2,id) > 0 && flag_bc == 0)
          while (data(i,j,2) - arry_yz(2,id) > 0) 
              id = id +1; 
          end
          idx_up = id; idx_low = id-1;
          while(data(i,j,1) - arry_yz(1,idx_up) > 0) 
              idx_up = idx_up +1;
          end
          while(data(i,j,1) - arry_yz(1,idx_low) < 0) 
              idx_low = idx_low -1;
          end
          delta_x = arry_yz(1,idx_up) - arry_yz(1,idx_up-1);
          delX = data(i,j,1) - arry_yz(1,idx_up-1);
          u1 = delX*(arry_vw(1,idx_up) - arry_vw(1,idx_up-1))/delta_x +  arry_vw(1,idx_up-1);
          v1 = delX*(arry_vw(2,idx_up) - arry_vw(2,idx_up-1))/delta_x +  arry_vw(2,idx_up-1);
          p1 = delX*(arry_pt(1,idx_up) - arry_pt(1,idx_up-1))/delta_x +  arry_pt(1,idx_up-1);
          t1 = delX*(arry_pt(2,idx_up) - arry_pt(2,idx_up-1))/delta_x +  arry_pt(2,idx_up-1);

          delta_x = arry_yz(1,idx_low+1) - arry_yz(1,idx_low);
          delX = data(i,j,1) - arry_yz(1,idx_low);
          u2 = delX*(arry_vw(1,idx_low+1) - arry_vw(1,idx_low))/delta_x +  arry_vw(1,idx_low);
          v2 = delX*(arry_vw(2,idx_low+1) - arry_vw(2,idx_low))/delta_x +  arry_vw(2,idx_low);
          p2 = delX*(arry_pt(1,idx_low+1) - arry_pt(1,idx_low))/delta_x +  arry_pt(1,idx_low);
          t2 = delX*(arry_pt(2,idx_low+1) - arry_pt(2,idx_low))/delta_x +  arry_pt(2,idx_low);

          delta_x = arry_yz(2,id) - arry_yz(2,id-1);
          delX = data(i,j,2) - arry_yz(2,id-1);
          data(i,j,3) = delX*(u1 - u2)/delta_x +  u2;
          data(i,j,4) = delX*(v1 - v2)/delta_x +  v2;
          data(i,j,5) = delX*(p1 - p2)/delta_x +  p2;
          data(i,j,6) = delX*(t1 - t2)/delta_x +  t2;


        elseif (data(i,j,2) - arry_yz(2,id) < 0 && flag_bc == 0)  
          while (data(i,j,2) - arry_yz(2,id) < 0) 
            id = id -1; 
          end
          idx_up = id+1; idx_low = id;
          while(data(i,j,1) - arry_yz(1,idx_up) > 0) 
              idx_up = idx_up +1;
          end
          while(data(i,j,1) - arry_yz(1,idx_low) < 0) 
              idx_low = idx_low -1;
          end
          delta_x = arry_yz(1,idx_up) - arry_yz(1,idx_up-1);
          delX = data(i,j,1) - arry_yz(1,idx_up-1);
          u1 = delX*(arry_vw(1,idx_up) - arry_vw(1,idx_up-1))/delta_x +  arry_vw(1,idx_up-1);
          v1 = delX*(arry_vw(2,idx_up) - arry_vw(2,idx_up-1))/delta_x +  arry_vw(2,idx_up-1);
          p1 = delX*(arry_pt(1,idx_up) - arry_pt(1,idx_up-1))/delta_x +  arry_pt(1,idx_up-1);
          t1 = delX*(arry_pt(2,idx_up) - arry_pt(2,idx_up-1))/delta_x +  arry_pt(2,idx_up-1);

          delta_x = arry_yz(1,idx_low+1) - arry_yz(1,idx_low);
          delX = data(i,j,1) - arry_yz(1,idx_low);
          u2 = delX*(arry_vw(1,idx_low+1) - arry_vw(1,idx_low))/delta_x +  arry_vw(1,idx_low);
          v2 = delX*(arry_vw(2,idx_low+1) - arry_vw(2,idx_low))/delta_x +  arry_vw(2,idx_low);
          p2 = delX*(arry_pt(1,idx_low+1) - arry_pt(1,idx_low))/delta_x +  arry_pt(1,idx_low);
          t2 = delX*(arry_pt(2,idx_low+1) - arry_pt(2,idx_low))/delta_x +  arry_pt(2,idx_low);

          delta_x = arry_yz(2,id+1) - arry_yz(2,id);
          delX = data(i,j,2) - arry_yz(2,id);
          data(i,j,3) = delX*(u1 - u2)/delta_x +  u2;
          data(i,j,4) = delX*(v1 - v2)/delta_x +  v2;
          data(i,j,5) = delX*(p1 - p2)/delta_x +  p2;
          data(i,j,6) = delX*(t1 - t2)/delta_x +  t2;
        end
    end
end

%% Write 

status = writenek('newvortex0.f00001',data,lr1,elmap,time,istep,fields,emode,wdsz,etag);

%% Function
function idx = binarySearch(arr,val)
    % arr must be sorted ascending
    lo = 1; 
    hi = numel(arr);
    while lo <= hi
        mid = floor((lo+hi)/2);
        if arr(mid) < val
            lo = mid + 1;
        elseif arr(mid) > val
            hi = mid - 1;
        else
            idx = mid;  % exact match
            return;
        end
    end
    % If not exact, return lower index (for interpolation)
    idx = max(1,hi);
end





