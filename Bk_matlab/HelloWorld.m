function BK_UnitTest
% BK_UnitTest   Compile, load, and test the BK_MATLAB library.
%    BK_UnitTest will make sure the wrapper compiles on the target 
%    platform and then exercises the library to look for silly bugs.

    function Assert(cond,msg)   % for older MATLAB without assert()
        if (exist('assert') == 5)
            if (nargin < 2)
                assert(cond);
            else
                assert(cond,msg);
            end
        elseif (~cond)
            if (nargin < 2)
                msg = 'Assertion failed';
            end
            error(msg);
        end
    end

    function D = Sparse2Dense(S)
        [i,j,s] = find(S);
        z = zeros(size(s,1),1);
        D = [i,j,z,s,s,z];  % [i,j,e00,e01,e10,e11]
    end

BK_BuildLib; disp('BuildLib PASSED');
BK_LoadLib;  disp('LoadLib PASSED');


% Test large scale, and incremental data costs DATA+SMOOTH costs
rand('twister', 987+4); % get the same random stream each time
wd = 128; ht = 96;
noise1 = rand([wd,ht])*20;
noise2 = rand([wd,ht])*20;
H = fspecial('disk',5);
noise1 = imfilter(noise1,H,'replicate');
noise2 = imfilter(noise2,H,'replicate');
noise = [noise1(:)';noise2(:)'];
nb = sparse(wd*ht,wd*ht);
for y=1:ht % set up a grid-like neighbourhood, arbitrarily
    for x=1:wd
        if (x < wd), nb((y-1)*wd+x,(y-1)*wd+x+1) = 1; end
        if (y < ht), nb((y-1)*wd+x, y   *wd+x  ) = 1; end
    end
end
distmap = [zeros(ht,wd/2) ones(ht,wd/2)];
distmap = bwdist(distmap);
distmap = distmap / max(distmap(:));
distmap1 = distmap;
distmap2 = flipdim(distmap,2);
distmap = [distmap1(:)'; distmap2(:)'];

hinc = BK_Create(wd*ht,2*wd*ht);
BK_SetNeighbors(hinc,nb);

time_inc = [];
time_new = [];

figure;

lambda = 16;
while (lambda >= 1)
    newdc = double(noise+lambda*distmap);
   
    BK_SetUnary(hinc,newdc);
    tic; 
        e_inc = BK_Minimize(hinc);
    time_inc(end+1) = toc;
    lab = BK_GetLabeling(hinc);
    
    imagesc(reshape(lab,[wd,ht]));
    drawnow;
    
    hnew = BK_Create(wd*ht,2*wd*ht);
    BK_SetNeighbors(hnew,nb);
    BK_SetUnary(hnew,newdc); 
    tic; 
        e_new = BK_Minimize(hnew);
    time_new(end+1) = toc;
    BK_Delete(hnew);
    
    Assert(abs(e_inc - e_new) < 1e-6);
    lambda = lambda*0.9;
end

BK_Delete(hinc);
fprintf('Dynamic PASSED (%.3fsec normal, %.3fsec dynamic)\n',sum(time_new),sum(time_inc)); 

figure; plot(time_new,'-or'); hold on; plot(time_inc,'-xb');


end

