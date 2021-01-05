function [idx, C, sumD, D] = kmeans3(X, k)
% n points in p dimensional space
[n, p] = size(X);

[Xsort,Xord] = sort(X,1);

D = NaN(n,k);   % point-to-cluster distances 
Del = NaN(n,k); % reassignment criterion     
m = zeros(k,1); % 

totsumDBest = Inf;% at first a very large number

C = double(X(randsample(n,k),:)); % Initial cluster centroid

changed = 1:k; % everything is newly assigned 
idx = zeros(n,1);
totsumD = Inf;

iter = 0;
Maxitr=1000;

for rep = 1:10
while true
    % Compute the distance from every point to each cluster centroid     
    D(:,changed) = distfun(X, C(changed,:));

    % Compute the total sum of distances for the current configuration.
    % Can't do it first time through, there's no configuration yet.
    if iter > 0
        totsumD = sum(D((idx-1)*n + (1:n)'));
        % Test for a cycle: if objective is not decreased, back out
        % the last step and move on to the single update phase
       
        if prevtotsumD <= totsumD
            idx = previdx;
            
            [C(changed,:), m(changed)] = gcentroids(X, idx, changed);
            iter = iter - 1;
            % break��1��
            break;%  
        end

        if iter >= Maxitr,               %maxiter
            % break��2��
            break; % while 
        end
    end

    % Determine closest cluster for each point and reassign points to clusters
    
    previdx = idx; % ��СΪn*1
    % totsumD ����ʼ��Ϊ����������ʾ�ܾ���
    prevtotsumD = totsumD;
    % ����ÿһ������С��Ԫ�أ�d�Ĵ�СΪn*1��nidxΪ��СԪ�������е�λ�ã����СΪn*1��DΪn*p
    [d, nidx] = min(D, [], 2);

    if iter == 0
        % iter==0,��ʾ��һ�ε���
        % Every point moved, every cluster will need an update
        % ÿһ������Ҫ�ƶ���ÿһ���ظ���
        moved = 1:n;
        idx = nidx;
        changed = 1:k;
    else
        % Determine which points moved ������һ�����ƶ�
        % �ҵ���һ�κ͵�ǰ��СԪ�ز�ͬ��λ��
        moved = find(nidx ~= previdx);
        if length(moved) > 0
            % Resolve ties in favor of not moving 
            % ���·���������ƶ� ��������һ���߼����� ȷ����Ҫ�ƶ����λ��
            moved = moved(D((previdx(moved)-1)*n + moved) > d(moved));
        end
        % ���û�в�ͬ�ģ�����ǰ������СԪ�أ�����ѭ�����õ���Ԫ���Ѿ��Ǹ��е���Сֵ
        if length(moved) == 0
            % break��3��
            break;
        end
        idx(moved) = nidx(moved);

        % Find clusters that gained or lost members �ҵ���õĻ��߶�ʧ�ĳ�Ա�ķִ�
        % �õ�idx(moved)��previdx(moved)�в��ظ����ֵ�����Ԫ�أ�������������
        changed = unique([idx(moved); previdx(moved)])';
    end

    % Calculate the new cluster centroids and counts. �����µķִ����ĺͼ���
    % C(changed,:)��ʾ�µľ������ģ�m(changed)��ʾ��������idx�г��ֵĴ���
    % sort ��Ԫ�ذ��������У�Xsort��ŶԵ�Ԫ�أ�Xord�д����Ԫ����ԭʼ�����е����ж�Ӧ�Ĵ�Сλ��
    [C(changed,:), m(changed)] = gcentroids(X, idx, changed);
    iter = iter + 1;        
end % phase one

Xmid = zeros([k,p,2]);
for i1 = 1:k
    if m(i1) > 0
        % Separate out sorted coords for points in i'th cluster,
        % and save values above and below median, component-wise
        % �ֽ����i����������ѡ�ĵ�����꣬���������ϣ�����λ��
        % reshape�Ѿ���ֽ�ΪҪ���������m*p
        % sort ��Ԫ�ذ��������У�Xord�д����Ԫ����ԭʼ�����е����ж�Ӧ�Ĵ�Сλ��
        Xsorted1 = reshape(Xsort(idx(Xord)==i1), m(i1), p);
        % floorȡ��ֵС���ߵ��ڵ������ֵ
        nn1 = floor(.5*m(i1));
        if mod(m(i1),2) == 0
            Xmid(i1,:,1:2) = Xsorted1([nn1, nn1+1],:)';
        elseif m(i1) > 1
            Xmid(i1,:,1:2) = Xsorted1([nn1, nn1+2],:)';
        else
            Xmid(i1,:,1:2) = Xsorted1([1, 1],:)';
        end
    end
end
% ------------------------------------------------------------------
% Begin phase two:  single reassignments 

% m�б������ÿһ������ĸ�����Ԫ�غ�Ϊn
% find(m' > 0)�õ�m'�д���0��Ԫ�ص�λ�ã�������
% ʵ�������Ĭ������£�changed=1:k
changed = find(m' > 0);
lastmoved = 0;
nummoved = 0;
iter1 = iter;
while iter < Maxitr
    % Calculate distances to each cluster from each point, and the
    % potential change in total sum of errors for adding or removing
    % each point from each cluster.  Clusters that have not changed
    % membership need not be updated.
   
    % Singleton clusters are a special case for the sum of dists
    % calculation.  Removing their only point is never best, so the
    % reassignment criterion had better guarantee that a singleton
    % point will stay in its own cluster.  Happily, we get
    % Del(i,idx(i)) == 0 automatically for them.
   
    for i2 = changed
            if mod(m(i2),2) == 0 % this will never catch singleton clusters
                ldist = Xmid(repmat(i2,n,1),:,1) - X;
                rdist = X - Xmid(repmat(i2,n,1),:,2);
                mbrs = (idx == i2);
                sgn = repmat(1-2*mbrs, 1, p); % -1 for members, 1 for nonmembers
                Del(:,i2) = sum(max(0, max(sgn.*rdist, sgn.*ldist)), 2);
            else
                Del(:,i2) = sum(abs(X - C(repmat(i2,n,1),:)), 2);
            end
    end
      % Determine best possible move, if any, for each point.  Next we
    % will pick one from those that actually did move.
   
    previdx = idx;
    prevtotsumD = totsumD;
    [minDel, nidx] = min(Del, [], 2);
    moved = find(previdx ~= nidx);
    if length(moved) > 0
        % Resolve ties in favor of not moving
        
        moved = moved(Del((previdx(moved)-1)*n + moved) > minDel(moved));
    end
    if length(moved) == 0
        % Count an iteration if phase 2 did nothing at all, or if we're
        % in the middle of a pass through all the points
        if (iter - iter1) == 0 || nummoved > 0
            iter = iter + 1;

        end

        break;
    end

    % Pick the next move in cyclic order
    moved = mod(min(mod(moved - lastmoved - 1, n) + lastmoved), n) + 1;

    % If we've gone once through all the points, that's an iteration
    if moved <= lastmoved
        iter = iter + 1;
        if iter >= Maxitr, break; end
        nummoved = 0;
    end
    nummoved = nummoved + 1;
    lastmoved = moved;

    oidx = idx(moved);
    nidx = nidx(moved);
    totsumD = totsumD + Del(moved,nidx) - Del(moved,oidx);

    % Update the cluster index vector, and rhe old and new cluster
    % counts and centroids
    idx(moved) = nidx;
    m(nidx) = m(nidx) + 1;
    m(oidx) = m(oidx) - 1;
    for i2 = [oidx nidx]
            % Separate out sorted coords for points in each cluster.
            % New centroid is the coord median, save values above and
            % below median.  All done component-wise.
            Xsorted2 = reshape(Xsort(idx(Xord)==i2), m(i2), p);
            nn2 = floor(.5*m(i2));
            if mod(m(i2),2) == 0
                C(i2,:) = .5 * (Xsorted2(nn2,:) + Xsorted2(nn2+1,:));
                Xmid(i2,:,1:2) = Xsorted2([nn2, nn2+1],:)';
            else
                C(i2,:) = Xsorted2(nn2+1,:);
                if m(i2) > 1
                    Xmid(i2,:,1:2) = Xsorted2([nn2, nn2+2],:)';
                else
                    Xmid(i2,:,1:2) = Xsorted2([1, 1],:)';
                end
            end
    end
    changed = sort([oidx nidx]);
end % phase two

% Calculate cluster-wise sums of distances
nonempties = find(m(:)'>0);
D(:,nonempties) = distfun(X, C(nonempties,:));
d = D((idx-1)*n + (1:n)');
sumD = zeros(k,1);

for i3 = 1:k
    sumD(i3) = sum(d(idx == i3));
end

% Save the best solution so far
if totsumD < totsumDBest
    totsumDBest = totsumD;
    idxBest = idx;
    Cbest = C;
    sumDBest = sumD;
    if nargout > 3
        Dbest = D;
    end
end
end

% Return the best solution
idx = idxBest;
C = Cbest;
sumD = sumDBest;
if nargout > 3
    D = Dbest;
end

function [centroids, counts] = gcentroids(X, index, clusts)
%GCENTROIDS Centroids and counts stratified by group.
p = size(X,2);
num = length(clusts);
centroids = NaN(num,p,'like',X);
counts = zeros(num,1,'like',X);
for i = 1:num
    members = (index == clusts(i));
    if any(members)
        counts(i) = sum(members);          
        % Separate out sorted coords for points in i'th cluster,
        % and use to compute a fast median, component-wise
        Xsorted = sort(X(members,:),1);
        nn = floor(.5*counts(i));
        if mod(counts(i),2) == 0
            centroids(i,:) = .5 * (Xsorted(nn,:) + Xsorted(nn+1,:));
        else
            centroids(i,:) = Xsorted(nn+1,:);
        end      
    end
end
end

function D = distfun(X, C)
%DISTFUN Calculate point to cluster centroid distances.
[n,p] = size(X);
D = zeros(n,size(C,1),'like',X);
nclusts = size(C,1);
for i = 1:nclusts
    D(:,i) = abs(X(:,1) - C(i,1));
    for j = 2:p
        D(:,i) = D(:,i) + abs(X(:,j) - C(i,j));
    end
    % D(:,i) = sum(abs(X - C(repmat(i,n,1),:)), 2);
end
end % function
end