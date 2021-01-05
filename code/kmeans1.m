function [center] = kmeans1( X,N )
[m,n]=size(X);
      pattern=zeros(m,n+1);      
      pattern(:,1:n)=X(:,:);
      center = double(X(randsample(m,N),:));
      while 1
          distance=zeros(1,N);
          num=zeros(1,N);
          new_center=zeros(N,n);
          
          for x=1:m
              for y=1:N
                  distance(y)=norm(X(x,:)-center(y,:));%Distance every point to cluster center
              end
              [~, temp]=min(distance);%
              pattern(x,n+1)=temp;
          end
          k=0;
          for y=1:N
              for x=1:m
                  if pattern(x,n+1)==y
                      new_center(y,:)=new_center(y,:)+pattern(x,1:n);
                      num(y)=num(y)+1;
                  end
              end
              new_center(y,:)=new_center(y,:)/num(y);
              if norm(new_center(y,:)-center(y,:))<0.05
                  k=k+1;
              end
          end
          if k==N
              break;
          else
              center=new_center;
          end
      end
end

