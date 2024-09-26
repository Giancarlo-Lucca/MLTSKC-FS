 function [x_g] = calc_x_gf(x,v,b)
n_examples = size(x,1);
x_e = [x,ones(n_examples,1)];
%eps = 10^(-10);
%eps = 0;
[k,d] = size(v); 

B = [];

p = 2;
q = 2;


for i=1:k
    
    v1 = repmat(v(i,:),n_examples,1);    
    bb = repmat(b(i,:),n_examples,1);
    
    
    v2 = exp(-(x-v1).^2./(2*bb.^2));  
    
    
    v3 = v2-v2;    
    [A,I]=sort(v2,2);
    
    
    q_t = 8;
    switch(q_t)
        case 1
            wt(:,i) = exp(-sum((x-v1).^2./(2*bb),2));       
        case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            wt(:,i) = min(v2,[],2); 
        case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            wt(:,i) = max(v2,[],2);            
        
        case 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medida uniforme  
            for j=1:d
                if j == 1
                    v3(:,j) = A(:,j); 
                    
                end
                if j > 1
                    v3(:,j) = (A(:,j) - A(:,j-1)).*(d-j+1)/d  + eps; 
                    
                    
                end
            end
            
            wt(:,i) = sum(v3,2);
        case 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medida relativa
            for j = 1:d
                if j == 1
                    v3(:, j) = A(:, j);
                    B(j,k) = 1;
                end
                
                if j > 1
                    PP = I;
                    QQ = I(:, j:d);
                    
                    B(j,k) = sum(QQ(j, :)) / sum(PP(j, :));
                    v3(:, j) = (A(:, j) - A(:, j - 1)) * B(j,k) + eps;
                end
            end
            %display(B);
            wt(:, i) = sum(v3, 2);
        case 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medida produto
            for j = 1:d
                if j == 1
                    v3(:, j) = A(:, j);
                    
                end
                
                if j > 1
                    
                    QQ = I(:, 1:j-1);
                    
                    v3(:, j) = (A(:, j) - A(:, j - 1))./prod(QQ,2) + eps;
                end
            end
            
            wt(:, i) = sum(v3, 2);
        case 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medida potencia
            for j=1:d
                if j == 1
                    v3(:,j) = A(:,j);         
                end
                if j > 1
                    v3(:,j) = (A(:,j) - A(:,j-1)).*((d-j+1)/d)^q + eps;    
                end
            end
            wt(:,i) = sum(v3,2);
        case 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% medida ponderada
            %QQ = rand(size(A));
            QQ = A;
            SQQ=sum(QQ,2);
            QQQ=QQ./SQQ;
            for j=1:d
                if j == 1
                    v3(:,j) = A(:,j);        
                end
                if j > 1
                    v3(:,j) = (A(:,j) - A(:,j-1)).*sum(QQQ(:,j:d),2) + eps;
                end
            end
            wt(:,i) = sum(v3,2); 
        
            
    end
    
   
end


wt2 = sum(wt,2);
 
% To avoid the situation that zeros are exist in the matrix wt2
ss = wt2==0;
wt2(ss,:) = eps;
wt = wt./repmat(wt2,1,k);





x_g = [];

for i=1:k
    wt1 = wt(:,i);
    wt2 = repmat(wt1,1,d+1);
    x_g = [x_g,x_e.*wt2];
end
end
