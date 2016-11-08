%Global robustness
count=0;
x=ones(834,1);
stdtrials=0:1:50;
matrix_trials2=zeros(length(stdtrials), 500)
for k=1:length(stdtrials)
    for j=1:500
        load('iMLTC834cdf_final.mat')
        for i=1:length(fbamodel.rxns)
        R1=normrnd(0,stdtrials(k));
        R2=normrnd(0,stdtrials(k));
        fbamodel.vmax(i)=fbamodel.vmax(i)+R1;
        if fbamodel.rev==1
            fbamodel.vmin(i)=fbamodel.vmin(i)+R2;
        end
        end

    tmp=evaluate_objective(x,2,numel(geni),fbamodel,geni,reaction_expression);

%b1=-130
%b2=17467
    A(j)=tmp(1);
    B(j)=tmp(2);
        
    
    
    end
    matrix_trials2(k,:)=A;
    matrix_trialsAAA(k,:)=B;

end
J=matrix_trialsAAA;
J=J/100;



stdv=0:1:50;
epsilon=0:0.0001:0.01;
matrix_robust=zeros(length(stdv),length(epsilon));
for k=2:size(J,1)
    for j=1:length(epsilon)
        C=abs(J(k,:)+1.42512381)./1.42512381;
        tmp=0;
        for i=1:length(C)
            if C(i)<=epsilon(j)
                tmp=tmp+1;
            end
        end
        matrix_robust(k-1,j)=(tmp/500)*100;
    end
    %matrix_robust(k,:)=countrobust_trials(j);
end
epsilon=epsilon*100;
stdv=((stdv.*3)./1000).*100; 
[X,Y]=meshgrid(epsilon,stdv);



figure;
surf(X,Y,matrix_robust)
xlabel('\epsilon (%)','FontSize',18)
ylabel('\sigma (%)', 'FontSize',18)
zlabel('GR','FontSize',14)
colorbar
