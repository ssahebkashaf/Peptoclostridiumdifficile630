%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -Optimisation de la distribution de l'eau sur un pÈrimËtre secondaire du-
% ---------------canal de Gignac-Saint AndrÈ et Ceyras---------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---Equation : J=C'x+1/2x'Qx et Ax >= or <= or = b --> J=f'x+1/2x'Hx et---
% ---------------------Aineqx <=bineq et Aeqx = beq------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dt is the time slot (mn). When dt is small, the time slot number
% increase(nb_dt). when nb_dt increases, the decision variables and
% constraints increase but the problems are still the same. It just gives more
% calculating precision when dt is small.
% When dt =5, i got the message as i ask you. if dt=10, it works normally.

dt=60; % pas de temps de calcul (mn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------------------DonnÈes--------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pr=[1 1 1 1 1 1 1 0 0 60 30;
    1 2 1 1 1 1 1 1 0 180 30;
    2 1 1 1 1 1 1 2 0 120 30;
    2 2 1 1 1 1 1 3 0 60 20];
% %     3 1 1 1 1 1 1 2 0 120 30;
% %     3 1 2 1 1 1 1 3 0 60 20;
% %     4 1 1 1 1 1 1 2 0 120 30;
% %     4 2 1 1 1 1 1 3 0 60 20];
% Pr=[1 1 1 1 1 1 1 0 0 60 30;
%     2 1 1 1 1 1 1 1 0 180 30;
%     2 2 1 1 1 1 1 2 0 120 30;
%     2 2 2 1 1 1 1 3 0 60 20;
%     2 3 1 1 1 1 1 4 0 60 30;
%     3 1 1 1 1 1 1 5 0 60 30;
%     4 1 1 1 1 1 1 6 0 120 50;
%     5 1 1 1 1 1 1 7 0 60 30;
%     5 1 2 1 1 1 1 8 0 60 30;
%     5 2 1 1 1 1 1 9 0 60 30;
%     6 1 1 1 1 1 1 10 0 60 35]; % Demandes des prises, format [i j k coef_cult coef_tech coef_irri td_j td_h td_mn dd qd]
nb_van=2;% Nombre des vannes
d_tour=0.5; % 0.5 jour du tour d'eau
Qr=[1 0 15;1 8 10;1 10 10;1 12 0]; % DÈbit en tÍte du canal secondaire [jour, heure, dÈbit(l/s)] 
t_inter=[1 0 2; 1 4 11]; % (h) Temps d'intervention [jour, ts, tf], ts=temps de commencement, tf= temps de fin du travail d'une intervention du jour
t_mo_max=30; % 30mn
bief2i=[1 1;2 2]; % [n? bief, valeur i qui correspond ? la limite du bief]
Qc_bief=[1 70;2 50]; %[n? bief, capacitÈ du bief l/s] 
Qci=[1 35;2 35]; % capacitÈ du canal secondaire, format [i Qci; i Qci;...]
Qcij=[1 1 35; 2 1 35]; % capacitÈ du canal tertiaire, format [i j Qci; i j Qcij;...]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------ParamËtres de calcul---------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w=[1/6 1/6 1/6 1/6 1/6 1/6]; % co?t d'optimisation
q_lmin=50/100; % facteur de limite du dÈbit minimum
q_lmax=125/100; % facteur de limite du dÈbit maximum 
d_lmin=1/q_lmax; % facteur de limite de la durÈe minimale
d_lmax=1/q_lmin; % facteur de la limite de la durÈe maximale
M=100; % grande valeur positive
%dt=10; % pas de temps de calcul (mn)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----Extraction des donnÈes au format de calcul et des coefficients -----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calcul des nombres de i

n=max(Pr(:,1));

% Calcul des nombres de j qui dÈpend de i
m=zeros(1,1);

for i=1:n
	xx=1;
	j2i=zeros(size(Pr,1),1);
        for prise=1:size(Pr,1)
            if Pr(prise,1)==i            
                j2i(xx)=Pr(prise,2);
                xx=xx+1;
            end        
        end
	m(i)=max(j2i);
end
clear j2i

% Calcul des nombres de k qui dÈpend de i et j
p=zeros(i,max(m));

for i=1:n
    for j=1:m(i)
    	xx=1;
        k2ij=zeros(size(Pr,1),1);
        for prise=1:size(Pr,1)
            if Pr(prise,1)==i && Pr(prise,2)==j                    
            	k2ij(xx)=Pr(prise,3);
            	xx=xx+1;
            end
       	end
        p(i,j)=max(k2ij);
    end
end
    clear k2ij

% Calcul des nombres de prise    
nb_prise=sum(sum(p));
            
% Calcul des nombres du pas de temps
nb_dt=floor(d_tour*24*60/dt+1); % Nombre de pas de temps

% Calcul du Qt (Q au pas de temps dt)
Qrt=zeros(nb_dt,1);

    for xx=1:size(Qr,1)-1;
        for t=1:nb_dt
            if t>=(Qr(xx,1)-1)*24*60/dt+Qr(xx,2)*60/dt && t<=(Qr(xx,1)-1)*24*60/dt+Qr(xx+1,2)*60/dt            
                Qrt(t)=Qr(xx,3);
            end
        end
    end

% Conversion des donnÈes au format du calcul et Calcul des coefficients 
td=zeros(nb_prise,1);
dd=zeros(nb_prise,1);
qd=zeros(nb_prise,1);
qdv=zeros(nb_prise,1);
ts_emax2=zeros(nb_prise,1);
d_emax2=zeros(nb_prise,1);
q_emax2=zeros(nb_prise,1);
qv_emax=zeros(nb_prise,1);
coef_t=zeros(nb_prise,1);
coef_q=zeros(nb_prise,1);
coef_qv=zeros(nb_prise,1);
alpha=zeros(nb_prise,1);
beta=zeros(nb_prise,1);
mphi=zeros(nb_prise,1);
cphi=zeros(nb_prise,1);
const=zeros(nb_prise,1);

for prise=1:nb_prise
        
    % Conversion des donnÈes au format du calcul
    td(prise)=floor(((Pr(prise,7)-1)*60*24+Pr(prise,8)*60 + Pr(prise,9))/dt+1);
    dd(prise)=ceil(Pr(prise,10)/dt);
    qd(prise)=Pr(prise,11);     

    % Calcul des volumes demandÈs
    qdv(prise)=qd(prise)*dd(prise);
        
    % Calcul des Ècarts maximums 
    ts_emax2(prise)=(max(nb_dt-td(prise)-ceil((d_lmin*dd(prise))),td(prise)-1))^2;
    d_emax2(prise)=(dd(prise)-ceil(d_lmax*dd(prise)))^2;
    q_emax2(prise)=(qd(prise)-ceil(q_lmin*qd(prise)))^2;       
    qv_emax(prise)=qdv(prise)-ceil(q_lmin*qdv(prise));
        
    % Calcul des coefficients prioritaires des prises sur t, q et qv
    coef_t(prise)=1/4*(Pr(prise,4)/max(Pr(:,4))+Pr(prise,5)/max(Pr(:,5))+Pr(prise,6)/max(Pr(:,6))+prise/nb_prise);
    coef_q(prise)=Pr(prise,5)/max(Pr(:,5));
    coef_qv(prise)=Pr(prise,4)/max(Pr(:,4));
        
    % Calcul des coeffients de la fonction d'objectif
    alpha(prise)=w(1)*coef_t(prise)/sum(ts_emax2);
    beta(prise)=w(2)/sum(d_emax2);
    mphi(prise)=w(3)*coef_q(prise)/sum(q_emax2);
    cphi(prise)=w(4)*coef_qv(prise)/sum(qv_emax);
    const(prise)=alpha(prise)*td(prise)^2+beta(prise)*dd(prise)^2+mphi(prise)*qd(prise)^2+cphi(prise)*qdv(prise); % constant pour la fonction d'objectif
end

% Calcul des temps de l'intervention
d_inter=zeros(size(t_inter,1),1);
t_mos=zeros(size(t_inter,1),1);
t_mof=zeros(size(t_inter,1),1);
mo_max=zeros(size(t_inter,1),1);

for inter=1:size(t_inter,1)                    
    d_inter(inter)=(t_inter(inter,3)-t_inter(inter,2))*60;    
    t_mos(inter)=floor((t_inter(inter,1)-1)*24*60/dt+t_inter(inter,2)*60/dt+1);
    t_mof(inter)=floor((t_inter(inter,1)-1)*24*60/dt+t_inter(inter,3)*60/dt+1);
    mo_max(inter)=floor(d_inter(inter)/t_mo_max);
end     
clear Pr   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---- CrÈer la matrice C et Q (coefficient de la fonction d'objectif)-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Matrice C

    % Variable t, d et qs
    Cts=zeros(nb_prise,1);
    Cdd=zeros(nb_prise,1);
    Cqs=zeros(nb_prise,1);
    
    n_lig=1; % NumÈro de ligne de la matrice C ? remplir des valeurs
    for prise=1:nb_prise;                           
        Cts(n_lig,1)=-2*alpha(prise)*td(prise) ;
        Cdd(n_lig,1)=-2*beta(prise)*dd(prise);
        Cqs(n_lig,1)=-2*mphi(prise)*qd(prise);                                                          
    n_lig=n_lig+1;                
    end
    
    % Variable q
    Cq=zeros(nb_prise*nb_dt,1);
    
    n_lig=1;
    for t=1:nb_dt
        for prise=1:nb_prise
           Cq(n_lig,1)=-cphi(prise);
    n_lig=n_lig+1;
        end
    end
    
    % Variable mo (Cmo)et pe
    Cmo=zeros(nb_van*nb_dt,1);
    Cpe=zeros(nb_van*nb_dt,1);
    
    n_lig=1;
    for t=1:nb_dt
        for v=1:nb_van
            Cmo(n_lig,1)=w(5)/sum(mo_max);
            Cpe(n_lig,1)=w(6)/sum(Qrt);
    n_lig=n_lig+1;
        end   
    end

    % Variable  pe et y--> Cpe=Cy=0;

    % SynthËse, format C(:,1)-->f pour Cplex sovler
    f=[Cts;Cdd;Cqs;Cq;Cmo;Cpe;zeros(nb_prise*nb_dt,1)];
    clear Cts Cdd Cqs Cq Cmo Cpe

% Matrice Q
Qts=zeros(nb_prise,nb_prise);
Qdd=zeros(nb_prise,nb_prise);
Qqs=zeros(nb_prise,nb_prise);

n_lig=1;
    for prise=1:nb_prise
        n_col=1;
        for prise1=1:nb_prise
            if prise1==prise
                Qts(n_lig,n_col)=alpha(prise);
                Qdd(n_lig,n_col)=beta(prise);                
                Qqs(n_lig,n_col)=mphi(prise);
            end
        n_col=n_col+1;
        end               
n_lig=n_lig+1;
    end

    % SynthËse Q --> H pour Cplex solver
    Q=[Qts sparse(size(Qts,1),(2*nb_prise+2*nb_prise*nb_dt+(2*nb_van)*nb_dt));
        sparse(size(Qdd,1),nb_prise) Qdd sparse(size(Qdd,1),(nb_prise+2*nb_prise*nb_dt+ (2*nb_van)*nb_dt));
        sparse(size(Qqs,1),2*nb_prise) Qqs sparse(size(Qqs,1),(2*nb_prise*nb_dt+ (2*nb_van)*nb_dt));
        sparse(2*nb_prise*nb_dt+(2*nb_van)*nb_dt,(3*nb_prise+2*nb_prise*nb_dt+(2*nb_van)*nb_dt))];
    clear Qts Qdd Qqs 
    H=(Q+Q');% format symÈtrique
    %clear Q

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------- CrÈer la matrice A, b et sense----------------------------% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Contraintes de tous les ijk (c2)
Adc2=zeros(nb_prise,nb_prise);

for prise=1:nb_prise                            
    n_col=1;    
    for prise1=1:nb_prise
        if prise1==prise
        Adc2(prise,n_col)=1;
        end
    n_col=n_col+1;  
    end
end
 
% Contrainte de tous les ijk qui dÈpend de t
bc2=zeros(nb_prise,1);
bc12=zeros(nb_prise,1);
bc13=zeros(nb_prise,1);

Ayc2=zeros(nb_prise,nb_prise*nb_dt);
Aqc12=zeros(nb_prise,nb_prise*nb_dt);
Aqc13=zeros(nb_prise,nb_prise*nb_dt);

for prise=1:nb_prise    
    
    % vecteur b
    bc2(prise,1)=0;
    bc12(prise,1)=qdv(prise);
    bc13(prise,1)=q_lmin*qdv(prise);
            
    % vecteur de sense            
%     sensec2(prise,1)=0;
%     sensec12(prise,1)=-1;
%     sensec13(prise,1)=1;
            
    n_col=1;
    for t=1:nb_dt
        for prise1=1:nb_prise
            if prise1==prise
                % Contrainte 2
                Ayc2(prise,n_col)=-1;
                % Contrainte 12
                Aqc12(prise,n_col)=1;
                % Contrainte 13
                Aqc13(prise,n_col)=1;
            end
            n_col=n_col+1;
        end
    end
end

% Contrainte de tous les t et ijk
bc3=zeros(nb_prise*nb_dt,1);
bc4=zeros(nb_prise*nb_dt,1);
bc5=zeros(nb_prise*nb_dt,1);
bc6=zeros(nb_prise*nb_dt,1);
bc7=zeros(nb_prise*nb_dt,1);

Atsc3=sparse(nb_prise*nb_dt,nb_prise);
Ayc3=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Atsc4=sparse(nb_prise*nb_dt,nb_prise);
Adc4=sparse(nb_prise*nb_dt,nb_prise);
Ayc4=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Aqsc5=sparse(nb_prise*nb_dt,nb_prise);
Aqc5=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Ayc5=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Aqsc6=sparse(nb_prise*nb_dt,nb_prise);
Aqc6=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Ayc6=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Aqc7=sparse(nb_prise*nb_dt,nb_prise*nb_dt);
Ayc7=sparse(nb_prise*nb_dt,nb_prise*nb_dt);

n_lig=1;
for t=1:nb_dt
    for prise=1:nb_prise
            
    % vecteur b            
    bc3(n_lig,1)=t+M;
    bc4(n_lig,1)=t-M+1;
    bc5(n_lig,1)=-M;
    bc6(n_lig,1)=M;
    bc7(n_lig,1)=0;
                            
    % vecteur de sense            
%     sensec3(n_lig,1)=-1;
%     sensec4(n_lig,1)=1;
%     sensec5(n_lig,1)=1;
%     sensec6(n_lig,1)=-1;
%     sensec7(n_lig,1)=-1;
            
    n_col=1;
        for prise1=1:nb_prise
            if prise1==prise
            % Contrainte 3                                      
            Atsc3(n_lig,n_col)=1;                                        
            Ayc3(n_lig,(t-1)*nb_prise+n_col)=M;
            % Contrainte 4
            Atsc4(n_lig,n_col)=1;
            Adc4(n_lig,n_col)=1;
            Ayc4(n_lig,(t-1)*nb_prise+n_col)=-M;
            % Contrainte 5
            Aqsc5(n_lig,n_col)=-1;
            Aqc5(n_lig,(t-1)*nb_prise+n_col)=1;
            Ayc5(n_lig,(t-1)*nb_prise+n_col)=-M;
            % Contrainte 6
            Aqsc6(n_lig,n_col)=-1;
            Aqc6(n_lig,(t-1)*nb_prise+n_col)=1;
            Ayc6(n_lig,(t-1)*nb_prise+n_col)=M;
            % Contrainte 7
            Aqc7(n_lig,(t-1)*nb_prise+n_col)=1;
            Ayc7(n_lig,(t-1)*nb_prise+n_col)=-M;
            end                            
        n_col=n_col+1;       
        end           
n_lig=n_lig+1;
    end
end

% Contrainte de tous les t qui dÈpend de t-1
bc8=zeros(nb_van,1);
bc9=zeros(nb_van,1);
bc10=zeros(nb_van*(nb_dt-1),1);
bc11=zeros(nb_van*(nb_dt-1),1);

Aqc8=zeros(nb_van,nb_prise*nb_dt);
Amoc8=zeros(nb_van,nb_van*nb_dt);
Amoc9=zeros(nb_van,nb_van*nb_dt);
Apec9=zeros(nb_van,nb_van*nb_dt);
Aqc10=zeros(nb_van*(nb_dt-1),nb_prise*nb_dt);
Apec10=zeros(nb_van*(nb_dt-1),nb_van*nb_dt);
Amoc10=zeros(nb_van*(nb_dt-1),nb_van*nb_dt);
Aqc11=zeros(nb_van*(nb_dt-1),nb_prise*nb_dt);
Apec11=zeros(nb_van*(nb_dt-1),nb_van*nb_dt);
Amoc11=zeros(nb_van*(nb_dt-1),nb_van*nb_dt);

n_lig=1;
for t=1:nb_dt
    if t==1
        v=1;
        for i=1:n            
            if m(i)>1
            
                % matrice b
                bc8(n_lig,1)=0;
                bc9(n_lig,1)=0;
                % matrice de sense            
%                 sensec8(n_lig,1)=-1;
%                 sensec9(n_lig,1)=-1;
            
                n_col=1;            
                for t1=1:nb_dt                
                    for i1=1:n
                        for j1=1:m(i1)
                            for k1=1:p(i1,j1)                             
                                if t1==t && i1==i                                
                                % contrainte 8                                    
                                Aqc8(n_lig,n_col)=1;
                                end
                n_col=n_col+1;
                            end
                        end
                    end
                end

                n_col=1;
                for t1=1:nb_dt
                    for v1=1:nb_van
                        if t1==t && v1==v
                        % contrainte 8
                        Amoc8(n_lig,n_col)=-M;                    
                        % contrainte 9
                        Amoc9(n_lig,n_col)=-M;
                        Apec9(n_lig,n_col)=1;
                        end
                n_col=n_col+1;
                    end
                end
v=v+1;
n_lig=n_lig+1;
            end       
        end

n_lig=1;    
    elseif t>=2
        v=1;
        for i=1:n            
            if m(i)>1
                % matrice b
                bc10(n_lig,1)=0;
                bc11(n_lig,1)=0;
            
                % matrice sense        
%                 sensec10(n_lig,1)=1;
%                 sensec11(n_lig,1)=-1;
            
                n_col=1;            
                for t1=1:nb_dt                
                    for i1=1:n
                        for j1=1:m(i1)
                            for k1=1:p(i1,j1)                              
                                if t1==t-1 && i1==i                                
                                % contrainte 10                                    
                                Aqc10(n_lig,n_col)=-1;                                    
                                % Contrainte 11
                                Aqc11(n_lig,n_col)=-1;
                n_col=n_col+1;
                                elseif t1==t && i1==i        
                                Aqc10(n_lig,n_col)=1;
                                Aqc11(n_lig,n_col)=1;
                n_col=n_col+1;  
                                else                                    
                                    Aqc10(n_lig,n_col)=0;
                                    Aqc11(n_lig,n_col)=0;
                n_col=n_col+1;
                                end
                            end
                        end
                    end              
                end
            
                n_col=1;
                for t1=1:nb_dt 
                    for v1=1:nb_van
                        if t1==t-1 && v1==v                     
                        Apec10(n_lig,n_col)=-1;
                        Apec11(n_lig,n_col)=-1;
                n_col=n_col+1;                       
                        elseif t1==t && v1==v
                        Amoc10(n_lig,n_col)=M;
                        Apec10(n_lig,n_col)=1;
                        Amoc11(n_lig,n_col)=-M;
                        Apec11(n_lig,n_col)=1;
                n_col=n_col+1;             
                        else                        
                        Amoc10(n_lig,n_col)=0;
                        Apec10(n_lig,n_col)=0;
                        Amoc11(n_lig,n_col)=0;
                        Apec11(n_lig,n_col)=0;
                n_col=n_col+1;  
                        end
                    end
                end
v=v+1;
n_lig=n_lig+1;
            end 
        end
    end
end           
                    
% Contrainte 15 (CapacitÈ de la ressource)
bc15=zeros(nb_dt,1);

Aqc15=zeros(nb_dt,nb_prise*nb_dt);
Apec15=zeros(nb_dt,nb_van*nb_dt);

n_lig=1;
for t=1:nb_dt
    % matrice b
    bc15(n_lig,1)=Qrt(t);
            
    % matrice sense        
%     sensec15(n_lig,1)=-1;
            
    n_col=1;
    for t1=1:nb_dt
        for prise=1:nb_prise
            if t1==t
            Aqc15(n_lig,n_col)=1;
            end
    n_col=n_col+1;                    
        end
    end

    n_col=1;
    for t1=1:nb_dt
        for v=1:nb_van
            if t1==t 
            Apec15(n_lig,n_col)=1;
            end
    n_col=n_col+1;  
        end
    end
n_lig=n_lig+1;
end
    
% Contrainte 16 ( contrainte de la capacitÈ du canal en fonction du bief)
bc16=zeros(1,1);

Aqc16=zeros(1,nb_prise*nb_dt);

n_lig=1;
for t=1:nb_dt
    for bief=1:max(bief2i(:,1))
        
        % matrice b
        bc16(n_lig,1)=Qc_bief(bief,2);
                
        % matrice sense
%         sensec16(n_lig,1)=-1;
        if bief==1
            n_col=1;        
            for t1=1:nb_dt
                for i=1:n
                    for j=1:m(i)
                        for k=1:p(i,j)             
                            if t1==t 
                                Aqc16(n_lig,n_col)=1;
                            end    
            n_col=n_col+1;               
                        end
                    end
                end
            end
n_lig=n_lig+1;
        else
            n_col=1;        
            for t1=1:nb_dt
                for i=1:n
                    for j=1:m(i)
                        for k=1:p(i,j)             
                            if t1==t && i>bief2i(bief-1,2) 
                                Aqc16(n_lig,n_col)=1;
                            end    
            n_col=n_col+1;               
                        end
                    end
                end
            end
n_lig=n_lig+1;
        end
    end
end
                            
% Contrainte 17 ( contrainte de la capacitÈ de branche i, j=m(i)>1 )
bc17=zeros(1,1);

Aqc17=zeros(1,nb_prise*nb_dt);
Apec17=zeros(1,nb_van*nb_dt);

if max(m)>1
    n_lig=1;
    for t=1:nb_dt
        for br2i=1:size(Qci,1)
           
            % matrice b
            bc17(n_lig,1)=Qci(br2i,2);
     
            % matrice sense
%         sensec17(n_lig,1)=-1;
    
            n_col=1;
        
            for t1=1:nb_dt
                for i=1:n
                    for j=1:m(i)
                        for k=1:p(i,j)             
                            if t1==t && i==Qci(br2i,1)
                                Aqc17(n_lig,n_col)=1;
                            end
            n_col=n_col+1;               
                        end
                    end
                end
            end
        
            n_col=1;
                for t1=1:nb_dt
                    for v=1:nb_van
                        if t1==t && v==br2i 
                            Apec17(n_lig,n_col)=1;
                        end
            n_col=n_col+1;
                    end
                end
 n_lig=n_lig+1;
        end
    end
end

% Contrainte 18 ( contrainte de la capacitÈ de branche i et j, k=p(i,j)>1 )
bc18=zeros(1,1);

Aqc18=zeros(1,nb_prise*nb_dt);

if max(max(p))>1
    n_lig=1;
    for t=1:nb_dt
        for br2ij=1:size(Qcij,1)
           
            % matrice b
            bc18(n_lig,1)=Qcij(br2ij,3);
     
            % matrice sense
%         sensec18(n_lig,1)=-1;
    
            n_col=1;
        
            for t1=1:nb_dt
                for i=1:n
                    for j=1:m(i)        
                        for k=1:p(i,j)             
                            if t1==t && i==Qcij(br2ij,1) && j==Qcij(br2ij,2)
                                Aqc18(n_lig,n_col)=1; 
                            end
            n_col=n_col+1;               
                        end
                    end
                end
            end
    n_lig=n_lig+1;
        end
    end
end

% Contrainte 19 (MO)
bc19=zeros(1,1);

Amoc19=zeros(1,nb_van*nb_dt);

n_lig=1; 
for inter=1:size(t_inter,1)-1
    for t=1:nb_dt
        if t<t_mos(inter) ||(t>t_mof(inter) && t<t_mos(inter+1))||t>t_mof(size(t_inter,1))
            for v=1:nb_van
                    
                % Matrice b
                bc19(n_lig,1)=0;
                % Matrice sense
%                 sensec19(n_lig,1)=0;
                n_col=1;
                for t1=1:nb_dt
                    for v1=1:nb_van
                        if v1==v && t1==t
                            Amoc19(n_lig,n_col)=1;
                        end
                n_col=n_col+1;
                    end                
               end
 n_lig=n_lig+1;
            end
        end
    end
end

% Contrainte 20
bc20=zeros(1,1);

Amoc20=zeros(1,nb_van*nb_dt);

n_lig=1;
for inter=1:size(t_inter,1)    
    for t=1:nb_dt  
        if t_mos(inter)<=t && t<t_mof(inter)                            
            bc20(n_lig,1)=ceil(dt/t_mo_max);
%           sensec20(n_lig,1)=-1;                  
            n_col=1;
            for t1=1:nb_dt
                for v1=1:nb_van
                    if t1>=t && t1<t+ceil(t_mo_max/dt) 
                        Amoc20(n_lig,n_col)=1;
                    end
                n_col=n_col+1;
                end
            end
n_lig=n_lig+1;
        end
    end
end

% SynthËse
%Ac1=[Atsc1 Adc1 sparse(size(Atsc1,1),nb_prise+2*nb_prise*nb_dt+(2*nb_van)*nb_dt)];
Ac2=[sparse(size(Adc2,1),nb_prise) Adc2 sparse(size(Adc2,1),nb_prise+nb_prise*nb_dt+(2*nb_van)*nb_dt) Ayc2];
Ac3=[Atsc3 sparse(size(Atsc3,1),2*nb_prise+nb_prise*nb_dt+(2*nb_van)*nb_dt) Ayc3]; 
Ac4=[Atsc4 Adc4 sparse(size(Atsc4,1),nb_prise+nb_prise*nb_dt+(2*nb_van)*nb_dt) Ayc4];
Ac5=[sparse(size(Aqsc5,1),2*nb_prise) Aqsc5 Aqc5 sparse(size(Aqsc5,1),(2*nb_van)*nb_dt) Ayc5];
Ac6=[sparse(size(Aqsc6,1),2*nb_prise) Aqsc6 Aqc6 sparse(size(Aqsc6,1),(2*nb_van)*nb_dt) Ayc6];
Ac7=[sparse(size(Aqc7,1),3*nb_prise) Aqc7 sparse(size(Aqc7,1),(2*nb_van)*nb_dt) Ayc7];
Ac8=[sparse(size(Aqc8,1),3*nb_prise) Aqc8 Amoc8 sparse(size(Aqc8,1),(nb_van)*nb_dt+nb_prise*nb_dt)];
Ac9=[sparse(size(Apec9,1),3*nb_prise+nb_prise*nb_dt) Amoc9 Apec9 sparse(size(Apec9,1),nb_prise*nb_dt)];
Ac10=[sparse(size(Aqc10,1),3*nb_prise) Aqc10 Amoc10 Apec10 sparse(size(Aqc10,1),nb_prise*nb_dt)];
Ac11=[sparse(size(Aqc11,1),3*nb_prise) Aqc11 Amoc11 Apec11 sparse(size(Aqc11,1),nb_prise*nb_dt)];
Ac12=[sparse(size(Aqc12,1),3*nb_prise) Aqc12 sparse(size(Aqc12,1),(2*nb_van+nb_prise)*nb_dt)];
Ac13=[sparse(size(Aqc13,1),3*nb_prise) Aqc13 sparse(size(Aqc13,1),(2*nb_van+nb_prise)*nb_dt)];
Ac15=[sparse(size(Aqc15,1),3*nb_prise) Aqc15 sparse(size(Aqc15,1),nb_van*nb_dt) Apec15 sparse(size(Aqc15,1),nb_prise*nb_dt)];
Ac16=[sparse(size(Aqc16,1),3*nb_prise) Aqc16 sparse(size(Aqc16,1),(2*nb_van+nb_prise)*nb_dt)];
Ac17=[sparse(size(Aqc17,1),3*nb_prise) Aqc17 sparse(size(Aqc17,1),nb_van*nb_dt) Apec17 sparse(size(Aqc17,1),nb_prise*nb_dt)];
Ac18=[sparse(size(Aqc18,1),3*nb_prise) Aqc18 sparse(size(Aqc18,1),(2*nb_van+nb_prise)*nb_dt)];
Ac19=[sparse(size(Amoc19,1),3*nb_prise+nb_prise*nb_dt) Amoc19 sparse(size(Amoc19,1),(nb_van+nb_prise)*nb_dt)];
Ac20=[sparse(size(Amoc20,1),3*nb_prise+nb_prise*nb_dt) Amoc20 sparse(size(Amoc20,1),(nb_van+nb_prise)*nb_dt)];
clear Adc2 Ayc2 Atsc3 Ayc3 Atsc4 Adc4 Ayc4 Aqsc5 Aqc5 Ayc5 Aqsc6 Aqc6 Ayc6 Aqc7 Ayc7 Aqc8 Amoc8 Amoc9 Apec9 Aqc10 Amoc10 Apec10 
clear Aqc11 Amoc11 Apec11 Aqc12 Aqc13 Aqc15 Apec15 Aqc16 Aqc17 Apec17 Aqc18 Amoc19 Amoc20

Aineq=[Ac3;-1*Ac4;-1*Ac5;Ac6;Ac7;Ac8;Ac9;-1*Ac10;Ac11;Ac20;Ac12;-1*Ac13;Ac15;Ac16;Ac17;Ac18];
bineq=[bc3;-1*bc4;-1*bc5;bc6;bc7;bc8;bc9;-1*bc10;bc11;bc20;bc12;-1*bc13;bc15;bc16;bc17;bc18];
Aeq=[Ac2;Ac19];
beq=[bc2;bc19];
%sense=[sensec2;sensec3;sensec4;sensec5;sensec6;sensec7;sensec19;sensec8;sensec9;sensec10;sensec11;sensec20;sensec12;sensec13;sensec15;sensec16;sensec17;sensec18];
clear Ac2 Ac3 Ac4 Ac5 Ac6 Ac7 Ac19 Ac8 Ac9 Ac10 Ac11 Ac20 Ac12 Ac13 Ac15 Ac16 Ac17 Ac18
clear bc2 bc3 bc4 bc5 bc6 bc7 bc19 bc8 bc9 bc10 bc11 bc20 bc12 bc13 bc15 bc16 bc17 bc18
clear sensec2 sensec3 sensec4 sensec5 sensec6 sensec7 sensec19 sensec8 sensec9 sensec10 sensec11 sensec20 sensec12 sensec13 sensec15 sensec16 sensec17 sensec18

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------- CrÈer la matrice lb, ub et ctype-----------------------% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variables : ts, d et qs
lbts=zeros(1,1);ubts=zeros(1,1);
lbd=zeros(1,1);ubd=zeros(1,1);
lbqs=zeros(1,1);ubqs=zeros(1,1);
% ctypets=zeros(1,1);
% ctyped=zeros(1,1);
% ctypeqs=zeros(1,1);

for prise=1:nb_prise                            
    % limites des variables
    lbts(prise,1)=1;ubts(prise,1)=nb_dt-floor(d_lmin*dd(prise));
    lbd(prise,1)=ceil(d_lmin*dd(prise));ubd(prise,1)=ceil(d_lmax*dd(prise));
    lbqs(prise,1)=floor(q_lmin*qd(prise));ubqs(prise,1)=ceil(q_lmax*qd(prise));
    % type des variables            
    ctypets(prise,1)='I';
    ctyped(prise,1)='I';
    ctypeqs(prise,1)='I';                             
end

% Variables : q et y
lbq=zeros(1,1);ubq=zeros(1,1);
lby=zeros(1,1);uby=zeros(1,1);
% ctypeq=zeros(1,1);
% ctypey=zeros(1,1);

n_lig=1;
for t=1:nb_dt
    for prise=1:nb_prise                           
        % limites des variables
        lbq(n_lig,1)=0;ubq(n_lig,1)=Qrt(t);
        lby(n_lig,1)=0;uby(n_lig,1)=1;
        % type des variables            
        ctypeq(n_lig,1)='I';
        ctypey(n_lig,1)='B';
n_lig=n_lig+1;
    end;            
end;
% Variables : mo & pe
lbmo=zeros(1,1);ubmo=zeros(1,1);
lbpe=zeros(1,1);ubpe=zeros(1,1);
% ctypemo=zeros(1,1);
% ctypepe=zeros(1,1);

n_lig=1;
for t=1:nb_dt
    for v=1:nb_van
        % limite 
        lbmo(n_lig,1)=0;ubmo(n_lig,1)=1; 
        lbpe(n_lig,1)=0;ubpe(n_lig,1)=Qrt(t);
        % type
        ctypemo(n_lig,1)='B';
        ctypepe(n_lig,1)='I';
n_lig=n_lig+1;
    end;
end;

% SynthËse
lb=[lbts;lbd;lbqs;lbq;lbmo;lbpe;lby];
ub=[ubts;ubd;ubqs;ubq;ubmo;ubpe;uby];
ctype=[ctypets;ctyped;ctypeqs;ctypeq;ctypemo;ctypepe;ctypey];
clear ctypets ctyped ctypeqs ctypeq ctypemo ctypepe ctypey
clear lbts lbd lbqs lbq lbmo lbpe lby
clear ubts ubd ubqs ubq ubmo ubpe uby


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------- Cplex solver Optimization -------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = cplexoptimset('cplex');
%options.Exportmodel = 'debug.lp';
options.diagnostics = 'on';

% solver
[x,obj,exitflag,output]=cplexmiqp(H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);

% Affichage

output
disp('X Optimum=')
x
disp('CritËre Optimum=')
obj+sum(const)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------Extraction et vÈrification-----------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VÈrification de la valeur de fonction d'objectif

obj1=f'*x+1/2*x'*Q*x+sum(const);

% Extraction des valeurs x aux variables correspondant  

ts=zeros(1,1); ds=zeros(1,1);qs=zeros(1,1);
OBJ1=zeros(1,1);OBJ2=zeros(1,1);OBJ3=zeros(1,1);OBJ4=zeros(1,1);

for prise=1:nb_prise
    ts(prise)=x(prise);
    ds(prise)=x(prise+nb_prise);
    qs(prise)=x(prise+2*nb_prise);
    % Valeur d'objectif de diffÈrents terms 
    OBJ1(prise)=alpha(prise)*(td(prise)-ts(prise))^2;
    OBJ2(prise)=beta(prise)*(dd(prise)-ds(prise))^2;
    OBJ3(prise)=mphi(prise)*(qd(prise)-qs(prise))^2;
    OBJ4(prise)=cphi(prise)*(qdv(prise)-ds(prise)*qs(prise));
end

q=zeros(1,1);y=zeros(1,1);

n_lig=1;
for t=1:nb_dt
    for prise=1:nb_prise
        q(prise,t)=x(nb_prise*3+n_lig);
        y(prise,t)=x(nb_prise*3+(nb_prise+2*nb_van)*nb_dt+n_lig);
n_lig=n_lig+1;
    end
end

mo=zeros(1,1);pe=zeros(1,1);
OBJ5=zeros(1,1);OBJ6=zeros(1,1);

n_lig=1;
for t=1:nb_dt
    for v=1:nb_van
        mo(v,t)=x(n_lig+3*nb_prise+nb_prise*nb_dt);
        pe(v,t)=x(n_lig+3*nb_prise+nb_prise*nb_dt+nb_van*nb_dt);
        OBJ5(n_lig)=w(5)/sum(mo_max)*mo(n_lig);
        OBJ6(n_lig)=w(6)/sum(Qrt)*pe(n_lig);
    n_lig=n_lig+1;
    end
end

% VÈrification de la valeur de fonction d'objectif ? partir de fonction
% d'origine
OBJ=sum(OBJ1)+sum(OBJ2)+sum(OBJ3)+sum(OBJ4)+sum(OBJ5)+sum(OBJ6);OBJ3)+sum(OBJ4)+s