clear;
M0 = 1.989e33;
R = [3e6];
N_grav = 3;
N = 50;
N_per = 10;
t_step = 1e-7;
tol_err = 5e-13;
stab_min = -0.95;
stab_max = 0.95;
PN_param = 6.67259e-8*M0./(R*2.998e10*2.998e10);

fndmu_max = @(mu)N_stab(mu,M0,stab_max);
mu_max = fzero(fndmu_max,0.5);
fndmu_min = @(mu)N_stab(mu,M0,stab_min);
mu_min = fzero(fndmu_min,0.5);
mu = linspace(mu_min,mu_max,25);
stab = zeros(length(mu),length(R));
stab_r = zeros(length(mu),length(R));
crit = zeros(length(mu),length(R));
x = zeros(length(mu),length(R),N*N_per,N_grav,3);
v = zeros(length(mu),length(R),N*N_per,N_grav,3);
t = zeros(length(mu),length(R),N*N_per);

for i = 1:length(mu)
    M1 = mu(i)*M0;
    M2 = (1-mu(i))*0.5*M0;
    M3 = M2;
    for j = 1:length(R)
        %^[x(i,j,:,:,:),v(i,j,:,:,:),t(i,j,:),stab(i,j),stab_r(i,j),crit(i,j)] = ...
            %test_f(M1,M2,M3,R(j));
        [x_aux,v_aux,t_aux,stab(i,j),stab_r(i,j),crit(i,j)] = ...
            test_f(M1,M2,M3,R(j),N_grav,N,N_per,t_step,tol_err);
        [x,v,t] = pass_f(x,v,t,x_aux,v_aux,t_aux,i,j,length(squeeze(x_aux(:,1,1))));
        aux_str = ['init',num2str(i),'.dat'];
        copyfile('init.dat',aux_str);
        disp('***Porgress***');
        disp((i+j)/(length(mu)+length(R)));
        disp('**************');
    end
end

function [x,v,t] = pass_f(x,v,t,x_aux,v_aux,t_aux,i,j,leng)
    leng_curr = length(squeeze(x(1,1,:,1,1)));
    if leng >= leng_curr
        x(i,j,:,:,:) = x_aux(1:leng_curr,:,:);
        v(i,j,:,:,:) = v_aux(1:leng_curr,:,:);
        t(i,j,:) = t_aux(1:leng_curr);
    else
        x = x(:,:,1:leng,:,:);
        v = v(:,:,1:leng,:,:);
        t = t(:,:,1:leng);

        x(i,j,:,:,:) = x_aux;
        v(i,j,:,:,:) = v_aux;
        t(i,j,:) = t_aux;
    end
end

function stab = N_stab(mu,M0,rez)
    M1 = mu*M0;
    M2 = (1-mu)*M0*0.5;
    M3 = M2;
    stab = ((M1*M2+M1*M3+M2*M3)*27/(M1+M2+M3)^2)-1-rez;
end

function [x,v,t,stab,stab_r,crit] = test_f(M1,M2,M3,R,N_grav,N,N_per,t_step,tol_err)
    G = 6.67259e-8;
    c = 2.998e10;

    S1 = [0 0 0];
    S2 = [0 0 0];
    S3 = [0 0 0];

    omega = sqrt(G*(M1+M2+M3)/R^3);
    omega_v = [0 0 omega];
    T = 2*pi/omega;
    t_end = N*T;
    t_wr = T/N_per;

    [x1,x2,x3,stab] = init_N(M1,M2,M3,R);
    [x1,x2,x3,v1,v2,v3,crit] = triangl_adj_PN(x1,x2,x3,S1,S2,S3,...
        M1,M2,M3,omega_v,R);
    crit = crit/(M1+M2+M3);
    create_init(N_grav,t_end,t_wr,t_step,tol_err,M1,M2,M3,...
            x1,x2,x3,v1,v2,v3,S1,S2,S3);
    system('a.exe');
    [x,v,stab_r,t] = numerical_simple('output_chks_.dat');
    function [] = create_init(N_grav,t_end,t_wr,t_step,tol_err,M1,M2,M3,...
            x1,x2,x3,v1,v2,v3,S1,S2,S3)
        if exist('output_chks_.dat','file')
            delete output_chks_.dat;
        end

        if exist('init.dat','file')
            delete init.dat;
        end
        %Init prepare
        fid=fopen('init.dat','wt');

        wr_init(fid,'N_grav',N_grav,'int');

        wr_init(fid,'t_end',t_end,'real');
        wr_init(fid,'t_wr',t_wr,'real');
        wr_init(fid,'t_step',t_step,'real');
        wr_init(fid,'tol_err',tol_err,'real');

        wr_init(fid,'M1',M1,'real');
        wr_init(fid,'M2',M2,'real');
        wr_init(fid,'M3',M3,'real');

        wr_init(fid,'S1(1)',S1(1),'real');
        wr_init(fid,'S1(2)',S1(2),'real');
        wr_init(fid,'S1(3)',S1(3),'real');
        wr_init(fid,'S2(1)',S2(1),'real');
        wr_init(fid,'S2(2)',S2(2),'real');
        wr_init(fid,'S2(3)',S2(3),'real');
        wr_init(fid,'S3(1)',S3(1),'real');
        wr_init(fid,'S3(2)',S3(2),'real');
        wr_init(fid,'S3(3)',S3(3),'real');


        wr_init(fid,'x1(1)',x1(1),'real');
        wr_init(fid,'x1(2)',x1(2),'real');
        wr_init(fid,'x1(3)',x1(3),'real');
        wr_init(fid,'v1(1)',v1(1),'real');
        wr_init(fid,'v1(2)',v1(2),'real');
        wr_init(fid,'v1(3)',v1(3),'real');
    
        wr_init(fid,'x2(1)',x2(1),'real');
        wr_init(fid,'x2(2)',x2(2),'real');
        wr_init(fid,'x2(3)',x2(3),'real');
        wr_init(fid,'v2(1)',v2(1),'real');
        wr_init(fid,'v2(2)',v2(2),'real');
        wr_init(fid,'v2(3)',v2(3),'real');

        wr_init(fid,'x3(1)',x3(1),'real');
        wr_init(fid,'x3(2)',x3(2),'real');
        wr_init(fid,'x3(3)',x3(3),'real');
        wr_init(fid,'v3(1)',v3(1),'real');
        wr_init(fid,'v3(2)',v3(2),'real');
        wr_init(fid,'v3(3)',v3(3),'real');

        fclose(fid);
        function wr_init(fid,varn,varv,type)
            fprintf(fid,[type,'(8), parameter :: ',varn,' = ',...
                num2str(varv),'\n']);
        end
    end

    function [x1,x2,x3,stab] = init_N(M1,M2,M3,R)
        r1 = sqrt(M2^2+M3^2+M2*M3)*R/(M1+M2+M3);
        r2 = sqrt(M1^2+M3^2+M1*M3)*R/(M1+M2+M3);
        r3 = sqrt(M1^2+M2^2+M1*M2)*R/(M1+M2+M3);

        x1(1) = r1;
        x1(2) = 0;
        x1(3) = 0;

        cos_a1 = (r1^2+r2^2-R^2)/(2*r1*r2);
        x2(1) = cos_a1*r2;
        x2(2) = sqrt(1-cos_a1^2)*r2;
        x2(3) = 0;

        cos_a2 = (r1^2+r3^2-R^2)/(2*r1*r3);
        x3(1) = cos_a2*r3;
        x3(2) = -sqrt(1-cos_a2^2)*r3;
        x3(3) = 0;

        stab = ((M1*M2+M1*M3+M2*M3)*27/(M1+M2+M3)^2) - 1;
    end

    function [x1,x2,x3,v1,v2,v3,crit] = triangl_adj_PN(x1,x2,x3,S1,S2,S3,...
        M1,M2,M3,omega_v,R)
        a = R;
        fndngl = @(phi)find_angl(phi,x1,x3);
        angl = fzero(fndngl,0);
        M = rot_m(angl);
        x1 = x1*M;
        x2 = x2*M;
        x3 = x3*M;

        if (x2(2)<x1(2))
            M = rot_m(pi);
            x1 = x1*M;
            x2 = x2*M;
            x3 = x3*M;
        end

        if (x1(1) < x3(1))
            error('x1 and x3 are on the wrong side');
        end
        v1 = cross(omega_v,x1);
        v2 = cross(omega_v,x2);
        v3 = cross(omega_v,x3);

        Y = r_h_pn_F(vertcat(x1,x2,x3),vertcat(v1,v2,v3),[M1,M2,M3],vertcat(S1,S2,S3));
        Y1 = Y(1,:);
        Y2 = Y(2,:);
        Y3 = Y(3,:);
        n21 = (x2-x1)/norm(x2-x1);
        n31 = (x3-x1)/norm(x3-x1);
        n23 = (x2-x3)/norm(x2-x3);

        delta21 = 2*dot((M1*Y1+M3*Y3-(M1+M3)*Y2),(n21+n31))/(9*M1*dot(omega_v,omega_v));
        delta23 = 2*dot((M1*Y1+M3*Y3-(M1+M3)*Y2),(n23-n31))/(9*M3*dot(omega_v,omega_v));
        delta13 = 2*dot((M2*Y2+M3*Y3-(M2+M3)*Y1),(-n31+n23))/(9*M3*dot(omega_v,omega_v));

        a1 = R + delta21;
        a2 = R + delta23;
        a3 = R + delta13;

        Y_t = -(M1*Y1+M2*Y2+M3*Y3)/(dot(omega_v,omega_v));

        S = h_eron(a,a,a);
        S1 = h_eron(a1,a2,a3);
        h = 2*S/a;
        h1 = 2*S1/a3;
        dh = h1-h;

        xsi1(2) = Y_t(2)/(M1+M2+M3)-M2*dh/(M1+M2+M3);
        xsi3(2) = xsi1(2);
        xsi2(2) = dh + xsi1(2);

        a1y = x1(2)+xsi1(2)-x2(2)-xsi2(2);
        a1x = sqrt(a1^2 - a1y^2);

        c_2x = x1(1)-x2(1)-a1x;
        c_3x = a-a3;

        xsi1(1) = (Y_t(1)-M2*c_2x-M3*c_3x)/(M1+M2+M3);
        xsi2(1) = xsi1(1)+c_2x;
        xsi3(1) = xsi1(1)+c_3x;

        xsi1(3) = 0;
        xsi2(3) = 0;
        xsi3(3) = 0;
        crit = norm(M1*xsi1+M2*xsi2+M3*xsi3-Y_t);
        x1 = x1+xsi1;
        x2 = x2+xsi2;
        x3 = x3+xsi3;
        v1 = cross(omega_v,x1);
        v2 = cross(omega_v,x2);
        v3 = cross(omega_v,x3);


%         plot([x1(1),x2(1),x3(1),x1(1)],[x1(2),x2(2),x3(2),x1(2)]);
%         hold all;

%         plot([x1(1),x2(1),x3(1),x1(1)],[x1(2),x2(2),x3(2),x1(2)]);

        function S = h_eron(a1,a2,a3)
            p = (a1+a2+a3)*0.5;
            S = sqrt(p*(p-a1)*(p-a2)*(p-a3));
        end

        function diff = find_angl(phi,x,y)
            U = rot_m(phi);
            x_t = x*U;
            y_t = y*U;
            diff = x_t(2)-y_t(2);
        end

        function M = rot_m(phi)
            M = zeros(3,3);
            M(1,1) = cos(phi);
            M(1,2) = -sin(phi);
            M(2,1) = sin(phi);
            M(2,2) = cos(phi);
            M(3,3) = 1;
        end
    end

    function [x,v,stab_r,t,t_step] = numerical_simple(dumpname)
        N_b = 3;
        data = load(dumpname);
        t = data(:,1);
        t_step = data(:,2);
        data(:,1) = [];
        data(:,1) = [];
        x = zeros(length(t),3,N_b);
        v = zeros(length(t),3,N_b);
        for i = 1:N_b
            x(:,:,i) = data(:,(i-1)*3+1:i*3);
            v(:,:,i) = data(:,N_b*3+(i-1)*3+1:N_b*3+i*3);
        end
        stab_r = stab_check(x);

        function rez = stab_check(x)
%             rez = 0;
%             R_aux = norm(squeeze(x(1,:,1))-squeeze(x(1,:,2))); 
%             for k = 1:length(x(:,1,1))
%                 stab12 = norm(squeeze(x(k,:,1))-squeeze(x(k,:,2)))/R_aux;
%                 stab23 = norm(squeeze(x(k,:,2))-squeeze(x(k,:,3)))/R_aux;
%                 stab31 = norm(squeeze(x(k,:,3))-squeeze(x(k,:,1)))/R_aux;
%                 rez = max([rez,stab12,stab23,stab31]);
%             end

            stab_t(1) = stab_couple(squeeze(x(1,:,1)),squeeze(x(end,:,1)),...
                squeeze(x(1,:,2)),squeeze(x(end,:,2)));
            stab_t(2) = stab_couple(squeeze(x(1,:,2)),squeeze(x(end,:,2)),...
                squeeze(x(1,:,3)),squeeze(x(end,:,3)));
            stab_t(3) = stab_couple(squeeze(x(1,:,1)),squeeze(x(end,:,1)),...
                squeeze(x(1,:,3)),squeeze(x(end,:,3)));
            rez = max(stab_t);

            function stab = stab_couple(x1_0,x1_end,x2_0,x2_end)
                x12_0 = norm(x1_0-x2_0);
                x12_end = norm(x1_end-x2_end);
                stab = abs(x12_end - x12_0)/max(x12_0,x12_end);
            end
        end

    end
end