function [] = vis_funct(x,mu,R,t,i,j,N,M0)
    close all;
    G = 6.67259e-8;
    x = squeeze(x(i,j,:,:,:));
    t = squeeze(t(i,j,:));
    M1 = mu(i)*M0;
    M2 = (1-mu(i))*0.5*M0;
    M3 = M2;
    R = R(j);
    omega = sqrt(G*(M1+M2+M3)/R^3);
    T = 2*pi/omega;
    
    
    figure;
    hold all;
    
    plot3(x(:,1,1),x(:,2,1),x(:,3,1),'linewidt',1);
    plot3(x(:,1,2),x(:,2,2),x(:,3,2),'linewidt',1);
    plot3(x(:,1,3),x(:,2,3),x(:,3,3),'linewidt',1);
    xlabel('X','color','blue','rotation',0);
    ylabel('Y','color','blue','rotation',0);
    %title(['Stability = ',num2str(N_stab(mu(i),M0)),'; R = ',num2str(R,'%.1e'),'(cm)'],'color','blue');
    title(['M1 = ',num2str(M1,'%.1e'),'(g); M2 = ',num2str(M2,'%.1e'),...
        '(g); M3 = ',num2str(M3,'%.1e'),'(g)'],'color','blue');
    %title(['Stability = ',num2str(stab),'; R = 1e7(cm)'],'color','blue');
    %legend('M1 = 2*M0','M1 = 0.03*M0','M1 = 0.03*M0');
    set(gca,'Fontsize',14);
    grid on;

    x12_t = zeros(1,length(t));
    x23_t = zeros(1,length(t));
    x13_t = zeros(1,length(t));

    for k = 1:length(t)
       x12_t(k) = norm(squeeze(x(k,:,1))-squeeze(x(k,:,2)))/R;
       x23_t(k) = norm(squeeze(x(k,:,2))-squeeze(x(k,:,3)))/R;
       x13_t(k) = norm(squeeze(x(k,:,1))-squeeze(x(k,:,3)))/R;
    end

    figure;
    hold all;
    plot(t/T,x12_t,'linewidt',1);
    plot(t/T,x23_t,'linewidt',1);
    plot(t/T,x13_t,'linewidt',1);
    xlabel('t/T','color','blue','rotation',0);
    ylabel('r/R','color','blue','rotation',0);
    title(['PN factor = ',num2str(6.67259e-8*M0/(R*2.998e10*2.998e10)),...
        '; R = ',num2str(R,'%.1e'),'(cm)'],'color','blue');
    %title(['Stability = ',num2str(stab),'; R = 1e7(cm)'],'color','blue');
    legend('r12(t)/R','r23(t)/R','r13(t)/R');
    set(gca,'Fontsize',14);
    xlim([0 N]);
    grid on;

    figure;
    hold all;
    plot(t/T,log10(abs(x12_t)),'linewidt',1);
    plot(t/T,log10(abs(x23_t)),'linewidt',1);
    plot(t/T,log10(abs(x13_t)),'linewidt',1);
    xlabel('t/T','color','blue','rotation',0);
    ylabel('lg(|r/R|)','color','blue','rotation',0);
    title(['Stability = ',num2str(N_stab(mu(i),M0))],'color','blue');
    legend('lg(|r12(t)/R|)','lg(|r23(t)/R|)','lg(|r13(t)/R|)');
    set(gca,'Fontsize',14);
    xlim([0 N]);
    grid on;
end
function stab = N_stab(mu,M0)
    M1 = mu*M0;
    M2 = (1-mu)*M0*0.5;
    M3 = M2;
    stab = ((M1*M2+M1*M3+M2*M3)*27/(M1+M2+M3)^2) - 1;
end
