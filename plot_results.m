% SangMok Han, Turbulence Lab, Yonsei University, August 2023
function plot_results(x,y,u,v,psi,omega,p)
    % Stream function
    figure;
    contourf(x,y,psi',15); colorbar;
    xlabel('x'); ylabel('y'); zlabel('\psi'); title('Stream Function'); set(gca,'fontsize',12);
%     contourf(x(round(end*2/3):end),y(1:round(end/3)),psi(round(end*2/3):end,1:round(end/3))',200); % Q2_4

    % Vorticity
    figure;
    contourf(x,y,omega',50); colorbar;
    xlabel('x'); ylabel('y'); zlabel('\phi'); title('Vorticity'); set(gca,'fontsize',12);
    
    % Pressure
    figure;
    contourf(x,y,p',30); colorbar;
    xlabel('x_p'); ylabel('y_p'); zlabel('pressure'); title('Pressure'); set(gca,'fontsize',12);

    % Centerline velocity u
    figure;
    plot(u(end/2,:),y,'Linewidth',3);
    xlabel('u/U'); ylabel('y/L'); title('Centerline Velocity u'); set(gca,'fontsize',12);

    % Centerline velocity v
    figure;
    plot(x,v(:,end/2),'Linewidth',3);
    xlabel('x/L'); ylabel('v/U'); title('Centerline Velocity v'); set(gca,'fontsize',12);
end