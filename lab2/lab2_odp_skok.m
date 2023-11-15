import sendControlsToG1AndDisturbance.*

%clear;

W1 = 50;
G1 = 28;    % Pkt. pracy 32.62
Z = 0;      %10, 30, 50

figure
k =1;

Z=0;

while (1)
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM9 % initialise com port
    %% obtaining measurements
    measurements = readMeasurements (1) ; % read measurements
    y(k) = measurements;
    x(k) = k;
    %% processing of the measurements
    disp ("T1 " +  measurements ) ; % process measurements

    %% sending new values of control signals
    sendControls ([ 1 , 2, 3, 4, 5, 6] , [W1, 0 , 0 , 0 , G1, 0]) ;
    sendControlsToG1AndDisturbance(G1, Z);
    
    if k==30
        Z = 50;
    end
    z(k) = Z;
    
    %% synchronising with the control process
    plot(x,y, 'Color', '#0072BD');
    hold on;
    plot(x, z, 'Color', '#D95319');
    xlabel('k');
    ylabel('y(k) | z(k)');
    legend('y(k)', 'z(k)', 'Location', 'best');
    drawnow;
    
    waitForNewIteration () ; % wait for new iteration
    
    k = k+1;
end