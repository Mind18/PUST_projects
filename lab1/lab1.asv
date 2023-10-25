import MinimalWorkingExample.*

%clear;

W1 = 50;
G1 = 28;    %15, 35, 50

figure
k =1;

while (1)
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM3 % initialise com port
    %% obtaining measurements
    measurements = readMeasurements (1) ; % read measurements
    y(k) = measurements;
    x(k) = k;
    %% processing of the measurements
    disp ("T1 " +  measurements ) ; % process measurements

    %% sending new values of control signals
    sendControls ([ 1 , 2, 3, 4, 5, 6] , [W1, 0 , 0 , 0 , G1, 0]) ;
    
    %% synchronising with the control process
    plot(x,y)
    drawnow;
    
    waitForNewIteration () ; % wait for new iteration
    
    k = k+1;
end