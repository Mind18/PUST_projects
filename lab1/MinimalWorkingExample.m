 


function MinimalWorkingExample()
    addpath ('D:\SerialCommunication') ; % add a path
    initSerialControl COM3 % initialise com port

    while (1)
    %% obtaining measurements
    measurements = readMeasurements (1:7) ; % read measurements

    %% processing of the measurements
    disp ( measurements ) ; % process measurements

    %% sending new values of control signals
    sendControls ([ 1 , 2 , 3 , 4 , 5 , 6] , [ 0 , 0 , 0 , 0 , 50 , 50]) ;

    %% synchronising with the control process
    waitForNewIteration () ; % wait for new iteration
    end
end