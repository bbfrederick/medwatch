function y = fcn()    
    persistent b;    %automatically initialized to []
    persistent hr;   %automatically initialized to []
    
    %% extrinsic declarations
    coder.extrinsic('ble');
    coder.extrinsic('characteristic');
    coder.extrinsic('read');
    if(isempty(b))
        b = ble("MED-WATCH009");
        b.Characteristics
        hr = characteristic(b, "6E400001-B5A3-F393-E0A9-E50E24DCCA9E", "6E400002-B5A3-F393-E0A9-E50E24DCCA9E");
    end
    
    %% initialize output
%     y = zeros(1,50);
    %% read data from BLE device
    while(1)
%     data = read(hr,'latest')
    data = read(hr)
    end
    %% post-process the data
%     y = data
end   