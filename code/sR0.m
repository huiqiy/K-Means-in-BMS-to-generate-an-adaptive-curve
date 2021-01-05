function sR0(block)
% Level-2 MATLAB file S-Function for inherited sample time demo.
%   Copyright 1990-2009 The MathWorks, Inc.

  setup(block);
  
%endfunction

function setup(block)
  
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 0;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = false;
  block.InputPort(1).SamplingMode = 'Sample';
 
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = false;
  block.InputPort(2).SamplingMode = 'Sample';
  
  block.NumDialogPrms     = 2;
  block.DialogPrmsTunable = {'Tunable','Tunable'};

  %% Set block sample time to inherited
  block.SampleTimes = [10 0];
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Update',                  @Update);  

  
%endfunction

function DoPostPropSetup(block)

  %% Setup Dwork
  block.NumDworks = 2;
  
  block.Dwork(1).Name = 'x1'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;
  
  block.Dwork(2).Name = 'x2'; 
  block.Dwork(2).Dimensions      = 1;
  block.Dwork(2).DatatypeID      = 0;
  block.Dwork(2).Complexity      = 'Real';
  block.Dwork(2).UsedAsDiscState = true;

%endfunction

function InitConditions(block)

  %% Initialize Dwork

  block.Dwork(1).Data = -1;
  block.Dwork(2).Data = -1;
%endfunction

function Output(block)
  persistent i SOC Param;
  if block.Dwork(1).Data==-1
      i=1;
      SOC=[];
      Param=[];
  else
      i=i+1;
      SOC=[SOC;block.Dwork(1).Data];
      Param=[Param;block.Dwork(2).Data];
  end
  X=[SOC Param];
  if isempty(X)~=1
  h=1/20;%Aufloesung
  if size(SOC)<1/h
      N=1;
      center=kmeans1( X,N );%nur ein zentral Punkt
  else
      N=ceil((1-min(SOC))/h);%die Zahl von zentralen Punkten
      [idx,center,sumD,~]=kmeans3(X, N); % mehr als ein zentral Punkt
  end
  [~, I] = sort(center(:,1));
  center = center(I,:);  
  figure(1);
  cla;
  plot(SOC,Param,'k*');
  hold on;
  plot(center(:,1),center(:,2),'r+','MarkerSize',15,'LineWidth',2);
  end
  
  
%endfunction

function Update(block)

%   block.Dwork(1).Data = block.InputPort(1).Data;
%   block.Dwork(2).Data = block.InputPort(2).Data;
  block.Dwork(1).Data = block.InputPort(1).Data;
  block.Dwork(2).Data = block.InputPort(2).Data;
%   figure(2);
%   hold on;
%   plot(block.Dwork(1).Data,block.Dwork(2).Data,'r+');
%endfunction



