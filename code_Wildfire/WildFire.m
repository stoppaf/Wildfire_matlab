clear all
close all
clc

% Modeling and Simulating Wildfire with MATLAB

% WILDFIRE PROPAGATION - TI - VALLE MAGGIA - GORDEVIO

% http://www.soms.ethz.ch/matlab
% Author: Anastasia Gavrilova, Daniele Speziali, Francesco Stoppa, Mattia Bacchetta-Cattori
% Semester: HS 2013

%Anhang: 'plusmat.m', 'Terreno.mat', 'Densita.mat', 'lista_case.mat'
%(this MFile and Matrix file are request to the program)

%% 0. Declaration
% 0.0. Figure
schermoint = get(0,'ScreenSize'); 
fig=figure('Position',schermoint);


% 0.1. Display - Questions -  parameters
 Questions = {'Wind direction:','Wind intensity:','Details:','Fire Starting point - x:','Fire Starting point - y:'};
   titolo = 'Parameters';
   num_lines= 1;
   def     = {'N, NE, E, SE, S, SW, W, NW','0 - 1','2 (min) - 4 (max)','0 - 1','0 - 1'};
   answ  = inputdlg(Questions,titolo,num_lines,def);
   
   if answ{1}=='N'
       wd=[1 0];
   elseif answ{1}=='S'
       wd=[-1 0];
   elseif answ{1}=='E'
       wd=[0 1];
   elseif answ{1}=='W'
       wd=[0 -1];
   elseif answ{1}=='NE'
       wd=[1 1];
   elseif answ{1}=='SW'
       wd=[-1 -1];
   elseif answ{1}=='SE'
       wd=[-1 1];
   elseif answ{1}=='NW'
       wd=[1 -1];
   end
   wi = str2num(answ{2});%2  
   dur=0.4;
   
   % Number of steps
   St=200;
   
   detail = str2num(answ{3});%5
   y0 = str2num(answ{4});%6
   x0 = str2num(answ{5});%7
  
   % Here can be changed the output
   ID=2;


% 0.2. Parameters
WindDir=wd;                 %Wind direction: is a vector with -1,0,1 as component;
WindInt=wi;                 %Wind intensity: from 0 to 1;
Propag=0.1;                 %Propagation rate of the fire: from 0 to 1;
Duration=dur;               %Burn-time of the trees: from 0 to 1 (0:need 1 step to burn; 1:burns for ever);
Steps=St;                   %Steps of the simulation;

SlopeInt=1;                 %Slope gradient used in chapter 2.2: 'im Normalfall' SlopeInt=1;
                            %   The size of the gradient is from 0 to 1, this value is moltiplicated by SlopeInt.
CorrFactorWind=0.2;         %Correction Factor for Wind, used in chapter 3.2.1;
CorrFactorSlope=0.2;        %Correction Factor for Slope, used in chapter 3.2.1;
                            %   CorrFactor allow to change the influence of Wind and Slope on the normal propagation 3.2.1;

CorrSlopeToWind=0.2;        %Influence of the slope on the wind-additional-propagation 3.2.2, 3.2.3;
CorrWindToSlope=0.2;        %Influence of the wind on the slope-additonal-propagation 3.2.4, 3.2.5;

%% 1. Create the forest
%Size of the forest - N
if detail == 0;
N=31;    
elseif detail == 1;
N=61;
elseif detail == 2;
N=121; 
elseif detail == 3;
N=241;
elseif detail == 4;
N=481;
end
  
map=zeros(N,N,2);           %Forest: if map(i,j)==0 the tree doesn't burn (green);
Color=ones(N,N,3);          %        if map(i,j)=={1,2,3,4,5,6} the tree is burning with different intensity(jellow to red);
                            %        if map(i,j)==7 the tree is already burned (black).
                            %Color is a 3D-Matrix, that stores for each cell 3 values. They represent the RGB color of the cells.
fire1=zeros(1,2);           %Localize the cells with active fire: fire(n,1)=i, fire(n,2)=j;
f=1;                        %The size of fire1 (+1);
burned=zeros(1,2);          %Locatize the cells with burned tree;
b1=1;                       %The size of burned (+1);
time=zeros(1);              %
tm=1;                       %

x0=round(x0*N);
y0=round(y0*N);

% 1.1. Stating point of the Fire (grid 4x4)                           
for i=0:4
    for j=0:4
        map(x0+i,y0+j,1)=1;
        fire1(f,1)=x0+i;
        fire1(f,2)=y0+j;
        f=f+1;
    end
end

%% 2. Create the terrain model
Surface=zeros(N,N,5);           %Surface(:,:,1) store the altitude of each cell;
                                %Surface(:,:,{2,3}) is a 2X1 vector that store the inclination-vectort of each cell;
                                %Surface(:,:,4) store the modulo of inclination-vector;
                                %Surface(:,:,5) store the density of the forest.
                                

% 2.1. Import: VALLE MAGGIA GEBIET
load('Terreno.mat','Terreno');
load('Densita.mat','Densita');
piano = Terreno;
density = Densita;

% 2.1.1. Enlarge the imput matrix (function plusmat)
if detail == 1;
Surface(:,:,1)=plusmat(piano);
Surface(:,:,5)=plusmat(density);
load('lista_case.mat','List1');
house=List1;
elseif detail == 2;
Surface(:,:,1)=plusmat(plusmat(piano));
Surface(:,:,5)=plusmat(plusmat(density));
load('lista_case.mat','List2');
house=List2;
elseif detail == 3;
Surface(:,:,1)=plusmat(plusmat(plusmat(piano)));
Surface(:,:,5)=plusmat(plusmat(plusmat(density)));
load('lista_case.mat','List3');
house=List3;
elseif detail == 4;
Surface(:,:,1)=plusmat(plusmat(plusmat(plusmat(piano))));
Surface(:,:,5)=plusmat(plusmat(plusmat(plusmat(density))));
load('lista_case.mat','List4');
house=List4;
end

% 2.2. Computation of the inclination of the Surface
xy=[1 0;0 1;-1 0;0 -1;1 1;-1 1;-1 -1;1 -1];                 %Neighbour: help-number to search the inclination-vector;
for i=2:N-1
    for j=2:N-1
        
        diff=0;
        for k=1:4                                           %k go through the xy-matrix;
            
            %If the new calculated inclination is bigger than the old one, than store the new one.
            if diff<(Surface(i+xy(k,1),j+xy(k,2),1)-Surface(i,j,1))
                diff=Surface(i+xy(k,1),j+xy(k,2),1)-Surface(i,j,1);            	%Store the modulo of the inclination-vector;
                kmax=k;                                                         %Store which vector of xy;
            end
            %Same as above, but the distance between two diagonal neighbour is no more 1, it is sqrt(2).
            if diff<((Surface(i+xy(k+4,1),j+xy(k+4,2),1)-Surface(i,j,1))/sqrt(2))
                diff=(Surface(i+xy(k+4,1),j+xy(k+4,2),1)-Surface(i,j,1))/sqrt(2);
                kmax=k+4;
            end
        end
        Surface(i,j,2)=xy(kmax,1);                          %Saving the first component of the inclination-vector;
        Surface(i,j,3)=xy(kmax,2);                          %Saving the second component of the inclination-vector;
        Surface(i,j,4)=diff;                                %Saving the modulo of the inclination-vector;
    end
end

maxDiff=max(max(Surface(:,:,4)));                           %Search the biggest inclination-vector;
for i=2:N-1
    for j=2:N-1
        Surface(i,j,4)=(Surface(i,j,4)/maxDiff)*SlopeInt;   %Normalization of inclination-vector: inclination-vector is now from (0 to 1)/SlopeInt;
    end
end

for i=1:size(house,1)                                       %Insert home in matrix Surface;
    map(house(i,1),house(i,2),1)=8;
end


% 2.3. Color of the map (see 3.3.7.)
MaxAltitude=max(max(Surface(:,:,1)));
MinAltitude=min(min(Surface(:,:,1)));
for i=1:N
  	for j=1:N
      	switch map(i,j,1)
         	case 0                      %If 0: not burned;
                hohe=(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude);        %Hohenkurven
            	hohe=floor(hohe*100*detail)/(100*detail);
                if hohe==0.1
                    col=0;
                elseif hohe==0.2
                    col=0;
                elseif hohe==0.3
                    col=0;
                elseif hohe==0.4
                    col=0;
                elseif hohe==0.5
                    col=0;
                elseif hohe==0.6
                    col=0;
            	elseif hohe==0.7
                    col=0;
                elseif hohe==0.8
                    col=0;
                elseif hohe==0.9
                    col=0;
                else 
                    col=1;
                end

                Color(i,j,1)=0.6*(1-col);           %Color: green(if col==0), brown(if col==1)
               	Color(i,j,2)=(0.9*(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude)*col) *(1-0.3*Surface(i,j,5)) + 0.3*(1-col);      %Low altitude has darker green than high altitude;
              	Color(i,j,3)=0;
          	case {1,2}                              %If 1 or 2: starting burning;
             	Color(i,j,1)=1;                     %Color: red;
                Color(i,j,2)=0;
                Color(i,j,3)=0;
            case 8
                Color(i,j,1)=0;
                Color(i,j,2)=0;
                Color(i,j,3)=1;
        end
    end
end


%% 3. Simulation


% 3.1. Preparation for simulating
xy1=[-1 -1;-1 0;-1 1;0 1;1 1;1 0;1 -1;0 -1];    %Neighbours;
xyWind2=[2 1;2 2; 1 2];                         %Additional neighbours in case of diagonal wind;
xyWind3=[2 0;3 0];                              %Additional neighbours in case of straight wind;

map(:,:,2)=map(:,:,1);                          %Next matrix: means the forest situation in the next step;
fire2=fire1;                                    %Next matrix;
b2=b1;                                          %Next number of burned tree;

% 3.2. Preparation movie
      aviobj = avifile('wildfire.avi');

% 3.3. Start the simulation
for t=1:Steps              	%Steps;
    
    map(:,:,1)=map(:,:,2);  %Next-Matrix becomes the Current-Matrix;
    
    fire1=fire2;            %Next-Matrix becomes the Current-Matrix
    fire2=zeros(1,2);       %Clear the matrix fire2
    f=1;
    b1=b2;
    
    for ij=1:size(fire1,1)              %Update only the cells that burns
            
        i=fire1(ij,1);
        j=fire1(ij,2);
        
            %If the tree in these cell burns, than it will in some case turn on his neighbour.
            if (map(i,j,1)==1 || map(i,j,1)==2 || map(i,j,1)==3 || map(i,j,1)==4 || map(i,j,1)==5 || map(i,j,1)==6)
                
                switch map(i,j,1)                   %Different probability of propagation: counting intensity of the fire, from 1 to 6;
                    case 1
                        Intensity=0.7;
                    case 2
                        Intensity=0.85;
                    case {3,4}
                        Intensity=1;
                    case 5
                        Intensity=0.85;
                    case 6
                        Intensity=0.7;
                    case 8
                        Intensity=0.4;
                end
                
                % 3.3.1. Lighting caused by proximity
                for k=1:8   %Take a vector from the xy-Matrix;
                    
                    %Adjustement of the probability of propagation in specific direction: counting on the wind direction and the wind intensity.
                    WindCorrection=(WindDir(1)*xy1(k,1)+WindDir(2)*xy1(k,2))*WindInt *CorrFactorWind;                     %RealProp1 must always be between 0 and 1: CorrFactorWind garanty this condition;
                    SlopeCorrection=(Surface(i,j,2)*xy1(k,1)+Surface(i,j,3)*xy1(k,2))*Surface(i,j,4) *CorrFactorSlope;    %RealProp1 must always be between 0 and 1: CorrFactorSlope garanty this condition;
                    RealProp1=(Propag + WindCorrection + SlopeCorrection) *Intensity *Surface(i,j,5);         
                    
                    if ((i+xy1(k,1)>=1 && i+xy1(k,1)<=N) && (j+xy1(k,2)>=1 && j+xy1(k,2)<=N))
                        
                        %If a random number is smaller than the RealPropag1 and the neighbour is not burning or is not already burned, than set on fire the tree.
                        if (rand<RealProp1 && map(i+xy1(k,1),j+xy1(k,2),1)==0)
                            map(i+xy1(k,1),j+xy1(k,2),2)=1;         %Position: current location[i,j] + distance to the neighbour[xy1(k,1),xy1(k,2)];
                    
                            fire2(f,1)=i+xy1(k,1);                  %Store the burning cells in fire2;
                            fire2(f,2)=j+xy1(k,2);
                            f=f+1;
                        end
                        if (rand<RealProp1 && map(i+xy1(k,1),j+xy1(k,2),1)==8)
                            map(i+xy1(k,1),j+xy1(k,2),2)=1;
                    
                            fire2(f,1)=i+xy1(k,1);
                            fire2(f,2)=j+xy1(k,2);
                            f=f+1;
                            
                            time(tm)=t;                             %
                            tm=tm+1;
                            
                        end
                        
                    end
                end
                
                % 3.3.2. Additional lighting in case of diagonal wind
                if (WindDir(1)^2==1 && WindDir(2)^2==1)
                    %Adjust the propagation probability;
                    
                    RealProp2=(Propag*WindInt + (WindDir(1)*Surface(i,j,2)+WindDir(2)*Surface(i,j,3))*Surface(i,j,4) *CorrSlopeToWind) *Intensity *Surface(i,j,5);
                    for k=1:3                                   %Position: current location[i,j] + distance to the neighbour[xyWind2(k,1),xyWind2(k,2)];
                        
                        %If the neighbour is inside the map, proceed;
                        if ((i+xyWind2(k,1)*WindDir(1)>=1 && i+xyWind2(k,1)*WindDir(1)<=N) && (j+xyWind2(k,2)*WindDir(2)>=1 && j+xyWind2(k,2)*WindDir(2)<=N))
                         	
                            %If a random number is smaller than the RealProp2 and the neighbour is not burned yet, set onfire the tree.
                            if (rand<RealProp2 && map(i+xyWind2(k,1)*WindDir(1),j+xyWind2(k,2)*WindDir(2),1)==0)
                                map(i+xyWind2(k,1)*WindDir(1),j+xyWind2(k,2)*WindDir(2),2)=1;
                                
                                fire2(f,1)=i+xyWind2(k,1)*WindDir(1);             	%Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind2(k,2)*WindDir(2);
                                f=f+1;
                                
                            end
                        end
                    end
                
                % 3.3.3. Additional lighting in case of straight wind
                elseif (WindDir(1)^2==1 && WindDir(2)^2==0)
                    
                    %Help value
                    dir=0;
                    if (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==0)
                        dir=2;
                    elseif (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==1)
                        dir=1;
                    end
                    
                    RealProp3=(Propag*WindInt + (WindDir(1)*Surface(i,j,2)*dir)*Surface(i,j,4) *CorrSlopeToWind) *Intensity *Surface(i,j,5);
                    for k=1:2
                        if ((i+xyWind3(k,1)*WindDir(1)>=1 && i+xyWind3(k,1)*WindDir(1)<=N) && (j+xyWind3(k,2)>=1 && j+xyWind3(k,2)<=N))
                            if (rand<RealProp3 && map(i+xyWind3(k,1)*WindDir(1),j+xyWind3(k,2),1) == 0)
                                map(i+xyWind3(k,1)*WindDir(1),j+xyWind3(k,2),2)=1;
                                
                                fire2(f,1)=i+xyWind3(k,1)*WindDir(1);             	%Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind3(k,2);
                                f=f+1;
                                
                            end
                        end
                    end
                elseif (WindDir(1)^2==0 && WindDir(2)^2==1)
                    
                    %Help value
                    dir=0;
                    if (Surface(i,j,2)^2==0 && Surface(i,j,3)^2==1)
                        dir=2;
                    elseif (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==1)
                        dir=1;
                    end
                    
                    RealProp4=(Propag*WindInt + (WindDir(2)*Surface(i,j,3)*dir)*Surface(i,j,4) *CorrSlopeToWind) *Intensity *Surface(i,j,5);
                    for k=1:2
                        if ((i+xyWind3(k,2)>=1 && i+xyWind3(k,2)<=N) && (j+xyWind3(k,1)*WindDir(2)>=1 && j+xyWind3(k,1)*WindDir(2)<=N))
                            if (rand<RealProp4 && map(i+xyWind3(k,2),j+xyWind3(k,1)*WindDir(2),1) == 0)
                                map(i+xyWind3(k,2),j+xyWind3(k,1)*WindDir(2),2)=1;
                                
                                fire2(f,1)=i+xyWind3(k,2);                          %Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind3(k,1)*WindDir(2);
                                f=f+1;
                                
                            end
                        end
                    end
                end
                
                % 3.3.4 Additionat lighting in case of diagonal slope
                %We use for the slope the same method as for the wind:
                %   instead using WindInt, that is constant for all the terrain,
                %   we use the size of the gradient, Surface(:,:,4),that
                %   change for each cell.
                %   Instead using the WindDir(:,{1,2}) we use the gradient
                %   direction, Surface(:,:,{2,3}).
                if (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==1)
                    %Adjust the propagation probability;
                    RealProp5=(Propag*Surface(i,j,4) + (WindDir(1)*Surface(i,j,2)+WindDir(2)*Surface(i,j,3))*WindInt *CorrWindToSlope) *Intensity *Surface(i,j,5);
                    for k=1:3                                   %Position: current location[i,j] + distance to the neighbour[xyWind2(k,1),xyWind2(k,2)];
                        
                        %If the neighbour is inside the map, proceed;
                        if ((i+xyWind2(k,1)*Surface(i,j,2)>=1 && i+xyWind2(k,1)*Surface(i,j,2)<=N) && (j+xyWind2(k,2)*Surface(i,j,3)>=1 && j+xyWind2(k,2)*Surface(i,j,3)<=N))
                         	
                            %If a random number is smaller than the RealProp5 and the neighbour is not burned yet, set on fire the tree.
                            if (rand<RealProp5 && map(i+xyWind2(k,1)*Surface(i,j,2),j+xyWind2(k,2)*Surface(i,j,3),1)==0)
                                map(i+xyWind2(k,1)*Surface(i,j,2),j+xyWind2(k,2)*Surface(i,j,3),2)=1;
                                
                                fire2(f,1)=i+xyWind2(k,1)*Surface(i,j,2);           	%%Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind2(k,2)*Surface(i,j,3);
                                f=f+1;
                                
                            end
                        end
                    end
                
                % 3.3.5 Additional lighting in case of straight slope
                elseif (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==0)
                    
                    %Help value
                    dir=0;
                    if (WindDir(1)^2==1 && WindDir(2)^2==0)
                        dir=2;
                    elseif (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==1)
                        dir=1;
                    end
                    
                    RealProp6=(Propag*Surface(i,j,4) + (WindDir(1)*Surface(i,j,2)*dir)*WindInt *CorrWindToSlope) *Intensity *Surface(i,j,5);
                    for k=1:2
                        if ((i+xyWind3(k,1)*Surface(i,j,2)>=1 && i+xyWind3(k,1)*Surface(i,j,2)<=N) && (j+xyWind3(k,2)>=1 && j+xyWind3(k,2)<=N))
                            if (rand<RealProp6 && map(i+xyWind3(k,1)*Surface(i,j,2),j+xyWind3(k,2),1) == 0)
                                map(i+xyWind3(k,1)*Surface(i,j,2),j+xyWind3(k,2),2)=1;
                                
                                fire2(f,1)=i+xyWind3(k,1)*Surface(i,j,2);           	%Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind3(k,2);
                                f=f+1;
                                
                            end
                        end
                    end
                elseif (Surface(i,j,2)^2==0 && Surface(i,j,3)^2==1)
                    
                    %Help value
                    dir=0;
                    if (Surface(i,j,2)^2==0 && Surface(i,j,3)^2==1)
                        dir=2;
                    elseif (Surface(i,j,2)^2==1 && Surface(i,j,3)^2==1)
                        dir=1;
                    end
                    
                    RealProp7=(Propag*Surface(i,j,4) + (WindDir(2)*Surface(i,j,3)*dir)*WindInt *CorrWindToSlope) *Intensity *Surface(i,j,5);
                    for k=1:2
                        if ((i+xyWind3(k,2)>=1 && i+xyWind3(k,2)<=N) && (j+xyWind3(k,1)*Surface(i,j,3)>=1 && j+xyWind3(k,1)*Surface(i,j,3)<=N))
                            if (rand<RealProp7 && map(i+xyWind3(k,2),j+xyWind3(k,1)*Surface(i,j,3),1) == 0)
                                map(i+xyWind3(k,2),j+xyWind3(k,1)*Surface(i,j,3),2)=1;
                                
                                fire2(f,1)=i+xyWind3(k,2);                              %Store the burning cells in fire2;
                                fire2(f,2)=j+xyWind3(k,1)*Surface(i,j,3);
                                f=f+1;
                                
                            end
                        end
                    end
                end
                
                % 3.3.6. Burning
                if (rand>Duration && (map(i,j,1)==1 || map(i,j,1)==2 || map(i,j,1)==3 || map(i,j,1)==4 || map(i,j,1)==5 || map(i,j,1)==6))
                    map(i,j,2)=map(i,j,1)+1;        %The cells continue to burn until they reach the value 7.
                    if map(i,j,2)==7;               %Store the burned cells in fire2;
                        burned(b2,1)=i;
                        burned(b2,2)=j;
                        b2=b2+1;
                    end
                end
                if (map(i,j,1)==1 || map(i,j,1)==2 || map(i,j,1)==3 || map(i,j,1)==4 || map(i,j,1)==5 || map(i,j,1)==6)
                 	fire2(f,1)=i;                   %Store the cells that are still burning;
                   	fire2(f,2)=j;
                  	f=f+1;
                end
            end
    end
    
    
    
    % 3.3.7. Setting the matrix Color: it rappresent the RGB-color of each cell
    for ij=1:size(fire2,1)                  %Update only the burning cells
        i=fire2(ij,1);
        j=fire2(ij,2);
            switch map(i,j,2)
                case {1,2}                  %If 1 or 2: starting burning;
                    Color(i,j,1)=1;         %Color: red;
                    Color(i,j,2)=0;
                    Color(i,j,3)=0;
                case 3                    	%If 3: burning medium;
                    Color(i,j,1)=1;         %Color: orange;
                    Color(i,j,2)=0.3;
                    Color(i,j,3)=0;
                case {4,5}                	%If 4 or 5: burning hard;
                    Color(i,j,1)=1;         %Color: jellow;
                    Color(i,j,2)=1;
                    Color(i,j,3)=0;
                case 6                  	%If 6: almost burned;
                    Color(i,j,1)=0.8;       %Color: orange;
                    Color(i,j,2)=0.2;
                    Color(i,j,3)=0;
            end
    end
    for ij=b1:(b2-1)                        %Update only the new burned cells
        i=burned(ij,1);
        j=burned(ij,2);
            switch map(i,j,2)
                case 7                   	%If 7: burned;          %Still plot the 'HohenKurven
                                            %Color: Black
                    
                  	hohe=(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude);
                    hohe=floor(hohe*100*detail)/(100*detail);
                    if hohe==0.1
                        col=0;
                    elseif hohe==0.2
                        col=0;
                    elseif hohe==0.3
                        col=0;
                    elseif hohe==0.4
                        col=0;
                    elseif hohe==0.5
                        col=0;
                    elseif hohe==0.6
                        col=0;
                    elseif hohe==0.7
                        col=0;
                    elseif hohe==0.8
                        col=0;
                    elseif hohe==0.9
                        col=0;
                    else 
                        col=1;
                    end
                    
                    Color(i,j,1)=0.3*(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude)*col + 0.6*(1-col);
                    Color(i,j,2)=0.3*(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude)*col + 0.3*(1-col);
                    Color(i,j,3)=0.3*(Surface(i,j,1)-MinAltitude)/(MaxAltitude-MinAltitude)*col;
            end
    end
    
    % 3.3.8. Number of burning tree
    burning=f-1;                                    %burned already calculated;
    notBurned=N*N-burning-(b2-1);
    burnedd=b2-1;
    
    % 3.3.9. PLOT AND MOVIE
    hold on
    
    clf                                             %Clear;
    
    if ID == 1;
    
    axis([0 N 0 N]);
    surf(Surface(:,:,1),Color,'EdgeColor','none')   %Plot the Surface with the color specified in the matrix Color;
    view([-90,90])
    grid off
    
    elseif ID == 2
        
    surf(Surface(:,:,1),Color,'EdgeColor','none')   %Plot the Surface with the color specified in the matrix Color;
    view([-70,25])
    grid off
    
    elseif ID == 3;
    axis([0 St 0 N*N]);
    plot(t,notBurned,'og',t,burning,'or',t,burnedd,'ok');
    legend('not burned','burning','burned','LOCATION','best');

    end
    
    drawnow                                  

% Video Frame
%F = getframe(fig);
%aviobj = addframe(aviobj,F);
    
    
    % 3.3.10. If the 80% of the forest is burned, than stop the simulation 
    if burned>0.8*N^2                                  %If more than 80% break: stop the for-loop;
        break;
    end 
end      
