
clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;

%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;

%Number of Nodes in the field
n=50;

%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules) 
%Initial Energy 
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=0.1;
%\alpha
a=1;

%maximum number of rounds
rmax=5

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%


%Range of temparature from a=0 to b=50

a=0;
b=50;

%Range of temparature

%Computation of do
doo=sqrt(Efs/Emp); 

%Creation of the random Sensor Network
figure(1);

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Randomly creating sensor data from 0 to 50
	S(i).temp= (b-a).*rand(1,1) + a;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%End of Randomly creating sensor data from 0 to 50
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
        %%%%plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a)
        S(i).ENERGY=1;
        %%%%plot(S(i).xd,S(i).yd,'+');
         hold on;
    end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
    
        
%First Iteration
figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
for r=0:1:rmax
    r
    packets_TO_BS_per_round=0;
  %Election Probability for Normal Nodes
  pnrm=( p/ (1+a*m) );
  %Election Probability for Advanced Nodes
  padv= ( p*(1+a)/(1+a*m) );
    
  %Operation for heterogeneous epoch
  if(mod(r, round(1/pnrm) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

 %Operations for sub-epochs
 if(mod(r, round(1/padv) )==0)
    for i=1:1:n
        if(S(i).ENERGY==1)
            S(i).G=0;
            S(i).cl=0;
        end
    end
  end

 
hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round

PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');


STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %Election of Cluster Heads for normal nodes
 if( ( S(i).ENERGY==0 && ( temp_rand <= ( pnrm / ( 1 - pnrm * mod(r,round(1/pnrm)) )) ) )  )

            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
			
			msg=['Normal node cluster head  ',num2str(C(cluster).id),'  with X  ',num2str(X(cluster)),'  with Y  ',num2str(Y(cluster))];
			disp(msg);
			msg1=['temp of the AN in Normal group ',num2str(S(C(cluster).id).temp)];
		    disp(msg1);
			
			
            cluster=cluster+1;
			
            
            %Calculation of Energy dissipated
            distance;
            if (distance>doo)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=doo)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     


 %Election of Cluster Heads for Advanced nodes
 if( ( S(i).ENERGY==1 && ( temp_rand <= ( padv / ( 1 - padv * mod(r,round(1/padv)) )) ) )  )
        
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';
            S(i).G=100;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
			X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
			
			
			msg=['Advanced node cluster head  ',num2str(C(cluster).id),'  with X  ',num2str(X(cluster)),'with Y  ',num2str(Y(cluster))];
			disp(msg);
			msg1=['temp of the AN in Advanced group ',num2str(S(C(cluster).id).temp)];
		    disp(msg1);
			
			
			
            cluster=cluster+1;
			
			

            
            %Calculation of Energy dissipated
            distance;
            if (distance>doo)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=doo)
               S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
            end
        end     
    
    end
  end 
end



STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;

%Election of Associated Cluster Head for Normal Nodes
count_nodes=0;
for i=1:1:n

   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
	   min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            if (min_dis>doo)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=doo)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %Energy dissipated
       if(min_dis>0)
           S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
           PACKETS_TO_CH(r+1)=n-dead-cluster+1; 
       end
        
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
	   
	   msg2=['for NOde-- ',num2str(i),'min dis calculated was  ',num2str(S(i).min_dis),'and Min distance cluster was ',num2str(S(i).min_dis_cluster), '----its X and Y Values were  ',num2str(S(i).xd),'  and ',num2str(S(i).yd),'  temp= ',num2str(S(i).temp)];
	   disp(msg2);
      % xlswrite('C:\Users\Admin\Downloads\SensorData.xlsx',num2str(S(i).xd),'Sheet1','A2');
        %xlswrite('C:\Users\Admin\Downloads\SensorData.xlsx',num2str(S(i).yd),'Sheet1','B2');
	   msg3=['THis is the Advanced node  ',num2str(C(min_dis_cluster).id)];
	   disp(msg3);
	   
	   %data_from_node(count_nodes++)=[i,S(i).min_dis_cluster,S(i).xd,S(i).yd,S(i).temp,C(min_dis_cluster).id];
	   
           
   end
 end
end
count_nodes
for i=1:count_nodes
data_from_node(i)
end
hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;
end


%Tn=T(:,:,3);
%xlswrite(Tn,'G:\CD\Book1.xlsx','Sheet',1,'Range','A1');
%writetable(T,'C:\Users\HP\Downloads\New folder\data.xlsx','Sheet',1,'Range','A1'),




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


