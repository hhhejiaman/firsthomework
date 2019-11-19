clc;
N=10;                 
T=50;                

par=struct([]);
for i=1:N
    par(i).v1=223.881+6.888*rand();           
    par(i).v2=223.881+6.888*rand();          
    par(i).v3=209.79+4.496*rand();           
    par(i).v4=200+4.082*rand(); 
    par(i).v5=209.79+4.496*rand(); 
    par(i).v6=202.703+4.194*rand(); 

    par(i).p1=0.8+0.7*rand();               
    par(i).p2=0.5+0.7*rand();                
    par(i).p3=0.05+0.35*rand();                
    par(i).p4=0.01+0.14*rand();
    par(i).p5=0.05+0.35*rand();
    par(i).p6=0+0.3*rand();
  
    par(i).vv1=-0.134+0.268*rand();        
    par(i).vv2=-0.134+0.268*rand();        
    par(i).vv3=-0.134+0.268*rand();        
    par(i).vv4=-0.134+0.268*rand();
    par(i).vv5=-0.134+0.268*rand();
    par(i).vv6=-0.134+0.268*rand();
    
    
    par(i).vp1=-0.001+0.002*rand();        
    par(i).vp2=-0.001+0.002*rand();        
    par(i).vp3=-0.001+0.002*rand();       
    par(i).vp4=-0.001+0.002*rand();
    par(i).vp5=-0.001+0.002*rand();
    par(i).vp6=-0.001+0.002*rand();
   
    par(i).bestv1=par(i).v1;               
    par(i).bestv2=par(i).v2;
    par(i).bestv3=par(i).v3;
    par(i).bestv4=par(i).v4;
    par(i).bestv5=par(i).v5;
    par(i).bestv6=par(i).v6;
    
    par(i).bestp1=par(i).p1;
    par(i).bestp2=par(i).p2;
    par(i).bestp3=par(i).p3;
    par(i).bestp4=par(i).p4;
    par(i).bestp5=par(i).p5;
    par(i).bestp6=par(i).p6;
    
    par(i).fit=0;
    par(i).bestfit=0;
end

par_best=par(1);

 for k=1:T
    for p=1:N
        [par(p),par_best]=ErPSOupdate(par(p),par_best);   
      
        par(p).fit
    end
    
    
end
ErPSOcompute_fitness(par_best)
ErPSOAverageGain(par_best)
ErPSOGainFlatness(par_best)
disp('finish');
    
    
    
