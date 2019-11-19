function [par,par_best]=ErPSOupdate(par,par_best)
par.v1=par.v1+par.vv1;
par.v2=par.v2+par.vv2;
par.v3=par.v3+par.vv3;
par.v4=par.v4+par.vv4;
par.v5=par.v5+par.vv5;
par.v6=par.v6+par.vv6;

par.p1=par.p1+par.vp1;
par.p2=par.p2+par.vp2;
par.p3=par.p3+par.vp3;
par.p4=par.p4+par.vp4;
par.p5=par.p5+par.vp5;
par.p6=par.p6+par.vp6;
par.fit=ErPSOcompute_fitness(par);
%par.gain=PSOcompute_gain(par);
w=0.8;    %
c1=1;    %
c2=1;    %
par.vv1=w*par.vv1+c2*rand()*(par_best.v1-par.v1)+c1*rand()*(par.bestv1-par.v1);
par.vv2=w*par.vv2+c2*rand()*(par_best.v2-par.v2)+c1*rand()*(par.bestv2-par.v2);
par.vv3=w*par.vv3+c2*rand()*(par_best.v3-par.v3)+c1*rand()*(par.bestv3-par.v3);
par.vv4=w*par.vv4+c2*rand()*(par_best.v4-par.v4)+c1*rand()*(par.bestv4-par.v4);
par.vv5=w*par.vv5+c2*rand()*(par_best.v5-par.v5)+c1*rand()*(par.bestv5-par.v5);
par.vv6=w*par.vv6+c2*rand()*(par_best.v6-par.v6)+c1*rand()*(par.bestv6-par.v6);

par.vp1=w*par.vp1+c2*rand()*(par_best.p1-par.p1)+c1*rand()*(par.bestp1-par.p1);
par.vp2=w*par.vp2+c2*rand()*(par_best.p2-par.p2)+c1*rand()*(par.bestp2-par.p2);
par.vp3=w*par.vp3+c2*rand()*(par_best.p3-par.p3)+c1*rand()*(par.bestp3-par.p3);
par.vp4=w*par.vp4+c2*rand()*(par_best.p4-par.p4)+c1*rand()*(par.bestp4-par.p4);
par.vp5=w*par.vp5+c2*rand()*(par_best.p5-par.p5)+c1*rand()*(par.bestp5-par.p5);
par.vp6=w*par.vp6+c2*rand()*(par_best.p6-par.p6)+c1*rand()*(par.bestp6-par.p6);
if (par.fit>par.bestfit)                           
        par.bestfit=par.fit;                       
        par.bestv1=par.v1;                         
        par.bestv2=par.v2;
        par.bestv3=par.v3;
        par.bestv4=par.v4;
        par.bestv5=par.v5;
        par.bestv6=par.v6;
        
        par.bestp1=par.p1;
        par.bestp2=par.p2;
        par.bestp3=par.p3;
        par.bestp4=par.p4;
        par.bestp5=par.p5;
        par.bestp6=par.p6;
        if (par.bestfit>par_best.bestfit)             
            par_best.bestfit=par.bestfit;             
            %par_best.bestgain=par.bestgain;
            par_best.v1=par.bestv1;                    
            par_best.v2=par.bestv2;
            par_best.v3=par.bestv3;
            par_best.v4=par.bestv4;
            par_best.v5=par.bestv5;
            par_best.v6=par.bestv6;
            
            par_best.p1=par.bestp1;
            par_best.p2=par.bestp2;
            par_best.p3=par.bestp3;
            par_best.p4=par.bestp4;
            par_best.p5=par.bestp5;
            par_best.p6=par.bestp6;
        end
end

