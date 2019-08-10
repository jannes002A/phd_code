%Skipt Homotopy Paper

load start_butan2.txt -ascii

M=1000;
%par_hom=0:0.1:0.4;
par_hom=0.8;
anzSamples=1000;
avg_hit=zeros(1,length(par_hom));
avg_t=zeros(1,length(par_hom));
avg_acc=zeros(1,length(par_hom));

for r=1:length(par_hom)
    r
    hit_total=0;
    t_total=0;
    acc_total=0;
    for i=1:M
        [hit,t,acc] = Homotopy_butan(start_butan2,par_hom(r),anzSamples); 
        hit_total=hit_total+hit;
        t_total=t_total+t;
        acc_total=acc_total+acc;
    end
    avg_hit(r) = hit_total/M;
    avg_t(r) =   t_total/M;
    avg_acc(r) = acc_total/M;
end
%%

load start_butan2.txt -ascii

M=1000;
par_hom=0:0.1:0.8;
%par_hom=0.5;
anzSamples=1000;
avg_hit1=zeros(1,length(par_hom));
avg_t1=zeros(1,length(par_hom));
avg_acc1=zeros(1,length(par_hom));

for r=1:length(par_hom)
    r
    hit_total=0;
    t_total=0;
    acc_total=0;
    for i=1:M
        [hit,t,acc] = Homotopy_butan(start_butan2,par_hom(r),anzSamples); 
        hit_total=hit_total+hit;
        t_total=t_total+t;
        acc_total=acc_total+acc;
    end
    avg_hit1(r) = hit_total/M;
    avg_t1(r) = t_total/M;
    avg_acc1(r) = acc_total/M;
end
% %%
% load start_butan2.txt -ascii

% M=100;
% par_hom=0:0.02:0.16;
% %par_hom=0.5;
% anzSamples=100;
% avg_hit2=zeros(1,length(par_hom));
% avg_t2=zeros(1,length(par_hom));
% avg_acc2=zeros(1,length(par_hom));
% 
% for r=1:length(par_hom)
%     r
%     hit_total=0;
%     t_total=0;
%     acc_total=0;
%     for i=1:M
%         [hit,t,acc] = Homotopy_butan(start_butan2,par_hom(r),anzSamples); 
%         hit_total=hit_total+hit;
%         t_total=t_total+t;
%         acc_total=acc_total+acc;
%     end
%     avg_hit2(r) = hit_total/M;
%     avg_t2(r) = t_total/M;
%     avg_acc2(r) = acc_total/M;
% end
