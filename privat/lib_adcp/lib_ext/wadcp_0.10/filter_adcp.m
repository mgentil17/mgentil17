function filter_adcp

global nbens theta WN WS WB IMM ens_int f a_t ouv Kc EC0 B SL0 R1 Lt binpos yrange
global time ve vn vu dir spd a1 a2 a3 a4 amean pitch roll heading c1 c2 c3 c4 cmean BI Mgl
global Batt Tb_adcp Depth_adcp Sb_adcp mode_exp beam_choice cor_threshold fchoice cellskip fidlog 


%% FILTRE DE SURFACE -- recherche de la dernière cellule avant la surface
%
% * plusieurs possibilités : 
%   (1) à partir du capteur de pression
%   (2) à partir du max de signal rétrodiffusé
%   (3) à partir d'un travail sur la corrélation des faisceaux



if mode_exp==0
bv=char('Pressure signal','Backscatter','Correlation');
[fchoice,vresult]=listdlg('PromptString','Processing option - filtering the sea surface','ListString',bv);
nbv=length(fchoice);
while nbv>1
    diplay('Carefull, choose only one filter'); 
    bv=char('Pressure signal','Backscatter','Correlation');
    [fchoice,vresult]=listdlg('PromptString','Processing option - filtering the sea surface','ListString',bv);
    nbv=length(fchoice);
end


prompt = {'How many cells below the detected sea surface to skip?'};
dlg_title = 'Processing options';
num_lines = 1;
def = {'0'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
cellskip=eval(char(tmp{1}));

end

    pos=zeros(nbens,1);

if fchoice==1

    
    
%     hh = waitbar(0,'Filtration [step 1] in progress...');
%     for i=1:nbens
%         waitbar(i/nbens,hh,strcat('Filtration [step 1] in progress...--',int2str(i/nbens*100),'%'));
%         tmp=find(binpos<=Depth_adcp(i));
%         if isempty(tmp)
%            pos(i)=NaN;
%            continue
%         else
%         pos(i)=tmp(length(tmp))-cellskip;
%             if pos(i)<=0
%             pos(i)=NaN;
%             end
%         end
%     end
% close(hh)

binpos2=repmat(binpos,nbens,1);
Depth2=repmat(Depth_adcp,1,length(binpos));


pos=ones(nbens,length(binpos));
pos(binpos2>Depth2)=0;

elseif fchoice==2
    
    if mode_exp==0
beam_txt=char('beam 1','beam 2','beam 3','beam 4 [RDI only]','averaged beams');
[beam_choice,presult]=listdlg('PromptString','Beam choice for sea surface detection from backscatter signal','ListString',beam_txt);
beam=beam_txt(beam_choice,:);
    end
if beam_choice==1
    ampl=a1;
elseif beam_choice==2
    ampl=a2;    
elseif beam_choice==3
    ampl=a3;
elseif beam_choice==4
    ampl=a4;    
elseif beam_choice==5
    ampl=amean;    
end
        for i=1:nbens
           tmp=find(ampl(i,:)==max(ampl(i,:)));
           pos(i)=tmp(1)-cellskip;
           if pos(i)<=0
               pos(i)=NaN;
           end
        end    

elseif fchoice==3
if mode_exp==0
beam_txt=char('beam 1','beam 2','beam 3','beam 4 [RDI only]','averaged beams');
[beam_choice,presult]=listdlg('PromptString','Beam choice for sea surface detection from backscatter signal','ListString',beam_txt);
beam=beam_txt(beam_choice,:);
end

if beam_choice==1
    cor=a1;
elseif beam_choice==2
    cor=a2;    
elseif beam_choice==3
    cor=a3;
elseif beam_choice==4
    cor=a4;    
elseif beam_choice==5
    cor=amean;    
end
    
    
    pos=zeros(nbens,1);
tmp_cor=zeros(nbens,WN-1);
if mode_exp==0
prompt = {'Differences in correlation signal to account for surface detection'};
dlg_title = 'Processing options';
num_lines = 1;
def = {'5'};
tmp = inputdlg(prompt,dlg_title,num_lines,def);
cor_threshold=eval(char(tmp{1}));
end

for i=1:nbens
   
   tmp_cor(i,:)=cor(i,2:WN)-cor(i,1:WN-1);
   

   
   for j=2:WN-2
       
       if (tmp_cor(i,j)<-cor_threshold && tmp_cor(i,j-1)>cor_threshold) || (tmp_cor(i,j-1)>cor_threshold && tmp_cor(i,j+1)>cor_threshold)
            pos(i)=j; % ou +1...à voir
            break
       end
       
       if j==WN-2
            pos(i)=WN-1;
       end
   end  
end
pos2=pos;
pos(pos==WN-1)=NaN;
delta_pos=pos(2:nbens)-pos(1:nbens-1);
% 
% for i=1:nbens-1
%     if abs(delta_pos(i))>=2
%     %k=1;
%     %while abs(pos(i+1+k)-pos(i))>2
%     %k=k+1;
%         if abs(pos(i+2)-pos(i))<=2
%         pos(i+1)=floor(0.5*(pos(i)+pos(i+2)));
%         elseif abs(pos(i+3)-pos(i))<=2
%         pos(i+1)=floor(0.5*(pos(i)+pos(i+3)));
%         else
%         pos(i+1)=floor(0.5*(pos(i)+pos(i+4)));
%         end
%     %pos(i+1)=floor(0.5*(pos(i)+pos(i+k+1)));
%     delta_pos=pos(2:nbens)-pos(1:nbens-1);
%     
%     end
% end
% 
% 
% nanpos=isnan(pos);
% 
% for i=2:nbens-1
%    if isnan(pos(i))
%        pos(i)=(pos(i+1)+pos(i-1))/2;
%    end
% 
%    %if isnan(pos(i)) & pos(i+1)>1
%    %pos(i+1)=NaN;
%    %end
%    
%    %if isnan(pos(i)) & pos(i-1)>1
%    %pos(i-1)=NaN;
%    %end
%        
% end

pos=pos-cellskip;
pos(pos<0)=0;

end

%% application du filtre choisi
%





% hh = waitbar(0,'Filtration [step 2] in progress...');
% for i=1:nbens
%             waitbar(i/nbens,hh,strcat('Filtration [step 2] in progress...--',int2str(i/nbens*100),'%'));
%      if isnan(pos(i))
% a1(i,:)=NaN;
% a2(i,:)=NaN;
% a3(i,:)=NaN;
% a4(i,:)=NaN;
% amean(i,:)=NaN;
% 
% ve(i,:)=NaN;
% vn(i,:)=NaN;
% vu(i,:)=NaN;
% 
% c1(i,:)=NaN;
% c2(i,:)=NaN;
% c3(i,:)=NaN;
% c4(i,:)=NaN;
% cmean(i,:)=NaN;
% 
% spd(i,:)=NaN;
% dir(i,:)=NaN;
% 
% %BI(:,i)=NaN;
% %Mgl(i,:)=NaN;
% 
%      else
% a1(i,binpos>binpos(pos(i)))=NaN;
% a2(i,binpos>binpos(pos(i)))=NaN;
% a3(i,binpos>binpos(pos(i)))=NaN;
% a4(i,binpos>binpos(pos(i)))=NaN;
% amean(i,binpos>binpos(pos(i)))=NaN;
% 
% ve(i,binpos>binpos(pos(i)))=NaN;
% vn(i,binpos>binpos(pos(i)))=NaN;
% vu(i,binpos>binpos(pos(i)))=NaN;
% 
% c1(i,binpos>binpos(pos(i)))=NaN;
% c2(i,binpos>binpos(pos(i)))=NaN;
% c3(i,binpos>binpos(pos(i)))=NaN;
% c4(i,binpos>binpos(pos(i)))=NaN;
% cmean(i,binpos>binpos(pos(i)))=NaN;
% 
% spd(i,binpos>binpos(pos(i)))=NaN;
% dir(i,binpos>binpos(pos(i)))=NaN;
% 
% % BI(yrange>binpos(pos(i)),i)=NaN;
% %Mgl(i,yrange>binpos(pos(i)))=NaN;
% 
%      end
%  end
% 
% % % figure
% % % plot(Depth_adcp)
% % % 
% % % for i=1:nbens
% % %             waitbar(i/nbens,hh,strcat('Filtration [step 2] in progress...--',int2str(i/nbens*100),'%'));
% % %      if Depth_adcp(i)<0.01
% % % a1(i,:)=NaN;
% % % a2(i,:)=NaN;
% % % a3(i,:)=NaN;
% % % a4(i,:)=NaN;
% % % amean(i,:)=NaN;
% % % 
% % % ve(i,:)=NaN;
% % % vn(i,:)=NaN;
% % % vu(i,:)=NaN;
% % % 
% % % c1(i,:)=NaN;
% % % c2(i,:)=NaN;
% % % c3(i,:)=NaN;
% % % c4(i,:)=NaN;
% % % cmean(i,:)=NaN;
% % % 
% % % spd(i,:)=NaN;
% % % dir(i,:)=NaN;
% % %      end
% % % %BI(:,i)=NaN;
% % % %Mgl(i,:)=NaN;
% % %      end
% % 
% close(hh)

a1(pos==0)=NaN;
a2(pos==0)=NaN;
a3(pos==0)=NaN;
a4(pos==0)=NaN;
amean(pos==0)=NaN;

ve(pos==0)=NaN;
vn(pos==0)=NaN;
vu(pos==0)=NaN;

c1(pos==0)=NaN;
c2(pos==0)=NaN;
c3(pos==0)=NaN;
c4(pos==0)=NaN;
cmean(pos==0)=NaN;

spd(pos==0)=NaN;
dir(pos==0)=NaN;
%     end
%BI(:,i)=NaN;
%Mgl(i,:)=NaN;


% a1(c1<100)=NaN;
% a2(c1<100)=NaN;
% a3(c1<100)=NaN;
% a4(c1<100)=NaN;
% amean(c1<100)=NaN;
% 
% ve(c1<100)=NaN;
% vn(c1<100)=NaN;
% vu(c1<100)=NaN;
% 
% 
% c2(c1<100)=NaN;
% c3(c1<100)=NaN;
% c4(c1<100)=NaN;
% cmean(c1<100)=NaN;
% 
% spd(c1<100)=NaN;
% dir(c1<100)=NaN;
% c1(c1<100)=NaN;
%BI(:,i)=NaN;
%Mgl(i,:)=NaN;

    fprintf(fidlog,'Postprocessing information\n');
    fprintf(fidlog,'*****************************************************\n');    
    fprintf(fidlog,'Filtering option (1 pressure, 2 correlation, 3 Backscatter) : fchoice : %1.0f\n', fchoice);    
    fprintf(fidlog,'Filtering option Beam choice : %1.0f\n', beam_choice);
    fprintf(fidlog,'Skipped cells : cellskip : %2.0f\n', cellskip);
    if fchoice==2
    fprintf(fidlog,'Correlation threshold %2.2f\n', cor_threshold);    
    end
    fprintf(fidlog,'\n');


