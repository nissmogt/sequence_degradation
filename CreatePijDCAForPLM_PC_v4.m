function CreatePijDCAForPLM_PC_v4(pseudocount,IsingMatrix,SequenceMatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    %all relavent parameters of sequence file
    load(SequenceMatrix);
    load(IsingMatrix);
    %load('../mc_Pi.mat');
    
    tic
    
    [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount,N,q);
    
    Pij_plm=zeros(N,N,q,q);
    
    StoreDI=zeros(N,N); %for DI and corrected DI.... i, j
    counter=0;
    for i=1:(N-1)
        for j=(i+1):N
            counter = counter+1;
            
            [PijPLM] = ComputeConditionalPij(pseudocount,N,i,j,q,B,weights,h,J,Pij_true,Pi_true);
            
            if 1-sum(sum(PijPLM))>0.00001
                error('Error for i=%d, j=%d',i,j);
            end
            
            Pij_plm(i,j,1:q,1:q)=PijPLM(1:q,1:q);
            
            
            
            DI=ComputeDI(PijPLM,Pi_true,i,j);
            
            StoreDI(i,j)=DI;
            StoreDI(j,i)=StoreDI(i,j);
            
            
            
            %%fprintf(fOut,'%d,%d,%f\n',i,j,DI);
        end
    end
    

        
        
    
    %calculate corrected
    APCCorrectedDI=zeros(N,N);
    DI_means=mean(StoreDI)*N/(N-1);
    DI_means_all=mean(mean(StoreDI))*N/(N-1);
    Corr_DI=StoreDI-DI_means'*DI_means/DI_means_all;
    output=[];
    for i=1:(N-1)
        for j=(i+1):N
            output=[output;[i,j,StoreDI(i,j),Corr_DI(i,j)]];
            APCCorrectedDI(i,j)=Corr_DI(i,j);
            APCCorrectedDI(j,i)=Corr_DI(i,j);
        end
    end
    FileOut=sprintf('DI_PC-%.2f.txt',pseudocount);
    
    dlmwrite(FileOut,output,'precision',5)
    
    %fclose(fOut);
    MatOut=sprintf('PijPLM_PC-%.2f.mat',pseudocount);
    save(MatOut,'Pij_plm');


    toc
    
    %print mean rawDI
    meanRawDI=mean(StoreDI);
    FileOut=sprintf('mean_rawDI_PC-%.2f.txt',pseudocount);
    fOut=fopen(FileOut,'w');
    for i=1:N
        fprintf(fOut,'%d,%f\n',i,meanRawDI(i));
    end
    fclose(fOut);


end

function [DI] = ComputeDI(PijPLM,Pi,i,j)
    tiny=10^-100;
    PiPj=Pi(i,:)'*Pi(j,:);
    DI = trace( PijPLM' * log( (PijPLM+tiny)./(PiPj+tiny) ) );
    %DI=sum(sum(PijPLM.*log(PijPLM./PiPj)));


end


function [PijPLM] = ComputeConditionalPij(pseudocount,N,i,j,q,B,weights,h,J,Pij,Pi) %for a particular i and j
    PijPLM=zeros(q,q);
    Counter=zeros(q,q);
    z=0;
    
    for LetterI=1:q
        for LetterJ=1:q            
            
            %Compute numerator
            %func = h(i,LetterI) + h(j,LetterJ) + J(i,j,LetterI,LetterJ);
            func = J(i,j,LetterI,LetterJ);
            
            temp=func; %denominator
            
            z=z+Pij(i,j,LetterI,LetterJ)*exp(temp); %denominator
            
            %Set distribution
            PijPLM(LetterI,LetterJ)=PijPLM(LetterI,LetterJ) + Pij(i,j,LetterI,LetterJ)*exp(func);
            
        end
    end
    


    PijPLM(1:q,1:q)=PijPLM(1:q,1:q)./z;  %/sum(sum(PijPLM));
    
    
    PijPLM=(1.-pseudocount)*PijPLM + pseudocount/(q*q);
    
    [mu1,mu2] = compute_mu(i,j,PijPLM,Pi,q);
    PijPLM=PijPLM.*(mu1'*mu2);
    PijPLM=PijPLM/sum(sum(PijPLM));
    
%     [mu1,mu2] = compute_mu(i,j,W,Pi,q);
%     W.*(mu1'*mu2);
    

end

function [Pij,Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight,N,q)
% adds pseudocount

    Pij = (1.-pseudocount_weight)*Pij_true + pseudocount_weight/q/q*ones(N,N,q,q);
    Pi = (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q*ones(N,q);

    scra = eye(q);

    for i=1:N
        for alpha = 1:q
            for beta = 1:q
               Pij(i,i,alpha,beta) =  (1.-pseudocount_weight)*Pij_true(i,i,alpha,beta) + pseudocount_weight/q*scra(alpha,beta);
            end
        end
    end 
end


function [mu1,mu2] = compute_mu(i,j,W,P1,q)

    epsilon=1e-4;
    diff =1.0;
    mu1 = ones(1,q)/q;
    mu2 = ones(1,q)/q;
    pi = P1(i,:);
    pj = P1(j,:);

    while ( diff > epsilon )

        scra1 = mu2 * W';
        scra2 = mu1 * W;

        new1 = pi./scra1;
        new1 = new1/sum(new1);

        new2 = pj./scra2;
        new2 = new2/sum(new2);

        diff = max( max( abs( new1-mu1 ), abs( new2-mu2 ) ) );

        mu1 = new1;
        mu2 = new2;

    end
end
