function [Xistruct] = multiD_Lambda(Thetalib,lambdavals)
% Run a multidimensional Lambda sweep
% Thetalib is a structure containing variables for data matrix
% lambdavals is a structure containing parameters for defining the sparsity
% threshold vector.
% plottag determins

Theta = Thetalib.Theta;
normTheta = Thetalib.normTheta;
dx = Thetalib.dx;
polyorder = Thetalib.polyorder;
usesine = Thetalib.usesine;

numlambda = lambdavals.numlambda;
lambdastart = lambdavals.lambdastart;
lambdaend = lambdavals.lambdaend;
if isfield(lambdavals,'method')
    sparsify_tag = lambdavals.method;
else
    sparsify_tag = 'threshold';
end
    

% make a vector of threshold values
Lambda = logspace(lambdastart,lambdaend, numlambda);

% get sizes and lengths of data and library
[ntimeseries, nstates] = size(dx);
[ntimeseries, nfunc] = size(Theta);

% declare stuff to be correct size
XiLambda_L2 = zeros(nfunc,numlambda,nstates);


% Run sparse-regression independently on each state variable and do a Least-squares
% regression on the output for each lambda.
for nn = 1:nstates % for each state variable
    for jj = 1:numlambda % and each threshold value
        if strcmp('elastic', sparsify_tag)
            XiLambda_L2(:,jj,nn) = sparsifyDynamics_elastic(Theta,dx(:,nn), Lambda(jj), 1);
        else
            XiLambda_L2(:,jj,nn) = sparsifyDynamics(Theta,dx(:,nn), Lambda(jj), 1);
        end
        
    end
end

% create a matrix of threshold values, so that all permutations of the
% recovered equations for each state are explored.
lambda_index = permn(1:numlambda,nstates);

% initialize the loop values
Xiold= 0;
nn = 1;
Xilib= [];

% calcualate validated error
for ll = 1:length(lambda_index) % for each sparse thershold, or model 
                                % defined by thepermuation matrix

    lambdanum_state = lambda_index(ll,:);
    % build the Xi matrix for this test (columns are coefficents for each
    % state variable, each row is a coefficient for a function in the library)
    Xibuild = [];
    for ii = 1:length(lambdanum_state) 
        if length(normTheta)>1 % check if the function library was normalized
            Xibuild = [Xibuild XiLambda_L2(:,lambdanum_state(ii),ii)./normTheta'];
        else 
            Xibuild = [Xibuild XiLambda_L2(:,lambdanum_state(ii),ii)];
        end
        Xibuild(abs(Xibuild)<1e-10) =0;
    end
    % put all coefficients into a vector to check if we have tried this Xi
    % before
    Xitest = reshape(Xibuild, nstates*nfunc,1);
    
    [junk, nlib] = size(Xilib); % get the size of the library
    
    % only run the rest of this for loop for unique Xi
    if (ll>1)&& any(all(abs((repmat(Xitest,1,nlib)-Xilib)<1e-10),1)) 
        continue
    else % if the coefficients are novel, add them to the library
        Xilib = [Xilib Xitest]; 
       
    end

    % find the number of non-zero coefficients in the current library
    % guess
    numcoeff(nn) = nnz(Xibuild);
    Xicomb{nn} = Xibuild;
     nn = nn+1;    
    
end
% store output into a structure

Xistruct.numcoeff = numcoeff;
Xistruct.Xicomb = Xicomb;
Xistruct.lambda = Lambda;
end


