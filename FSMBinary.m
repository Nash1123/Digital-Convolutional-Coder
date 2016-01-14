classdef FSMBinary < handle
% FSM: Finite state machine encoder and decoder
    properties
        % State transition and output maps:
        % x[t+1] = stateMap(x[t],b[t]+1)
        % c[t] = outMap(x[t],b[t]+1)
        stateMap;
        outMap;
        % Tail bit parameters
        ntail; % number of tail bits to be appended to the input 
        bitsIn_TailBitsAppended; % bitsIn appended with tail bits
        bitsOut; %bitOuts
        x; % states sequence
        gamma;
        alpha;
        beta;
        llrApp;
    end

    methods
        % Constructor. Empty for now
        function fsm = FSMBinary()
        end
        
        % Sets FSM for a recursive convolution code [1, (1+D+D^2)/(1+D^2)]
        function setMapConv(fsm)
            nx = 4; % number of states
            nin = 2; % number of inputs
            fsm.stateMap = zeros(nx,nin);
            fsm.outMap = zeros(nx,nin);
            fsm.ntail = 2;
            for iin = 1:nin
                ib = iin-1;
                for ix=1:nx
                    % Get state values from the two bits of the current state
                    x1 = bitget(ix-1,2); % z[t-1]
                    x2 = bitget(ix-1,1); % z[t-2]                    
                    % Set the state map to the index of the x[t+1] given
                    % x[t] = ix-1 and b[t] = iin-1. This should be an expression
                    % in terms of x1, x2 and b.
                    fsm.stateMap(ix,iin) = 2*mod(x1+x2+ib,2)+x1;
                    % Set the values of the output given x[t]=ix and
                    % b[t]=iin.
                    fsm.outMap(ix,iin) = 2*ib+mod(x1+ib,2);  
                end
            end
        end
        
        % Encode: Adds tail bits and generate a sequence of output bits.
        function [bitsOut,x] = encode(fsm,bitsIn)            
            n = size(bitsIn,1);
            b = zeros(n+2,1);
            for t = 1:n;
                b(t) = bitsIn(t); % Generate a vector with the tail bits appended
            end
            nt = size(b,1); % Get the vector length
            fsm.x = zeros(nt+1,1);
            fsm.bitsOut = zeros(2*nt,1);
            fsm.bitsIn_TailBitsAppended = b;
            cout = zeros(nt,1);
            for t = 1:nt
                cout(t) = fsm.outMap(fsm.x(t)+1,b(t)+1);  % Use fsm.outMap to get the output bits
                fsm.x(t+1) = fsm.stateMap(fsm.x(t)+1,b(t)+1); % Use fsm.stateMap to get the states sequence    
            end
            cout = dec2bin(cout)-'0';
            for t = 1:nt;
                fsm.bitsOut(2*t-1) = cout(t,1); % Parse cout to create a single stream of bits           
                fsm.bitsOut(2*t) = cout(t,2);
            end
            bitsOut = fsm.bitsOut;
            x = fsm.x;
        end
        
        % BCJR soft-decoder
        function [llrApp,gamma,alpha,beta] = bcjr(fsm,llrPri,llrExt)
            nt = size(fsm.bitsIn_TailBitsAppended,1); % Length of the trellis incl. tail bits
            nout = 2; % number of output bits per step (should be 2 for this case)
            [nx,nin] = size(fsm.outMap);
            gamma = zeros(nx,nin,nt);
            % add the prior probabilities
            for ix = 1:nx
                for ib = 1:nin
                    for t = 1:nt;
                        gamma(ix,ib,t) = (ib-1)*llrPri(t);    % llrPri is a 1*nt vector
                    end
                end
            end
            % add the external probabilities
            for ix = 1:nx
                for ib = 1:nin
                    for t = 1:nt
                        for j = 1:nout
                            gamma(ix,ib,t) = gamma(ix,ib,t) + bitget(fsm.outMap(ix,ib),3-j)*llrExt(j,t);    % llrExt is a 2*nt matrix
                        end
                    end
                end
            end    
            % Assign all the branch metrics to very large negative values
            largeNeg = -100;
            gamma(:,2:nin,nt-fsm.ntail+1:nt) = largeNeg;
            fsm.gamma = gamma;
        
            % Compute the forward metric
            alpha = largeNeg*ones(nx,nt+1); % Initialize all values to large negative
            alpha(1,1) = 0; % values except
            for t = 1:nt
                for ix0 = 1:nx
                    for iin = 1:nin
                        ix = fsm.stateMap(ix0,iin)+1;
                        alpha1 = alpha(ix0,t) + gamma(ix0,iin,t);
                        alpha(ix,t+1) = FSMBinary.maxstar(alpha(ix,t+1), alpha1);
                    end
                end
            end
            fsm.alpha = alpha;
            
            % Compute the reverse metrics
            beta = largeNeg*ones(nx,nt+1);
            beta(:,nt+1) = 0;
            for t = nt:-1:1
                for ix = 1:nx
                    for iin = 1:nin
                        ix1 = fsm.stateMap(ix,iin)+1;
                        beta1 = beta(ix1,t+1) + gamma(ix,iin,t);
                        beta(ix,t) = FSMBinary.maxstar(beta(ix,t), beta1); 
                    end
                end
            end 
            fsm.beta = beta;
            
            % APP likelihood
            llrApp = largeNeg*ones(nt,nin);
            for t = 1:nt
                for iin = 1:nin
                    for ix = 1:nx
                        ix1 = fsm.stateMap(ix,iin)+1;                        
                        llrApp(t,iin) = FSMBinary.maxstar(llrApp(t,iin),alpha(ix,t)+gamma(ix,iin,t)+beta(ix1,t+1));
                    end
                end
            end         
            llrApp = llrApp(:,2)-llrApp(:,1);
            fsm.llrApp = llrApp;
        end     
    end
    
    methods (Static)
        function z = maxstar(a,b)
            z = max(a,b) + log(1 + exp(-abs(b-a)));
        end
    end
    
end