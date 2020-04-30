classdef MyPerceptron < handle
    properties
        neuronNum;
        connVec;
    end
    
    methods
        
        function obj = MyPerceptron(neuronNum, connVec)
            if(nargin < 2)
                obj.connVec = randn(1,neuronNum+1);
            else
                % extra connection for constant displacement term
                obj.connVec = [connVec, mean(connVec)];
            end
            obj.neuronNum = neuronNum;
        end
        
        function success = Train(obj, s, sigma, eta, epochNumMax, verbose)
            
            if(nargin < 6)
                verbose = 0;
            end
            
            % code by: Shaul Druckmann
            errorNum = 1;
            epochNum = 1;
            
            % constant displacement term
            s = [s; ones(1, size(s,2))];

            while (errorNum > 0 && epochNum < epochNumMax )
              errorNum = 0;
              for pp=1:size(s,2)
                %  Single trial
                outputAct = obj.connVec*s(:,pp);
                if sign(outputAct) ~= sign(sigma(pp)) % Mistake made
                  obj.connVec = obj.connVec + 2*eta*(s(:,pp)*sign(sigma(pp)))';
                  errorNum = errorNum + 1;
                end
            %     neuronConnectivityMat(inputNeuronNum,1:(inputNeuronNum+1)) = connVec;
              end
              if (verbose)
                  disp(['In epoch ' num2str(epochNum) ' There were ' ...
                    num2str(errorNum) ' mistakes out of ' num2str(size(s,2)) ' patterns']);
              end
              
              epochNum = epochNum + 1;
            end
            
            success = (errorNum == 0);
        end
    end
    
end