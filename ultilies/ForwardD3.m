function [Dux,Duy] = ForwardD3(U)
            Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
            Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
        end