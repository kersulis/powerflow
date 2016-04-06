function [Vmag,Vang,Pcalc,Qcalc,ty] = powerflow(Vmag,Vang,Pinj,Qinj,ty,Vmax)

global tol maxiter nbus G B Y

converged = 0;

for mm=1:maxiter
    Pcalc = zeros(nbus,1); 
    Qcalc = zeros(nbus,1); 
    
  
   for i=1:nbus
      for k = find(Y(i,:))
         Pcalc(i) = Pcalc(i) + Vmag(i)*Vmag(k)*(G(i,k)*cos(Vang(i)-Vang(k)) + ...
                                                B(i,k)*sin(Vang(i)-Vang(k)));
         Qcalc(i) = Qcalc(i) + Vmag(i)*Vmag(k)*(G(i,k)*sin(Vang(i)-Vang(k)) - ...
                                                B(i,k)*cos(Vang(i)-Vang(k)));
      end
      if ty(i) == 2 && Vmax(i) > 0
         if Qcalc(i) > Qinj(i)
            ty(i) = 1;
         end
      end
   end
                
   dPda = spalloc(nbus,nbus,nnz(Y));
   dPdV = spalloc(nbus,nbus,nnz(Y));
   dQda = spalloc(nbus,nbus,nnz(Y));
   dQdV = spalloc(nbus,nbus,nnz(Y));
   Pmis = zeros(1,nbus);
   Qmis = zeros(1,nbus);
   
   for i=1:nbus
      Vmi = Vmag(i);
      for k = find(Y(i,:))
         if i==k
            dPda(i,i) = -B(i,i)*Vmi*Vmi - Qcalc(i);
            dPdV(i,i) = G(i,i)*Vmi + Pcalc(i)/Vmi;
            dQda(i,i) = -G(i,i)*Vmi*Vmi + Pcalc(i);
            dQdV(i,i) = -B(i,i)*Vmi + Qcalc(i)/Vmi;
         else
            Vmk = Vmag(k);
            Vai = Vang(i);
            Vak = Vang(k);
            dPda(i,k) = Vmi*Vmk*(G(i,k)*sin(Vai-Vak)-B(i,k)*cos(Vai-Vak));
            dPdV(i,k) = Vmi*(G(i,k)*cos(Vai-Vak)+B(i,k)*sin(Vai-Vak));
            dQda(i,k) = -Vmi*Vmk*(G(i,k)*cos(Vai-Vak)+B(i,k)*sin(Vai-Vak));
            dQdV(i,k) = Vmi*(G(i,k)*sin(Vai-Vak)-B(i,k)*cos(Vai-Vak));
         end
      end
   end

   for i=1:nbus
      if ty(i) == 1
         Pmis(i) = Pcalc(i) - Pinj(i);
         Qmis(i) = Qcalc(i) - Qinj(i);
      else                                
         Pmis(i) = Pcalc(i) - Pinj(i);
         Qmis(i) = 0;
         dQda(i,:) = zeros(1,nbus);
         dQdV(i,:) = zeros(1,nbus);
         dPdV(:,i) = zeros(nbus,1); % PV bus: V constant
         dQdV(:,i) = zeros(nbus,1);
         dQdV(i,i) = 1;
      end
      if ty(i) == 3
         Pmis(i) = 0;
         dPda(i,:) = zeros(1,nbus);
         dPdV(i,:) = zeros(1,nbus);
         dPda(:,i) = zeros(nbus,1);
         dQda(:,i) = zeros(nbus,1);
         dPda(i,i) = 1;
      end
   end

   mis = [Pmis Qmis]';                
   %keyboard
   if max(abs(mis))<tol
      converged = 1;
      %disp('Reached convergence')
      break
   else                                  
      J = [dPda dPdV;dQda dQdV];              
      update = -J\mis;                    %full operator and result (not sparse)    

      Vang = Vang + update(1:nbus);
      Vmag = Vmag + update(nbus+1:2*nbus);
      
      for i=1:nbus
         if ty(i) == 1 && Vmax(i) > 0
            if Vmag(i) > Vmax(i)
               Vmag(i) = Vmax(i);
               ty(i) = 2;
               disp(['Bus ' num2str(i) ' has hit voltage limit when Q=' ...
                     num2str(Qinj(i)*10)])
            end
         end 
      end
   end
end

if ~converged
   disp('Powerflow not converged')
end

% Jacobian gives sensitivities. Compare to DC power flow.