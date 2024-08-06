function [F] = Fitness_Function(PvsV,Torque,pUvars)
// Evaluates fitness function and returns the fitness value
// Inputs: PvsV: a two column vector, first column is the flow speeds, second column is the power values
//         pUvars: vector which contains the appropriate variables needed for the selected probability distribution, pUvars  is defined in Main.m
//         Torque: vector, the torque values
//         
// Output: F: scalar, the fitness value 

global PitCtrl SpdCtrl Prated OptMethod

V = PvsV(:,1);  // flow speed (m/s)
P = PvsV(:,2);  // rotor power (kW)

if OptMethod == 0    // Optimize Efficiency

    if SpdCtrl == 0 & PitCtrl == 0;      // Fixed-Speed Fixed-Pitch
       F = SimpInt(V,abs(Prated - P),1,length(V));  
    elseif SpdCtrl == 1 & PitCtrl == 0;  // Variable-Speed Fixed-Pitch
       TrqMax = max(Torque); // Maximum Torque experienced over the range of flow speeds
       F = SimpInt(V,abs(Prated - P),1,length(V)) + TrqMax(1);        
       
//         Area = SimpInt(V,abs(Prated - P),1,length(V));
//         Vrated = FindInd(P,Prated); // Index of the rated velocity
//         Trated = Torque(Vrated);
//         Texcess = Torque - Trated;
//         Texcess(Texcess<0) = 0;       
//         Tarea2 = SimpInt(V,Texcess,1,length(V));
//         
//         if Vrated < length(V)
//            Tarea1 = SimpInt(V,Torque,Vrated,length(V));        
//         else // rated power was never acheived
//            Tarea1 = 0;
//         end
//         F = Area*(1 + Tarea2/(Tarea1 + Tarea2));       
       
    elseif PitCtrl ~= 0;                  // Fixed-Speed Variable-Pitch or Variable-Speed Variable-Pitch
       Vrated = FindInd(P,Prated); // Index of the rated velocity
       F = SimpInt(V,Prated-P,1,Vrated);    
    end

elseif OptMethod == 1 // Optimize Annual Energy Production
    
    if SpdCtrl == 0 & PitCtrl == 0;      // Fixed-Speed Fixed-Pitch
       AEP = get_AEP(PvsV,pUvars);
       Pexcess = PvsV(:,2) - Prated;
       Pexcess(Pexcess<0) = 0;
       Area2 = SimpInt(V,Pexcess,1,length(V));
       Vrated = FindInd(P,Prated); // Index of the rated velocity
       if Vrated < length(V)
          Area1 = SimpInt(V,P,Vrated,length(V));
          Penalty = Area2/(Area1 + Area2);
       else
          Area1 = 0; 
          Penalty = 0;
       end
//         [Area1 Area2 Penalty]
       F = -1*AEP*(1 - Penalty);
       
    elseif SpdCtrl == 1 & PitCtrl == 0;  // Variable-Speed Fixed-Pitch
       AEP = get_AEP(PvsV,pUvars);
       Pexcess = PvsV(:,2) - Prated;
       Pexcess(Pexcess<0) = 0;
       Area2 = SimpInt(V,Pexcess,1,length(V));
       Vrated = FindInd(P,Prated); // Index of the rated velocity
       Trated = Torque(Vrated);
       Texcess = Torque - Trated;
       Texcess(Texcess<0) = 0;       
       Tarea2 = SimpInt(V,Texcess,1,length(V));
       
       if Vrated < length(V)
          Area1 = SimpInt(V,P,Vrated,length(V));
          Tarea1 = SimpInt(V,Torque,Vrated,length(V));
          Penalty1 = Area2/(Area1 + Area2);
          Penalty2 = Tarea2/(Tarea1 + Tarea2);
          if isnan(Penalty1) | isnan(Penalty2)
             Penalty1 = 0;
             Penalty2 = 1;
          end
       else // rated power was never acheived
          Area1 = 0;
          Tarea1 = 0;
          Penalty1 = 0;
          Penalty2 = 0;
       end
//         [Area1 Area2 Tarea1 Tarea2 Penalty1 Penalty2]
       F = -1*AEP*(1 - Penalty1 - Penalty2);       
       
    elseif PitCtrl ~= 0;                  // Fixed-Speed Variable-Pitch or Variable-Speed Variable-Pitch
       Vrated = FindInd(P,Prated); // Index of the rated velocity
       AEP = get_AEP(PvsV(1:Vrated,:),pUvars); // only send the part of the power curve below Vrated
       F = -1*AEP;   
    end
    
end
endfunction
