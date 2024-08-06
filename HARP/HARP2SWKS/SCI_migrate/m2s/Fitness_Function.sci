function [F] = Fitness_Function(PvsV,Torque,pUvars)

// Output variables initialisation (not found in input variables)
F=[];

// Display mode
mode(0);

// Display warning for floating point exception
ieee(1);

//Evaluates fitness function and returns the fitness value
//Inputs: PvsV: a two column vector, first column is the flow speeds, second column is the power values
//        pUvars: vector which contains the appropriate variables needed for the selected probability distribution, pUvars  is defined in Main.m
//        Torque: vector, the torque values
//        
//Output: F: scalar, the fitness value 

global("PitCtrl","SpdCtrl","Prated","OptMethod")

V = PvsV(:,1);//flow speed (m/s)
P = PvsV(:,2);//rotor power (kW)

if mtlb_logic(OptMethod,"==",0) then //Optimize Efficiency

 %v1_2 = %f; if mtlb_logic(SpdCtrl,"==",1) then %v1_2 = mtlb_logic(PitCtrl,"==",0);end;
 if %v03 then //Fixed-Speed Fixed-Pitch
  // !! L.17: Unknown function SimpInt not converted, original calling sequence used.
  F = SimpInt(V,abs(mtlb_s(Prated,P)),1,max(size(V)));
 elseif %v1_2 then //Variable-Speed Fixed-Pitch
  %v2_2 = Torque; TrqMax = mtlb_max(%v2_2,firstnonsingleton(%v2_2)); //Maximum Torque experienced over the range of flow speeds
  // !! L.20: Unknown function SimpInt not converted, original calling sequence used.
  F = mtlb_a(SimpInt(V,abs(mtlb_s(Prated,P)),1,max(size(V))),mtlb_e(TrqMax,1));
 
  //        Area = SimpInt(V,abs(Prated - P),1,length(V));
  //        Vrated = FindInd(P,Prated); %Index of the rated velocity
  //        Trated = Torque(Vrated);
  //        Texcess = Torque - Trated;
  //        Texcess(Texcess<0) = 0;       
  //        Tarea2 = SimpInt(V,Texcess,1,length(V));
  //        
  //        if Vrated < length(V)
  //           Tarea1 = SimpInt(V,Torque,Vrated,length(V));        
  //        else %rated power was never acheived
  //           Tarea1 = 0;
  //        end
  //        F = Area*(1 + Tarea2/(Tarea1 + Tarea2));       
 
 elseif mtlb_logic(PitCtrl,"~=",0) then //Fixed-Speed Variable-Pitch or Variable-Speed Variable-Pitch
  Vrated = FindInd(P,Prated); //Index of the rated velocity
  // !! L.38: Unknown function SimpInt not converted, original calling sequence used.
  F = SimpInt(V,mtlb_s(Prated,P),1,Vrated);
 end;

elseif mtlb_logic(OptMethod,"==",1) then //Optimize Annual Energy Production

 %v4_2 = %f; if mtlb_logic(SpdCtrl,"==",1) then %v4_2 = mtlb_logic(PitCtrl,"==",0);end;
 if %v33 then //Fixed-Speed Fixed-Pitch
  // !! L.44: Unknown function get_AEP not converted, original calling sequence used.
  AEP = get_AEP(PvsV,pUvars);
  Pexcess = mtlb_s(PvsV(:,2),Prated);
  Pexcess = mtlb_i(Pexcess,mtlb_logic(Pexcess,"<",0),0);
  // !! L.47: Unknown function SimpInt not converted, original calling sequence used.
  Area2 = SimpInt(V,Pexcess,1,max(size(V)));
  Vrated = FindInd(P,Prated); //Index of the rated velocity
  if mtlb_logic(Vrated,"<",max(size(V))) then
    // !! L.50: Unknown function SimpInt not converted, original calling sequence used.
    Area1 = SimpInt(V,P,Vrated,max(size(V)));
    Penalty = Area2/mtlb_a(Area1,Area2);
  else
    Area1 = 0;
    Penalty = 0;
  end;
  //        [Area1 Area2 Penalty]
  F = -(1*AEP)*mtlb_s(1,Penalty);
 
 elseif %v4_2 then //Variable-Speed Fixed-Pitch
  // !! L.60: Unknown function get_AEP not converted, original calling sequence used.
  AEP = get_AEP(PvsV,pUvars);
  Pexcess = mtlb_s(PvsV(:,2),Prated);
  Pexcess = mtlb_i(Pexcess,mtlb_logic(Pexcess,"<",0),0);
  // !! L.63: Unknown function SimpInt not converted, original calling sequence used.
  Area2 = SimpInt(V,Pexcess,1,max(size(V)));
  Vrated = FindInd(P,Prated); //Index of the rated velocity
  Trated = mtlb_e(Torque,Vrated);
  Texcess = mtlb_s(Torque,Trated);
  Texcess = mtlb_i(Texcess,mtlb_logic(Texcess,"<",0),0);
  // !! L.68: Unknown function SimpInt not converted, original calling sequence used.
  Tarea2 = SimpInt(V,Texcess,1,max(size(V)));
 
  if mtlb_logic(Vrated,"<",max(size(V))) then
    // !! L.71: Unknown function SimpInt not converted, original calling sequence used.
    Area1 = SimpInt(V,P,Vrated,max(size(V)));
    // !! L.72: Unknown function SimpInt not converted, original calling sequence used.
    Tarea1 = SimpInt(V,Torque,Vrated,max(size(V)));
    Penalty1 = Area2/mtlb_a(Area1,Area2);
    Penalty2 = Tarea2/mtlb_a(Tarea1,Tarea2);
    %v55 = %t;  if ~isnan(Penalty1) then %v55 = isnan(Penalty2);end;
    if %v55 then
      Penalty1 = 0;
      Penalty2 = 1;
    end;
  else //rated power was never acheived
   Area1 = 0;
   Tarea1 = 0;
   Penalty1 = 0;
   Penalty2 = 0;
  end;
  //        [Area1 Area2 Tarea1 Tarea2 Penalty1 Penalty2]
  F = -(1*AEP)*mtlb_s(mtlb_s(1,Penalty1),Penalty2);
 
 elseif mtlb_logic(PitCtrl,"~=",0) then //Fixed-Speed Variable-Pitch or Variable-Speed Variable-Pitch
  Vrated = FindInd(P,Prated); //Index of the rated velocity
  // !! L.90: Unknown function get_AEP not converted, original calling sequence used.
  AEP = get_AEP(PvsV(mtlb_imp(1,Vrated),:),pUvars); //only send the part of the power curve below Vrated
  F = -1*AEP;
 end;

end;
endfunction
