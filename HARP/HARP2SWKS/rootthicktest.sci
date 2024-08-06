

for i = 1:hubRoot
    figure();
    temp = cell2mat(airfoils(i));AAA=zeros(1,82);
    plot(temp(:,1),temp(:,2),'b+');
    [ccx, ccy] = chordcentreline(temp(:,1),temp(:,2));
    plot(ccx, ccy,'g+');
    Xmid = min(temp(:,1))+(max(temp(:,1))-min(temp(:,1)))/2;
    Xoff = abs(temp(:,1) - Xmid);
    Xoff = Xoff-min(Xoff);
    Xoff = abs(Xoff-max(Xoff));
    Zoff = 0.03*abs(temp(1,3)-minRad);
    
    for j = 1:length(temp(:,1))
        Yctr = ccy(find(ccx(:)==temp(j,1),1));
        if temp(j,2) > Yctr then
            temp(j,2) = temp(j,2) + Zoff*(Xoff(j)*0.12)^.36;
            AAA(j)= Zoff*(Xoff(j)*0.12)^0.37;            
        else
            temp(j,2) = temp(j,2) - Zoff*(Xoff(j)*0.12)^.35;
            AAA(j)= -Zoff*(Xoff(j)*0.12)^0.35;            
        end
    end
    plot(temp(:,1),AAA(:),'r+');
    plot(temp(:,1),temp(:,2),'r');
    
    //break;pause
end
