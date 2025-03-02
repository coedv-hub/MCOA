
function o =  Skin_darkening_or_lightening(Xbest,  X, SearchAgents_no)
 darkening= [0.0, 0.4046661];
 lightening= [0.5440510 , 1.0];

    dark1= darkening(1) + (darkening (2) - darkening(1))* rand();
    dark2= darkening(1) + (darkening (2) - darkening(1))* rand();
    light1= lightening(1) + (lightening(2)-lightening(1))* rand();
    light2= lightening(1) + (lightening(2)-lightening(1))* rand();

  [r1, r2, r3, r4]= R(SearchAgents_no);

    if(getBinary)
         o= Xbest +   light1*sin((X(r1,:)-X(r2,:))/2) - ((-1)^getBinary) * light2*sin((X(r3,:)-X(r4,:))/2);
    else
         o= Xbest +    dark1*sin((X(r1,:)-X(r2,:))/2) - ((-1)^getBinary) * dark2*sin((X(r3,:)-X(r4,:))/2);
    end
end

