function [num] = addzeros(iframe)
%This function adds zeros to iframe

if (iframe<10)

      num = strcat('000',num2str(iframe));

    elseif (iframe<100)
    
      num = strcat('00',num2str(iframe));

    elseif (iframe<1000)
    
      num = strcat('0',num2str(iframe));

    elseif (iframe<10000)
    
      num = strcat(num2str(iframe));

    end %if

end

