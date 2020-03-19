function [cm acc fm ] = confusionMatrix( y_actual,y_predicted )
tp=0;tn=0;fp=0;fn=0;
for i=1:length(y_actual)
if y_actual(i)>0
    if y_actual(i)==y_predicted(i)
        tp=tp+1;
    else
        fn=fn+1;
    end
end
if y_actual(i)<0
    if y_actual(i)==y_predicted(i)
        tn=tn+1;
    else
        fp=fp+1;
    end
end
end
cm=[tn fp; fn tp];
acc=(tp+tn)/(tp+fn+fp+tn);
sens=tp/(tp+fn);prec=tp/(tp+fp);
fm=(2*prec*sens)/(prec+sens);
end

