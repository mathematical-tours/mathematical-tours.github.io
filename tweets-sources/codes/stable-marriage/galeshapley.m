function stablematch = galeshapley(N, men_pref, women_pref)

%---------GALE - SHAPLEY ALGORITHM (Men propose) --------%

men_free = zeros(N,1);
women_suitor = zeros(N,N);
women_partner = zeros(N,1);
rank = zeros(N,N);


for i = 1:N
    for j = 1:N
        for k = 1:N
        if(women_pref(i,k) == j)
            rank(i,j) = k;
        end
        end
    end
end

while (min(women_partner) == 0)
    for i = 1:N
        if (men_free(i,1) == 0)
            next = find(men_pref(i,:) > 0, 1);
            women_suitor(men_pref(i,next),i) = i;
            men_pref(i,next) = 0;
        end
    end
    for i = 1:N
        for j = 1:N
            if(women_suitor(i,j) ~= 0)
                if(women_partner(i,1) == 0)
                    women_partner(i,1) = women_suitor(i,j);
                    men_free(j,1) = 1;
                end
                if(women_partner(i,1) ~= 0)
                if(rank(i,women_suitor(i,j)) < rank(i,women_partner(i,1)))
                    men_free(women_partner(i,1),1) = 0;
                    women_partner(i,1) = women_suitor(i,j);
                    men_free(j,1) = 1;
                    
                end
                end
            end
        end
    end
end

stablematch = women_partner;



          
        
        
        
        