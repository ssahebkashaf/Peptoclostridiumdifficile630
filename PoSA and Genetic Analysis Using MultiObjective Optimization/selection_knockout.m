%% funzione che mi genera quali geni selezionare

function f = selection_knockout

% posso in modo random sopprimere un numero max di geni.
value= randi(5,1);
sol=false;
f=0;
if value==0
    return
end
while sol==false
    for value1=1:5
        if value == value1 % scelgo il numero di knockout in modo random
            % posso scegliere un solo numero da 1 a 632    
            for j=1:value
                f(j) = randi(632,1);
            end
            sol=true;
            continue
        end
    end
end
