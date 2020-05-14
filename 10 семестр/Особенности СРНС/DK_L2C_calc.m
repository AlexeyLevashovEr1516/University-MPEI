%% Функция расчета дальномерного кода для сигнала GPS L2C
function DK_L2C_out = DK_L2C_calc( ISRS,L )

DK_L2C_out01(1) = ISRS(27);

for k = 2:L
    % Выход
    DK_L2C_out01(k) = ISRS(27);
    
    % Обратная связь
    ISRS(3)  = xor( ISRS(27),ISRS(3) );
    ISRS(6)  = xor( ISRS(27),ISRS(6) );
    ISRS(8)  = xor( ISRS(27),ISRS(8) );
    ISRS(11) = xor( ISRS(27),ISRS(11) );
    ISRS(14) = xor( ISRS(27),ISRS(14) );
    ISRS(16) = xor( ISRS(27),ISRS(16) );
    ISRS(18) = xor( ISRS(27),ISRS(18) );
    ISRS(21) = xor( ISRS(27),ISRS(21) );
    ISRS(22) = xor( ISRS(27),ISRS(22) );
    ISRS(23) = xor( ISRS(27),ISRS(23) );
    ISRS(24) = xor( ISRS(27),ISRS(24) );
    
    % Сдвиг
    ISRS_copy = ISRS;
    ISRS(1) = ISRS_copy(27);
    for i = 2:27
        ISRS(i) = ISRS_copy(i-1);
    end
    clear ISRS_copy
end

% Получение массива +-1
for k = 1:L
    if DK_L2C_out01(k)
        DK_L2C_out(k) = -1;
    else
        DK_L2C_out(k) = +1;
    end
end

end