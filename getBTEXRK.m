% set the Butcher tableau for explicit Runge-Kutta method
function bt = getBTEXRK(order)

switch order
    % Euler forward
    case 1
        bt.nstages = 1;
        bt.order = 1;
        bt.A = 0;
        bt.b = 1;
        bt.c = 0;
        
    % TVD, second order
    case 2
        bt.nstages = 2;
        bt.order = 2;
        bt.A = [0, 0; 1, 0];
        bt.b = [0.5; 0.5];
        bt.c = [0; 1];
    
    % TVD, third order
    case 3
        bt.nstages = 3;
        bt.order = 3;
        bt.A = [0, 0, 0; 1, 0, 0; 0.25, 0.25, 0];
        bt.b = [1 / 6; 1 / 6; 2 / 3];
        bt.c = [0; 1; 0.5];
        
    % classical fourth order
    case 4
        bt.nstages = 4;
        bt.order = 4;
        bt.A = [  0,   0,   0, 0;
                0.5,   0,   0, 0;
                  0, 0.5,   0, 0;
                  0,   0,   1, 0];
        bt.b = [1 / 6; 1 / 3; 1 / 3; 1 / 6];
        bt.c = [0; 0.5; 0.5; 1];
       
    % fifth order
    case 5
        bt.nstages = 6;
        bt.order = 5;
        bt.A = [       0,       0,         0,        0,          0, 0;
                     0.5,       0,         0,        0,          0, 0;
                    0.25,    0.25,         0,        0,          0, 0; 
                       0,      -1,         2,        0,          0, 0;
                  7 / 27, 10 / 27,         0,   1 / 27,          0, 0;
                28 / 625,  -1 / 5, 546 / 625, 54 / 625, -378 / 625, 0];
        bt.b = [1 / 24; 0; 0; 5 / 48; 27 / 56; 125 / 336];
        bt.c = [0; 0.5; 0.5; 1; 2 / 3; 0.2];
        
    otherwise
        error('Not implemented explicit Runge-Kutta time integration method')
end

