function cfl = setCFLNumber(k)

switch k
    case 1
        cfl = 0.3;
    case 2
        cfl = 0.18;
    otherwise
        error('Only polynomial degrees 1 and 2 are supported')
end

end